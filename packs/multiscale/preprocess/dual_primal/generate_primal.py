from networkx.convert_matrix import _generate_weighted_edges
import numpy as np
from packs.utils import utils_old
import networkx as nx
import collections
from packs.multiscale.ms_utils.multiscale_functions import print_mesh_volumes_data

def create_dual_and_primal(centroids_volumes: np.ndarray, volumes_dimension: np.ndarray, centroids_nodes: np.ndarray, cr: np.ndarray, adjacencies_internal_faces: np.ndarray, **kwargs):
    # import pdb; pdb.set_trace()
    M = kwargs.get('m_object')
    
    ##########################
    ## around the vector in order to fix the gmsh centroid location problem
    centroids_volumes2 = np.around(centroids_volumes, decimals=8)
    ## comment this block for real cases
    ######################
    
    # ###########
    # #### uncomment this block for real cases
    # centroids_volumes2 = centroids_volumes
    # ###########
    
    dimension = get_dimension_of_problem(centroids_volumes2)
    L = get_l_total(centroids_nodes)
    primal_coarse_ids, all_separated = create_primal(centroids_volumes2, volumes_dimension, cr, dimension)
    dual_ids = create_dual(centroids_volumes2, cr, dimension, L, volumes_dimension, adjacencies_internal_faces, all_separated)
    # M.data['DUAL_1'][:] = dual_ids
    # M.data['GID_1'][:] = primal_coarse_ids
    # M.data.update_variables_to_mesh()
    # print_mesh_volumes_data(M, 'results/dual_test.vtk')
    
    # import pdb; pdb.set_trace()
    test_coarse_vertices(primal_coarse_ids, dual_ids)
    
    return primal_coarse_ids, dual_ids

def get_l_total(centroids_nodes):
    min_nodes, max_nodes = get_min_and_max(centroids_nodes)
    L = max_nodes - min_nodes
    return L

def create_primal(centroids_volumes: np.ndarray, volumes_dimension: np.ndarray , cr: np.ndarray, dimension: int):
    
    all_separated = get_all_separated(dimension, centroids_volumes, cr)
    all_separated_limits = get_all_separated_limits(all_separated, centroids_volumes)
    primals = get_coarse_volumes(all_separated_limits, centroids_volumes, volumes_dimension)
    return primals, all_separated

def get_dimension_volumes(volumes_nodes, centroids_nodes):
    centroids = centroids_nodes[volumes_nodes]
    mins = centroids.min(axis=1)
    maxs = centroids.max(axis=1)
    dimensions = maxs - mins   
    return dimensions    

def get_min_and_max(vector: np.ndarray):
    
    mins = vector.min(axis=0)
    maxs = vector.max(axis=0)
    
    return mins, maxs
    
def get_dimension_of_problem(centroids_volumes: np.ndarray):
    dimension = 0
    
    unique_cents = ([
        np.unique(centroids_volumes[:, 0]),
        np.unique(centroids_volumes[:, 1]),
        np.unique(centroids_volumes[:, 2])
    ])
    
    for cent in unique_cents:
        if len(cent) > 1:
            dimension += 1
    
    return dimension

def get_all_separated(dimension: int, centroids_volumes: np.ndarray, cr: np.ndarray):
    
    all_separated = []
    for dim in range(dimension):
        unique_centroids_volumes = np.unique(centroids_volumes[:, dim])
        n_unique = len(unique_centroids_volumes)
        n_coarse_volumes = n_unique//cr[dim]
        resto = n_unique % cr[dim]
        n_fine_total = n_coarse_volumes*cr[dim]
        separated = list(unique_centroids_volumes[0:n_fine_total].reshape((n_coarse_volumes, cr[dim])))
        if resto > 0:
            separated[-1] = np.concatenate([separated[-1], unique_centroids_volumes[n_fine_total:n_unique]])
        all_separated.append(np.array(separated))
    
    n = len(all_separated)
    if n < 3:
        for i in range(3 - n):
            all_separated.append(np.array([]))
            
    return np.array(all_separated)

def get_all_separated_limits(all_separated: np.ndarray, centroids_volumes: np.ndarray):
    all_separated_limits = []
    
    for i, sep in enumerate(all_separated):
        if len(sep) > 0:
            lim_min = []
            lim_max = []
            for separated_dimension in sep:
                lim_min.append(separated_dimension.min())
                lim_max.append(separated_dimension.max())
            lim_min = np.array(lim_min)
            lim_max = np.array(lim_max)
            gg = np.vstack([lim_min, lim_max]).T
            all_separated_limits.append(gg)
        else:
            all_separated_limits.append(np.array([np.array([centroids_volumes[:,i].min(), centroids_volumes[:,i].max()])]))
    
    return np.array(all_separated_limits)

def get_coarse_volumes(all_separated_limits: np.ndarray, centroids_volumes: np.ndarray, volumes_dimension: np.ndarray):
    dh = volumes_dimension.min()/4
    
    coarse_ids = np.repeat(-1, len(centroids_volumes))
    coarse_id = 0

    for limits0 in all_separated_limits[0]:
        
        for limits1 in all_separated_limits[1]:
                
            for limits2 in all_separated_limits[2]:
                
                limites = np.array([
                    np.array([limits0[0]-dh, limits1[0]-dh, limits2[0]-dh]),
                    np.array([limits0[1]+dh, limits1[1]+dh, limits2[1]+dh])                    
                ])
                
                indexes = utils_old.get_box(centroids_volumes, limites)
                coarse_ids[indexes] = coarse_id
                coarse_id += 1
    
    return coarse_ids

def create_dual(centroids_volumes: np.ndarray, cr: np.ndarray, dimension: int, l_total: np.ndarray, volumes_dimension: np.ndarray, adjacencies_internal_faces: np.ndarray, all_separated: np.ndarray):
    
    vertices_ids = define_vertices(centroids_volumes, cr, l_total, volumes_dimension, all_separated)
    # identify_vertices_id = np.repeat(0, len(vertices_ids))
    # identify_vertices(vertices_ids, centroids_volumes, identify_vertices_id)
    
    dual_ids = np.repeat(0, len(centroids_volumes))
    # dual_ids0 = np.repeat(0, len(centroids_volumes))
    
    # g = nx.Graph()
    # g.add_edges_from(adjacencies_internal_faces)
    
    # create_dual_entities(vertices_ids, dual_ids0, centroids_volumes, dimension, volumes_dimension, g, identify_vertices_id)
    create_dual_entities2(vertices_ids, dual_ids, centroids_volumes, volumes_dimension, l_total)
    
    # print(np.allclose(dual_ids, dual_ids0))
    # import pdb; pdb.set_trace()
    return dual_ids

def define_vertices(centroids_volumes: np.ndarray, cr: np.ndarray, l_total: np.ndarray, volumes_dimension, all_separated):
    # gids = np.arange(len(centroids_volumes))
    centroids_min = centroids_volumes.min(axis=0)
    centroids_max = centroids_volumes.max(axis=0)
    vertices_by_direction = []
    for dim in range(3):
        vi = define_vertices_by_direction(centroids_min[dim], centroids_max[dim], cr[dim], l_total[dim], centroids_volumes[:, dim], volumes_dimension, all_separated[dim])
        vertices_by_direction.append(vi)
    
    vertices_by_direction = np.array(vertices_by_direction)
    vv = vertices_by_direction
    all_vertices = np.array(np.meshgrid(vv[0], vv[1], vv[2]))
    all_vertices = all_vertices.reshape((all_vertices.shape[0], all_vertices.shape[1], all_vertices.shape[2])).T
    all_vertices = all_vertices.reshape((all_vertices.shape[0]*all_vertices.shape[1], all_vertices.shape[2]))
    indexes = []
    
    for vertice in all_vertices:
        dist = np.linalg.norm(centroids_volumes - vertice, axis=1)
        index = np.argwhere(dist <= dist.min())[0][0]
        indexes.append(index)
    
    return np.array(indexes)
    
def define_vertices_by_direction(vmin: float, vmax: float, cr_direction: int, l_dimension: float, centroids_volumes_dimension: np.ndarray, volumes_dimension, all_separated_dimension: np.ndarray):
    # import pdb; pdb.set_trace()
    # min_dimensions_by_axis = volumes_dimension.min(axis=0)
    if vmin == vmax:
        return np.array([vmin])
    else:
        unique_centroids_volumes = np.unique(centroids_volumes_dimension)
        # n_blocks = len(unique_centroids_volumes)
        # # n_coarse_blocks = n_blocks//cr_direction
        # n_coarse_blocks = len(all_separated_dimension)
        # # l_coarse = l_dimension/n_coarse_blocks
        # # coarse_ls = [[l_coarse*i, l_coarse*(i+1)] for i in range(n_coarse_blocks-1)]
        # # coarse_ls.append([coarse_ls[-1][1], l_dimension])
        # # coarse_ls = np.array(coarse_ls)
        # # mean = np.mean(coarse_ls, axis=1)
        mean = []
        for cents in all_separated_dimension:
            mean.append(np.mean(cents))
        mean = np.array(mean)
        vertices = [vmin]
        for j in mean[1:-1]:
            dist = unique_centroids_volumes - j
            abs_dist = np.absolute(dist)
            min_dis = abs_dist.min()
            index = np.argwhere(abs_dist == min_dis)[0]
            index = index[0]
            # if len(index) > 1:
            #     ## mais de um volume proximo do centro
            #     index=index[0]
            # else:
            #     index = index[0]
            vertices.append(unique_centroids_volumes[index])
        vertices.append(unique_centroids_volumes[-1])
        vertices = np.array(vertices)
        
        return vertices

def identify_vertices(vertices_ids, centroids_volumes, identify_vertices_id):
    centroids_vertices = centroids_volumes[vertices_ids]
    cents_min = centroids_vertices.min(axis=0)
    cents_max = centroids_vertices.max(axis=0)
    
    # # all_vertices = np.array(np.meshgrid(vv[0], vv[1], vv[2]))
    # all_vertices = np.array(
    #     np.meshgrid(
    #         np.unique([cents_min[0], cents_max[0]]),
    #         np.unique([cents_min[1], cents_max[1]]),
    #         np.unique([cents_min[2], cents_max[2]])
    #     )
    # )
    
    # all_vertices = all_vertices.reshape((all_vertices.shape[0], all_vertices.shape[1], all_vertices.shape[2])).T
    # all_vertices = all_vertices.reshape((all_vertices.shape[0]*all_vertices.shape[1], all_vertices.shape[2]))
    
    ######
    # identify faces
    faces = []
    for i in [cents_min[0], cents_max[0]]:
        
        indexes = np.argwhere(centroids_vertices[:, 0] == i).flatten()
        faces.append(set(indexes))
        
    for i in [cents_min[1], cents_max[1]]:
        indexes = np.argwhere(centroids_vertices[:, 1] == i).flatten()
        faces.append(set(indexes))
            
    for i in [cents_min[2], cents_max[2]]:
        indexes = np.argwhere(centroids_vertices[:, 2] == i).flatten()
        faces.append(set(indexes))
    
    edges = []
    nf = len(faces)        
    for i in range(nf):
        face1 = faces[i]
        for face2 in faces[i+1:nf]:
            edge = face1 & face2
            if len(edge) > 1:
                edges.append(edge)
    
    edges = np.unique(edges)
    
    nedges = len(edges)
    
    verts = []
    for i in range(nedges):
        edge1 = edges[i]
        for edge2 in edges[i+i: nedges]:
            vert = edge1 & edge2
            if len(vert) == 1:
                verts.append(vert)
    
    verts = np.unique(verts)
    
    faces = np.array([np.array(list(f)) for f in faces]).flatten()
    edges = np.array([np.array(list(f)) for f in edges]).flatten()
    verts = np.array([np.array(list(f)) for f in verts]).flatten()
    
    for f in faces:
        identify_vertices_id[f] = 1
    
    for f in edges:
        identify_vertices_id[f] = 2
    
    for f in verts:
        identify_vertices_id[f] = 3

def create_dual_entities(vertices_ids: np.ndarray, dual_ids: np.ndarray, centroids_volumes: np.ndarray, dimension: int, volumes_dimension: np.ndarray, g: nx.Graph, identify_vertices_ids):
    
    nvertices = pow(2, dimension)
    edges_ids = get_edges_by_vertices(vertices_ids, centroids_volumes, nvertices, g, identify_vertices_ids, dimension)
    faces_ids = get_faces(vertices_ids, centroids_volumes, volumes_dimension)
    
    dual_ids[faces_ids] = 1
    dual_ids[edges_ids] = 2
    dual_ids[vertices_ids] = 3
  
def create_dual_entities2(vertices_ids: np.ndarray, dual_ids: np.ndarray, centroids_volumes: np.ndarray, volumes_dimension: np.ndarray, l_total):
    centroids_vertices = centroids_volumes[vertices_ids]
    edges = []
    faces = [[], [], []]
    delta = volumes_dimension.min()/4
    mins = centroids_vertices.min(axis=0)
    maxs = centroids_vertices.max(axis=0)
    unique_centroids_vertices = np.array([
        np.unique(centroids_vertices[:, 0]),
        np.unique(centroids_vertices[:, 1]),
        np.unique(centroids_vertices[:, 2])
    ])
    volumes_ids = np.arange(len(centroids_volumes))
    
    
    for i, centroids_vertices_direction in enumerate(unique_centroids_vertices):
        # get_faces_entities(centroids_vertices_direction, i, mins, maxs, delta, faces, centroids_volumes, l_total)
        get_faces_entities_v2(centroids_vertices_direction, i, mins, maxs, delta, faces, centroids_volumes, l_total, volumes_ids)
    
    faces = np.array(faces)
    get_edges_entities(faces, edges)
    edges = np.unique(np.concatenate(edges))
    faces = np.unique(np.concatenate([
        np.concatenate(faces[0]),
        np.concatenate(faces[1]),
        np.concatenate(faces[2])
    ]))
    
    dual_ids[faces] = 1
    dual_ids[edges] = 2
    dual_ids[vertices_ids] = 3
    
def get_faces_entities(centroids_vertices_direction, direction, mins, maxs, delta, faces, centroids_volumes, l_total):
    directions = np.array([0, 1, 2])
    directions2 = np.setdiff1d(directions, direction)
    limites = np.zeros((2, 3))
    
    for k in centroids_vertices_direction:
        limites[:] = 0
        limites[0, direction] = k - delta
        limites[1, direction] = k + delta
        limites[0, directions2[0]] = mins[directions2[0]] - delta
        limites[1, directions2[0]] = maxs[directions2[0]] + delta
        limites[0, directions2[1]] = mins[directions2[1]] - delta
        limites[1, directions2[1]] = maxs[directions2[1]] + delta        
        indexes = utils_old.get_box(centroids_volumes, limites)
        faces[direction].append(indexes)
    
    faces[direction] = np.array(faces[direction])

def get_faces_entities_v2(centroids_vertices_direction, direction, mins, maxs, delta, faces, centroids_volumes, l_total, volumes_ids):
    
    for k in centroids_vertices_direction:
        indexes = (centroids_volumes[:, direction] < k + delta) & (centroids_volumes[:, direction] > k - delta)
        faces[direction].append(volumes_ids[indexes])
    
    faces[direction] = np.array(faces[direction])
        
def get_edges_entities(faces, edges):
    directions = np.array([0, 1, 2])
    
    for direction in range(3):
        directions2 = np.setdiff1d(directions, direction)
        faces2 = faces[directions2]
        faces2 = np.unique(np.concatenate([np.concatenate(faces2[0]), np.concatenate(faces2[1])]))
        conj_edges = np.intersect1d(np.concatenate(faces[direction]), faces2)
        if conj_edges.shape[0] > 1:
            edges.append(conj_edges)
        
    
def get_edges_by_vertices(vertices_ids: np.ndarray, centroids_volumes: np.ndarray, nvertices: int, g: nx.Graph, identify_vertices_ids, dimension: int):
    
    centroids_vertices = centroids_volumes[vertices_ids]
    vertices_graph = nx.Graph()
    
    ########
    ## dict for toget (dimension, id_vert)
    
    # dict_toget = {
    #     (1, 3): nvertices-1, # 1
    #     (1, 2): nvertices, # 2
    #     (2, 3): nvertices - 2, # 2
    #     (2, 2): nvertices - 1, # 3
    #     (2, 1): nvertices, # 4
    #     (3, 3): nvertices - 5, # 3
    #     (3, 2): nvertices - 4, # 4
    #     (3, 1): nvertices - 3, # 5
    #     (3, 0): nvertices - 2, # 6
    # }
    
    dict_toget = {
        (1, 3): 1, # 1
        (1, 2): 2, # 2
        (2, 3): 2, # 2
        (2, 2): 3, # 3
        (2, 1): 4, # 4
        (3, 3): 3, # 3
        (3, 2): 4, # 4
        (3, 1): 5, # 5
        (3, 0): 6, # 6
    }    
    ########
    
    all_edges = []
    
    for i, centv in enumerate(centroids_vertices):
        id_vert = identify_vertices_ids[i]
        dist = np.linalg.norm(centroids_vertices - centv, axis=1)
        toget = dict_toget[(dimension, id_vert)] + 1
        dist_max = np.sort(dist)[0:toget].max()
        indexes = np.argwhere(dist <= dist_max).flatten()
        vertice_id = vertices_ids[i]
        other_vertices = np.setdiff1d(vertices_ids[indexes], vertice_id)
        for j in other_vertices:
            if vertices_graph.has_edge(vertice_id, j):
                continue
            path = np.array(nx.shortest_path(g, source=vertice_id, target=j))
            # edges = np.setdiff1d(path, [vertice_id, j])
            vertices_graph.add_edge(vertice_id, j)
            all_edges.append(path)
    
    all_edges = np.unique(np.concatenate(all_edges))
    return all_edges
            
def get_faces(vertices_ids, centroids_volumes, volumes_dimension):
    
    centroids_vertices = centroids_volumes[vertices_ids]
    
    delta = volumes_dimension.min()/4
    
    cents = np.array([
        np.unique(centroids_vertices[:,0]),
        np.unique(centroids_vertices[:,1]),
        np.unique(centroids_vertices[:,2])    
    ])
    
    faces = []
    
    for dim in range(3):
        for x in cents[dim]:
            indexes = np.argwhere(((centroids_volumes[:,dim] < x+delta) & (centroids_volumes[:,dim] > x-delta))).flatten()
            faces.append(indexes)
    
    faces = np.unique(np.concatenate(faces))
    
    return faces

def qualquer():
    pass

def test_coarse_vertices(primal_ids, dual_ids):
    # for cid in np.unique(primal_ids):
    #     local_dual_ids = dual_ids[primal_ids == cid]
    #     counts = collections.Counter(local_dual_ids)

    #     if counts[3] > 1:
    #         print('mais que um vertice numa dual')
    #         raise ValueError
    n_coarse_volumes = len(np.unique(primal_ids))
    n_vertices = len(dual_ids[dual_ids == 3])
    if n_coarse_volumes != n_vertices:
        raise ValueError('Numero de vertices diferente do numero de coarse volumes')
            
        