import numpy as np
from packs.utils import utils_old
import networkx as nx

def create_dual_and_primal(centroids_volumes: np.ndarray, volumes_dimension: np.ndarray, centroids_nodes: np.ndarray, cr: np.ndarray, adjacencies_internal_faces: np.ndarray):
    dimension = get_dimension_of_problem(centroids_volumes)
    L = get_l_total(centroids_nodes)
    primal_coarse_ids = create_primal(centroids_volumes, volumes_dimension, cr, dimension)
    dual_ids = create_dual(centroids_volumes, cr, dimension, L, volumes_dimension, adjacencies_internal_faces)
    
    return primal_coarse_ids, dual_ids

def get_l_total(centroids_nodes):
    min_nodes, max_nodes = get_min_and_max(centroids_nodes)
    L = max_nodes - min_nodes
    return L

def create_primal(centroids_volumes: np.ndarray, volumes_dimension: np.ndarray , cr: np.ndarray, dimension: int):
    
    all_separated = get_all_separated(dimension, centroids_volumes, cr)
    all_separated_limits = get_all_separated_limits(all_separated, centroids_volumes)
    primals = get_coarse_volumes(all_separated_limits, centroids_volumes, volumes_dimension)
    return primals

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
        separated = unique_centroids_volumes[0:n_fine_total].reshape((n_coarse_volumes, cr[dim]))
        if resto > 0:
            import pdb; pdb.set_trace()
            separated[-1] = np.concatenate([separated[-1], unique_centroids_volumes[n_fine_total:n_unique]])
        all_separated.append(separated)
    
    n = len(all_separated)
    if n < 3:
        for i in range(3 - n):
            all_separated.append(np.array([]))
    
    return np.array(all_separated)

def get_all_separated_limits(all_separated: np.ndarray, centroids_volumes: np.ndarray):
    all_separated_limits = []
    
    for i, sep in enumerate(all_separated):
        if len(sep) > 0:
            lim_min = sep.min(axis=1)
            lim_max = sep.max(axis=1)
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

def create_dual(centroids_volumes: np.ndarray, cr: np.ndarray, dimension: int, l_total: np.ndarray, volumes_dimension: np.ndarray, adjacencies_internal_faces: np.ndarray):
    vertices_ids = define_vertices(centroids_volumes, cr, l_total)
    identify_vertices_id = np.repeat(0, len(vertices_ids))
    identify_vertices(vertices_ids, centroids_volumes, identify_vertices_id)
    
    dual_ids = np.repeat(0, len(centroids_volumes))
    
    g = nx.Graph()
    g.add_edges_from(adjacencies_internal_faces)
    
    create_dual_entities(vertices_ids, dual_ids, centroids_volumes, dimension, volumes_dimension, g, identify_vertices_id)
    
    return dual_ids

def define_vertices(centroids_volumes: np.ndarray, cr: np.ndarray, l_total: np.ndarray):
    gids = np.arange(len(centroids_volumes))
    centroids_min = centroids_volumes.min(axis=0)
    centroids_max = centroids_volumes.max(axis=0)
    vertices_by_direction = []
    for dim in range(3):
        vi = define_vertices_by_direction(centroids_min[dim], centroids_max[dim], cr[dim], l_total[dim], centroids_volumes[:, dim])
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
    
def define_vertices_by_direction(vmin: float, vmax: float, cr_direction: int, l_dimension: float, centroids_volumes_dimension: np.ndarray):
    
    if vmin == vmax:
        return np.array([vmin])
    else:
        unique_centroids_volumes = np.unique(centroids_volumes_dimension)
        n_blocks = len(unique_centroids_volumes)
        n_coarse_blocks = n_blocks//cr_direction
        l_coarse = l_dimension/n_coarse_blocks
        coarse_ls = [[l_coarse*i, l_coarse*(i+1)] for i in range(n_coarse_blocks-1)]
        coarse_ls.append([coarse_ls[-1][1], l_dimension])
        coarse_ls = np.array(coarse_ls)
        mean = np.mean(coarse_ls, axis=1)
        vertices = [vmin]
        for j in mean[1:-1]:
            dist = unique_centroids_volumes - j
            abs_dist = np.absolute(dist)
            min_dis = abs_dist.min()
            index = np.argwhere(abs_dist == min_dis)[0][0]
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
