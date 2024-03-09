from pymoab import core, types, rng, topo_util
import numpy as np
from scipy.sparse.csgraph import shortest_path
from scipy import sparse as sp


def get_box_dep0(all_centroids, limites):
    '''
    all_centroids->coordenadas dos centroides do conjunto
    limites-> diagonal que define os volumes objetivo (numpy array com duas coordenadas)
    Retorna os indices cujo centroide está dentro de limites
    '''

    inds0 = np.where(all_centroids[:, 0] > limites[0, 0])[0]
    inds1 = np.where(all_centroids[:, 1] > limites[0, 1])[0]
    inds2 = np.where(all_centroids[:, 2] > limites[0, 2])[0]
    c1 = set(inds0) & set(inds1) & set(inds2)
    inds0 = np.where(all_centroids[:, 0] < limites[1, 0])[0]
    inds1 = np.where(all_centroids[:, 1] < limites[1, 1])[0]
    inds2 = np.where(all_centroids[:, 2] < limites[1, 2])[0]
    c2 = set(inds0) & set(inds1) & set(inds2)
    inds_vols = list(c1 & c2)
    return inds_vols


def get_box(all_centroids, limites):
    ids = np.arange(len(all_centroids))

    '''
    all_centroids->coordenadas dos centroides do conjunto
    limites-> diagonal que define os volumes objetivo (numpy array com duas coordenadas)
    Retorna os indices cujo centroide está dentro de limites
    '''

    inds_vols = ids[
        (all_centroids[:, 0] > limites[0, 0]) & (all_centroids[:, 1] > limites[0, 1]) &
        (all_centroids[:, 2] > limites[0, 2]) & (all_centroids[:, 0] < limites[1, 0]) &
        (all_centroids[:, 1] < limites[1, 1]) & (all_centroids[:, 2] < limites[1, 2])]
    return inds_vols

def get_box_v2(all_centroids, limites):
    inds_vols = \
        (all_centroids[:, 0] > limites[0, 0]) & (all_centroids[:, 1] > limites[0, 1]) & \
        (all_centroids[:, 2] > limites[0, 2]) & (all_centroids[:, 0] < limites[1, 0]) & \
        (all_centroids[:, 1] < limites[1, 1]) & (all_centroids[:, 2] < limites[1, 2])
    return inds_vols

def get_box_v2_2d(all_centroids, limites):
    inds_vols = \
        (all_centroids[:, 0] > limites[0, 0]) & (all_centroids[:, 1] > limites[0, 1]) & \
        (all_centroids[:, 0] < limites[1, 0]) & (all_centroids[:, 1] < limites[1, 1])
    return inds_vols

def getting_tag(mb, name, n, t1, t2, create, entitie, tipo, tags, tags_to_infos):
    types_data = ['handle', 'integer', 'array', 'double']
    entities = ['nodes', 'edges', 'faces', 'volumes', 'root_set', 'intern_faces', 'boundary_faces', 'vols_viz_face',
                'coarse_volumes_lv1', 'coarse_volumes_lv2', 'coarse_volumes']

    assert tipo in types_data, f'tipo nao listado: {tipo}'

    if entitie not in entities:
        raise NameError(f'\nA entidade {entitie} nao esta na lista\n')

    tag = mb.tag_get_handle(name, n, t1, t2, create)
    tag_to_infos = dict(zip(['entitie', 'type', 'n'], [entitie, tipo, n]))
    tags[name] = tag
    tags_to_infos[name] = tag_to_infos


def Min_Max(e, M1):
    verts = M1.mb.get_connectivity(e)
    coords = M1.mb.get_coords(verts).reshape(len(verts), 3)
    xmax, ymax, zmax = coords.max(axis=0)
    xmin, ymin, zmin = coords.min(axis=0)
    return ([xmin, xmax, ymin, ymax, zmin, zmax])


def add_topology(conj_vols,tag_local,lista, mb, mtu, ID_reordenado_tag):
    all_fac=np.uint64(mtu.get_bridge_adjacencies(conj_vols, 2 ,2))
    all_int_fac=np.uint64([face for face in all_fac if len(mb.get_adjacencies(face, 3))==2])
    adjs=np.array([mb.get_adjacencies(face, 3) for face in all_int_fac])
    adjs1=mb.tag_get_data(tag_local,np.array(adjs[:,0]),flat=True)
    adjs2=mb.tag_get_data(tag_local,np.array(adjs[:,1]),flat=True)
    adjsg1=mb.tag_get_data(ID_reordenado_tag,np.array(adjs[:,0]),flat=True)
    adjsg2=mb.tag_get_data(ID_reordenado_tag,np.array(adjs[:,1]),flat=True)
    Gids=mb.tag_get_data(ID_reordenado_tag,conj_vols,flat=True)
    lista.append(Gids)
    lista.append(all_int_fac)
    lista.append(adjs1)
    lista.append(adjs2)
    lista.append(adjsg1)
    lista.append(adjsg2)

def mount_graph(adjacencies, edges_dists, n_nodes):

    new_adjacencies = adjacencies[adjacencies[:, 1] != -1]
    graph = sp.csc_matrix((edges_dists, (new_adjacencies[:, 0], new_adjacencies[:, 1])), shape=(n_nodes, n_nodes))
    return graph

def get_shortest_path(graph, from_node, to_node):

    D, Pr = shortest_path(graph, directed=False, method='D', return_predecessors=True)
    path = get_Path(Pr, from_node, to_node)
    return path

def get_Path(Pr, from_node, to_node):
    i = from_node
    j = to_node
    path = [j]
    k = j
    while Pr[i, k] != -9999:
        path.append(Pr[i, k])
        k = Pr[i, k]
    return np.array(path[::-1])

def get_local_shortest_path_for_create_dual_edges(adjacencies, dists, coarse_gid, coarse_adj_fine, fine_vertice, fine_edge, node_to_remove: int=-1):
    
    test = (coarse_adj_fine[:, 0] == coarse_gid) & (coarse_adj_fine[:, 1] == coarse_gid)
    if node_to_remove == -1:
        test3 = test
        local_adjacencies = adjacencies[test3]
    else:
        test2 = (adjacencies[:, 0] == node_to_remove) | (adjacencies[:, 1] == node_to_remove)
        test2 = ~test2
        test3 = test2 & test
        local_adjacencies = adjacencies[test3]
    
    faces = np.unique(local_adjacencies)
    n_faces = len(faces)
    local_faces = np.arange(n_faces)
    
    local_map = np.repeat(-1, faces.max() + 1)
    local_map[faces] = local_faces
    
    local_adjacencies_remapped = local_map[local_adjacencies]
    lr = local_adjacencies_remapped
    
    dists_local = dists[test3].sum(axis=1)
    
    vertice_r = local_map[fine_vertice][0]
    edge_r = local_map[fine_edge]
    
    local_graph = sp.csr_matrix((dists_local, (lr[:, 0], lr[:, 1])), shape=(n_faces, n_faces))
    D, Pr = shortest_path(local_graph, directed=False, method='D', return_predecessors=True)
    path = get_Path(Pr, vertice_r, edge_r)
    path = np.setdiff1d(path, [vertice_r, edge_r])
    path = faces[np.isin(local_faces, path)]
    
    return path
     
    
def mount_graph(adjacencies, dists, coarse_gid, coarse_adj_fine, fine_vertice, fine_edge):
            
            test = (coarse_adj_fine[:, 0] == coarse_gid) & (coarse_adj_fine[:, 1] == coarse_gid)
            local_adjacencies = adjacencies[test]
            
            faces = np.unique(local_adjacencies)
            n_faces = len(faces)
            local_faces = np.arange(n_faces)
            
            local_map = np.repeat(-1, faces.max() + 1)
            local_map[faces] = local_faces
            
            local_adjacencies_remapped = local_map[local_adjacencies]
            lr = local_adjacencies_remapped
            
            dists_local = dists[test].sum(axis=1)
            
            vertice_r = local_map[fine_vertice][0]
            edge_r = local_map[fine_edge][0]
            
            local_graph = sp.csr_matrix((dists_local, (lr[:, 0], lr[:, 1])), shape=(n_faces, n_faces))
            D, Pr = shortest_path(local_graph, directed=False, method='D', return_predecessors=True)
            path = get_Path(Pr, vertice_r, edge_r)
            path = np.setdiff1d(path, [vertice_r, edge_r])
            path = faces[np.isin(local_faces, path)]
            
            return path
    
