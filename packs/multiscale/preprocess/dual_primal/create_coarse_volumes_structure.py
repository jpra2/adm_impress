from pdb import Pdb
import numpy as np
import scipy.sparse as sp

def create_coarse_structure(centroids_fine_volumes, primal_id, adjacencies_internal_faces, fine_volumes_nodes, centroids_fine_nodes, fine_internal_faces, fine_boundary_faces, adjacencies_boundary_faces, **kwargs):
        
    """Get the coarse structure mesh

    Args:
        centroids_fine_volumes (np.ndarray((len(fine_volumes), 3))): centroids of fine volumes
        primal_id (np.ndarray(len(fine_volumes))): map of primal id volumes in fine mesh 
        adjacencies_internal_faces (np.ndarray((len(internal_faces), 2))): volumes adjacencies of internal fine faces
        fine_volumes_nodes (np.ndarray((len(fine_volumes), 8))): map fine volumes in yours fine nodes
        centroids_fine_nodes (np.ndarray((len(fine_nodes), 3))): centroid of fine nodes
        fine_internal_faces (np.ndarray(len(internal_faces))): fine internal_faces
        fine_boundary_faces (np.ndarray(len(boundary_faces))): fine boundary faces
        adjacencies_boundary_faces (np.ndarray((len(boundary_faces), 1 or 2))): volumes adjacencies of boundary fine faces

    Returns:
        resp [dict]: coarse structure mesh
    """
    resp = get_coarse_volumes_info(primal_id, centroids_fine_nodes, fine_volumes_nodes, centroids_fine_volumes)
    resp.update(get_coarse_faces_info(fine_internal_faces, fine_boundary_faces, primal_id, adjacencies_internal_faces, get_fine_boundary_adjacencies(adjacencies_boundary_faces), resp['gids']))
    
    return resp

def get_cents_tuple(centroids_nodes: np.ndarray, **kwargs):
        
    centroids_in_primal_id = centroids_nodes.reshape(
        (
            centroids_nodes.shape[0]*centroids_nodes.shape[1],
            centroids_nodes.shape[2]
        )
    )
    
    cmin = centroids_in_primal_id.min(axis=0)
    cmax = centroids_in_primal_id.max(axis=0)
    cents = np.vstack([cmin, cmax]).T
    cents = np.array(np.meshgrid(cents[0], cents[1], cents[2]))
    cents = cents.reshape(
        cents.shape[3]*cents.shape[2]*cents.shape[1],
        cents.shape[0]
    )
    
    cents_tuple: list = [tuple(cent) for cent in cents]
    
    return cents_tuple

def get_all_coarse_nodes_reordenated(all_coarse_nodes):
    dtype_org = [('x', np.float64), ('y', np.float64), ('z', np.float64)]
    
    all_coarse_nodes_org = np.zeros(len(all_coarse_nodes), dtype=dtype_org)
    all_coarse_nodes_org['x'] = all_coarse_nodes[:,0]
    all_coarse_nodes_org['y'] = all_coarse_nodes[:,1]
    all_coarse_nodes_org['z'] = all_coarse_nodes[:,2]
    
    all_coarse_nodes_org = np.sort(all_coarse_nodes_org, order = ['x', 'y', 'z'])
    all_coarse_nodes[:,0] = all_coarse_nodes_org['x']
    all_coarse_nodes[:,1] = all_coarse_nodes_org['y']
    all_coarse_nodes[:,2] = all_coarse_nodes_org['z']
    
    return all_coarse_nodes
    
def get_coarse_gid_nodes(all_coarse_nodes, coarse_volumes_nodes):
    
    cents_tuple = [tuple(cent) for cent in all_coarse_nodes]
    coarse_gid_nodes = np.arange(len(cents_tuple))
    map_gid = dict(zip(cents_tuple, coarse_gid_nodes))
    
    for coarse_nodes in coarse_volumes_nodes:
        for i, node in enumerate(coarse_nodes):
            coarse_nodes[i] = map_gid[node]
    
    return coarse_gid_nodes, np.array(coarse_volumes_nodes)

def get_coarse_volumes_info(primal_id, centroids_fine_nodes, fine_volumes_nodes, centroids_fine_volumes):
    
    unique_pid = np.unique(primal_id)
    all_coarse_nodes = set()
    coarse_volumes_nodes = []
    coarse_centroids = np.zeros((len(unique_pid), 3))
    
    update_coarse_data(unique_pid, primal_id, centroids_fine_nodes, fine_volumes_nodes, coarse_volumes_nodes, all_coarse_nodes, centroids_fine_volumes, coarse_centroids)
    
    all_coarse_nodes = get_all_coarse_nodes_reordenated(np.array(list(all_coarse_nodes)))
    coarse_gid_nodes, coarse_volumes_nodes = get_coarse_gid_nodes(all_coarse_nodes, coarse_volumes_nodes)
    
    resp = {
        'centroids_nodes': all_coarse_nodes,
        'nodes': coarse_gid_nodes,
        'volumes_nodes': coarse_volumes_nodes,
        'gids': unique_pid,
        'centroids_volumes': coarse_centroids
    }
    
    return resp

def update_coarse_data(unique_pid, primal_id, centroids_fine_nodes, fine_volumes_nodes, coarse_volumes_nodes, all_coarse_nodes, centroids_fine_volumes, coarse_centroids):
    for pid in unique_pid:
        gids = np.argwhere(primal_id == pid).flatten()
        cents_tuple = get_cents_tuple(centroids_fine_nodes[fine_volumes_nodes[gids]])
        coarse_volumes_nodes.append(cents_tuple)
        all_coarse_nodes.update(set(cents_tuple))
        coarse_centroids[pid] = np.mean(centroids_fine_volumes[gids], axis=0)

def get_coarse_faces_info(fine_internal_faces, fine_boundary_faces, primal_id, adjacencies_internal_faces, adjacencies_boundary_faces, unique_pid):
    
    all_fine_faces = np.sort(np.concatenate([fine_internal_faces, fine_boundary_faces]))
    map_fine_internal_faces = np.repeat(-1, len(all_fine_faces))
    map_fine_internal_faces[fine_internal_faces] = np.arange(len(fine_internal_faces))
    adjs_coarse = primal_id[adjacencies_internal_faces]
    adjs_coarse_ordenated = np.sort(adjs_coarse, axis=1)
    unique_coarse_adjacencies = np.unique(adjs_coarse_ordenated, axis=0)
    test = ~(unique_coarse_adjacencies[:, 0] == unique_coarse_adjacencies[:, 1])
    unique_coarse_adjacencies = unique_coarse_adjacencies[test]
    del test
    coarse_internal_faces = np.arange(len(unique_coarse_adjacencies))    
    
    fine_adj_coarse_boundary_faces = primal_id[adjacencies_boundary_faces].flatten()
    unique_coarse_boundary_adjacencies = np.unique(fine_adj_coarse_boundary_faces)
    coarse_boundary_faces = np.arange(len(coarse_internal_faces), len(coarse_internal_faces)+len(unique_coarse_boundary_adjacencies))
    
    primal_id_coarse_faces = np.repeat(-1, len(all_fine_faces))
    
    define_primal_id_coarse_internal_faces(unique_coarse_adjacencies, adjs_coarse, fine_internal_faces, primal_id_coarse_faces)
    
    define_primal_id_coarse_boundary_faces(coarse_boundary_faces, unique_coarse_boundary_adjacencies, fine_adj_coarse_boundary_faces, fine_boundary_faces, primal_id_coarse_faces)
    
    all_coarse_faces, all_coarse_faces_adjs = get_coarse_faces_adjs(coarse_internal_faces, coarse_boundary_faces, unique_coarse_adjacencies, unique_coarse_boundary_adjacencies)
    
    coarse_volumes_faces = get_coarse_volumes_faces(unique_pid, all_coarse_faces, unique_coarse_adjacencies, coarse_internal_faces, unique_coarse_boundary_adjacencies, coarse_boundary_faces)
    
    resp = {
        'adjacencies_faces': all_coarse_faces_adjs,
        'boundary_faces': coarse_boundary_faces,
        'internal_faces': coarse_internal_faces,
        'faces_id': all_coarse_faces,
        'primal_id_coarse_faces': primal_id_coarse_faces,
        'volumes_faces': coarse_volumes_faces
    }
    
    return resp

def define_primal_id_coarse_internal_faces(unique_coarse_adjacencies, adjs_coarse, fine_internal_faces, primal_id_coarse_faces):
    for i, adj in enumerate(unique_coarse_adjacencies):
        test1 = (adjs_coarse[:,0] == adj[0]) & (adjs_coarse[:,1] == adj[1])
        test2 = (adjs_coarse[:,0] == adj[1]) & (adjs_coarse[:,1] == adj[0])
        test = test1 | test2
        fine_faces = fine_internal_faces[test]
        primal_id_coarse_faces[fine_faces] = i

def define_primal_id_coarse_boundary_faces(coarse_boundary_faces, unique_coarse_boundary_adjacencies, fine_adj_coarse_boundary_faces, fine_boundary_faces, primal_id_coarse_faces):
    for i, adj in zip(coarse_boundary_faces, unique_coarse_boundary_adjacencies):
        test = fine_adj_coarse_boundary_faces == adj
        boundary_faces = fine_boundary_faces[test]
        primal_id_coarse_faces[boundary_faces] = i

def get_coarse_faces_adjs(coarse_internal_faces, coarse_boundary_faces, unique_coarse_adjacencies, unique_coarse_boundary_adjacencies):
     
    all_coarse_faces = np.concatenate([coarse_internal_faces, coarse_boundary_faces])    
    n_coarse_faces = len(all_coarse_faces)
    all_coarse_faces_adjs = np.zeros((n_coarse_faces, 2), dtype=int)
    all_coarse_faces_adjs[coarse_internal_faces] = unique_coarse_adjacencies
    all_coarse_faces_adjs[coarse_boundary_faces, 0] = unique_coarse_boundary_adjacencies
    all_coarse_faces_adjs[coarse_boundary_faces, 1] = -1
    
    return all_coarse_faces, all_coarse_faces_adjs

def get_coarse_volumes_faces(unique_pid, all_coarse_faces, unique_coarse_adjacencies, coarse_internal_faces, unique_coarse_boundary_adjacencies, coarse_boundary_faces):
    n_coarse_volumes = len(unique_pid)
    n_coarse_faces = len(all_coarse_faces)
    coarse_volumes_faces = sp.lil_matrix((n_coarse_volumes, n_coarse_faces), dtype=bool).tocsc()
    coarse_volumes_faces[unique_coarse_adjacencies[:, 0], coarse_internal_faces] = True
    coarse_volumes_faces[unique_coarse_adjacencies[:, 1], coarse_internal_faces] = True
    coarse_volumes_faces[unique_coarse_boundary_adjacencies, coarse_boundary_faces] = True
    coarse_volumes_faces = sp.find(coarse_volumes_faces)
    coarse_volumes_faces = np.array([coarse_volumes_faces[0], coarse_volumes_faces[1]], dtype=int).T
    
    coarse_volumes_faces = get_coarse_volumes_faces_resp(unique_pid, coarse_volumes_faces)
    
    return coarse_volumes_faces
    
def get_coarse_volumes_faces_resp(unique_pid, coarse_volumes_faces):
    coarse_volumes_faces_resp = []
    for pid in unique_pid:
        coarse_volumes_faces_resp.append(coarse_volumes_faces[:, 1][coarse_volumes_faces[:,0] == pid])
    
    coarse_volumes_faces_resp = np.array(coarse_volumes_faces_resp, dtype='O')
    return coarse_volumes_faces_resp

def get_fine_boundary_adjacencies(adjacencies_boundary_faces: np.ndarray):
    if len(adjacencies_boundary_faces.shape) > 1:
        set_test = np.unique(adjacencies_boundary_faces[:, 0])
        if len(set([-1]) & set(set_test)) > 0:
            col = 1
        else:
            col = 0
        return adjacencies_boundary_faces[:, col].flatten()
    else:
        return adjacencies_boundary_faces.flatten()

def get_fine_volumes_by_coarse_information(coarse_structures, level, info: str):
    
    all_infos = ['primal_id', 'dual_id']
    if info not in all_infos:
        raise NotImplementedError(f'Unspected. info must be in {all_infos}')
    
    if level == 1:
        return coarse_structures[0][info]
    else:
        level_info = level-1
        n_total = len(coarse_structures)
        info_level = coarse_structures[level_info][info].copy()
        for i in range(level_info):
            n_level = n_total - (i+1)
            id_level_ant_info = coarse_structures[n_level-1]['gids']
            primal_info = coarse_structures[n_level-1]['primal_id']
            n_fine = len(primal_info)
            fine_info = np.repeat(-1, n_fine)
            for c_gid in id_level_ant_info:
                fine_info[primal_info==c_gid] = info_level[c_gid]
            
            info_level = fine_info.copy()
            
        return fine_info


    
    
    
    
    