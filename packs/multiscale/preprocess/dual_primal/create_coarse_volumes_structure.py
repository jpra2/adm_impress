import numpy as np

'''
centroids_volumes: np.ndarray, 
volumes_dimension: np.ndarray, 
centroids_nodes: np.ndarray, 
cr: np.ndarray, 
adjacencies_internal_faces: np.ndarray
'''

def create_coarse_volumes(fine_volumes, centroids_fine_volumes, primal_id, adjacencies_internal_faces, fine_nodes, fine_volumes_nodes, centroids_fine_nodes):
    coarse_nodes, coarse_gid_nodes, coarse_nodes_centroids, coarse_centroids = define_coarse_nodes(primal_id, fine_volumes_nodes, centroids_fine_nodes, centroids_fine_volumes)
    import pdb; pdb.set_trace()

def define_coarse_nodes(primal_id, fine_volumes_nodes, centroids_fine_nodes, centroids_fine_volumes):
    unique_pid = np.unique(primal_id)
    all_coarse_nodes = set()
    coarse_volumes_nodes = []
    coarse_centroids = []
    dtype_org = [('x', np.float), ('y', np.float), ('z', np.float)]
    
    for pid in unique_pid:
        gids = np.argwhere(primal_id == pid).flatten()
        centroids_in_primal_id = centroids_fine_nodes[fine_volumes_nodes[gids]]
        centroids_in_primal_id = centroids_in_primal_id.reshape(
            (
                centroids_in_primal_id.shape[0]*centroids_in_primal_id.shape[1],
                centroids_in_primal_id.shape[2]
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
        
        cents_tuple = [tuple(cent) for cent in cents]
        coarse_volumes_nodes.append(cents_tuple)
        all_coarse_nodes.update(set(cents_tuple))
        coarse_centroid = np.mean(centroids_fine_volumes[gids], axis=0)
        coarse_centroids.append(coarse_centroid)
    
    coarse_centroids = np.array(coarse_centroids)
    all_coarse_nodes = np.array(list(all_coarse_nodes))
    all_coarse_nodes_org = np.zeros(len(all_coarse_nodes), dtype=dtype_org)
    all_coarse_nodes_org['x'] = all_coarse_nodes[:,0]
    all_coarse_nodes_org['y'] = all_coarse_nodes[:,1]
    all_coarse_nodes_org['z'] = all_coarse_nodes[:,2]
    
    all_coarse_nodes_org = np.sort(all_coarse_nodes_org, order = ['x', 'y', 'z'])
    all_coarse_nodes[:,0] = all_coarse_nodes_org['x']
    all_coarse_nodes[:,1] = all_coarse_nodes_org['y']
    all_coarse_nodes[:,2] = all_coarse_nodes_org['z']
    
    cents_tuple = [tuple(cent) for cent in all_coarse_nodes]
    coarse_gid_nodes = np.arange(len(cents_tuple))
    map_gid = dict(zip(cents_tuple, coarse_gid_nodes))
    
    for coarse_nodes in coarse_volumes_nodes:
        for i, node in enumerate(coarse_nodes):
            coarse_nodes[i] = map_gid[node]
    
    return coarse_volumes_nodes, coarse_gid_nodes, all_coarse_nodes, coarse_centroids
    
    
            
            
    
    
    
    
    
    
        
        
        
        
        
        
        
        
        
        
    import pdb; pdb.set_trace()
    
    