import numpy as np

'''
centroids_volumes: np.ndarray, 
volumes_dimension: np.ndarray, 
centroids_nodes: np.ndarray, 
cr: np.ndarray, 
adjacencies_internal_faces: np.ndarray
'''

def create_coarse_volumes(fine_volumes, centroids_fine_volumes, primal_id, adjacencies_internal_faces, fine_nodes, fine_volumes_nodes, centroids_fine_nodes):
    coarse_nodes = define_coarse_nodes(primal_id, fine_volumes_nodes, centroids_fine_nodes)

def define_coarse_nodes(primal_id, fine_volumes_nodes, centroids_fine_nodes):
    unique_pid = np.unique(primal_id)
    for pid in unique_pid:
        gids = np.argwhere(primal_id == pid).flatten()
        centroids_in_primal_id = centroids_fine_nodes[gids]
        
        import pdb; pdb.set_trace()
    
    