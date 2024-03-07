from packs.manager.meshmanager import MeshProperty
from packs.utils.utils_old import get_box_v2_2d
import numpy as np
from shapely.geometry import Point, Polygon

def _define_fine_volumes_in_coarse_volumes_step1(all_fine_faces_centroids, all_fine_faces_ids, coarse_points_face, coarse_nodes_centroids):
    coarse_centroids_of_face = coarse_nodes_centroids[coarse_points_face]
    p0 = coarse_centroids_of_face.min(axis=0)
    p1 = coarse_centroids_of_face.max(axis=0)
    limites = np.array([p0, p1])
    fine_faces_step1 = all_fine_faces_ids[get_box_v2_2d(all_fine_faces_centroids, limites)]
    return fine_faces_step1

def _get_faces_in_coarse_volume(fine_faces_step1, all_fine_faces_centroids, coarse_points_face, coarse_nodes_centroids):
    centroids_coarse_points_face = coarse_nodes_centroids[coarse_points_face]    
    centroids_fine_faces = all_fine_faces_centroids[fine_faces_step1]
    
    poly = Polygon(centroids_coarse_points_face)
    points_list = [Point(i) for i in centroids_fine_faces]
    test = np.array([poly.contains(i) for i in points_list])
    fine_faces_in_coarse_face = fine_faces_step1[test]
    return fine_faces_in_coarse_face

def verify_fine_adjacencies(fine_primal_ids, coarse_faces_ids, fine_faces_of_nodes, fine_faces_ids, fine_adjacencies):
    fine_coarse_adjacencies = fine_primal_ids[fine_adjacencies]
    
    for coarse_face in coarse_faces_ids:
        test1 = fine_coarse_adjacencies[:, 0] == coarse_face
        test2 = fine_coarse_adjacencies[:, 1] == coarse_face
        test3 = fine_coarse_adjacencies[:, 0] != coarse_face
        test4 = fine_coarse_adjacencies[:, 1] != coarse_face
        
        test5 = test1 & test4
        test6 = test2 & test3
        
        test7 = test6 | test5        
        test_fine_faces = fine_adjacencies[test7][fine_adjacencies[test7] != -1]
        
        
        import pdb; pdb.set_trace()
        
        
        
        
        
        
        
        
    

def create_coarse_volumes(faces_id_level0, faces_centroids_level0, faces_ids_level1, nodes_centroids_level1, nodes_of_faces_level1, faces_of_nodes_level0, adjacencies_level0):
    """Insert the 'primal_fine_ids' tag in mesh_properties_level0

    Args:
        mesh_properties_level0 (MeshProperty): level0 mesh_properties
        mesh_properties_level1 (MeshProperty): level1 mesh_properties

    Raises:
        NotImplementedError: _description_

    Returns:
        _type_: _description_
    """
    
    fine_primal_ids = np.repeat(-1, faces_id_level0.shape[0])
    # mesh_delta = mesh_properties_level0['h_dist'].min()/4
    
    for coarse_face in faces_ids_level1:
        coarse_points_face = nodes_of_faces_level1[coarse_face]
        fine_faces_step1 = _define_fine_volumes_in_coarse_volumes_step1(
            all_fine_faces_centroids=faces_centroids_level0,
            all_fine_faces_ids=faces_id_level0,
            coarse_points_face=coarse_points_face,
            coarse_nodes_centroids=nodes_centroids_level1
        )
        
        fine_faces_in_coarse_face = _get_faces_in_coarse_volume(
            fine_faces_step1,
            faces_centroids_level0,
            coarse_points_face=coarse_points_face,
            coarse_nodes_centroids=nodes_centroids_level1
        )
        fine_primal_ids[fine_faces_in_coarse_face] = coarse_face
    
    test = fine_primal_ids == -1
    if np.any(test):
        # TODO : criar funcao para verificar as faces que nao foram definidas em algum coarse volume
        raise NotImplementedError
    
    verify_fine_adjacencies(
        fine_primal_ids=fine_primal_ids,
        coarse_faces_ids=faces_ids_level1,
        fine_faces_of_nodes=faces_of_nodes_level0,
        fine_faces_ids=faces_id_level0,
        fine_adjacencies=adjacencies_level0
    )
    
    return fine_primal_ids
        