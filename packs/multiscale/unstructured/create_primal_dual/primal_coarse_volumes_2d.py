# from packs.manager.meshmanager import MeshProperty
from packs.manager import MeshProperty, MeshData
from packs.utils.calculate_face_properties import sort_radial_sweep

from packs.utils.utils_old import get_box_v2_2d
import numpy as np
from shapely import geometry

# from shapely.geometry import Point, Polygon
from packs import defnames
import scipy.sparse as sp
from scipy.sparse.csgraph import shortest_path
from packs.utils import utils_old
import time
# from packs.utils.to_compile_funcs.funcs import nbspatial
# from matplotlib.path import Path

def _define_fine_volumes_in_coarse_volumes_step1(all_fine_faces_centroids, all_fine_faces_ids, coarse_points_face, coarse_nodes_centroids):
    coarse_centroids_of_face = coarse_nodes_centroids[coarse_points_face]
    p0 = coarse_centroids_of_face.min(axis=0)
    p1 = coarse_centroids_of_face.max(axis=0)
    
    limites = np.array([p0, p1])
    fine_faces_step1 = all_fine_faces_ids[get_box_v2_2d(all_fine_faces_centroids, limites)]
    
    return fine_faces_step1

def _get_fine_faces_in_triangle(centroids_coarse_points_face, centroids_fine_faces, fine_faces_step1):
    # poly = Path(centroids_coarse_points_face, closed=True)
    # test = poly.contains_points(centroids_fine_faces)
    # fine_faces_in_coarse_face = fine_faces_step1[test]

    index = sort_radial_sweep(centroids_coarse_points_face, np.arange(centroids_coarse_points_face.shape[0]))
    centroids2 = centroids_coarse_points_face[index]

    poly = geometry.Polygon(centroids2)    
    # points_list = [geometry.Point(i[0], i[1]) for i in centroids_fine_faces]
    points_list = geometry.MultiPoint(centroids_fine_faces)
    test = np.array([poly.contains(i) for i in points_list.geoms])
    fine_faces_in_coarse_face = fine_faces_step1[test]

    # polygon = centroids2
    # result = np.array([nbspatial.ray_tracing(point[0], point[1], polygon) for point in centroids_fine_faces])
    # fine_faces_in_coarse_face = fine_faces_step1[result]

    return fine_faces_in_coarse_face

def _get_fine_faces_in_polygon(centroids_coarse_points_face, centroids_fine_faces, fine_faces_step1, coarse_centroid):
    nvec = centroids_coarse_points_face.shape[0]
    v1 = np.arange(nvec + 1)
    v1[-1] = -1
    fine_faces_in_coarse_face = []

    for i in range(nvec):
        coarse_points = np.array([
            centroids_coarse_points_face[v1[i]],
            centroids_coarse_points_face[v1[i+1]],
            coarse_centroid
        ])
        v2 = _get_fine_faces_in_triangle(
            coarse_points,
            centroids_fine_faces,
            fine_faces_step1
        )
        fine_faces_in_coarse_face.append(v2)
    
    fine_faces_in_coarse_face = np.unique(np.concatenate(fine_faces_in_coarse_face))
    return fine_faces_in_coarse_face





def _get_faces_in_coarse_volume(fine_faces_step1, all_fine_faces_centroids, coarse_points_face, coarse_nodes_centroids, coarse_centroid):
    centroids_coarse_points_face = coarse_nodes_centroids[coarse_points_face]
    centroids_fine_faces = all_fine_faces_centroids[fine_faces_step1]
    
    # poly = geometry.Polygon(centroids_coarse_points_face)
    # points_list = [geometry.Point(i[0], i[1]) for i in centroids_fine_faces]
    # test = np.array([poly.contains(i) for i in points_list])
    # fine_faces_in_coarse_face = fine_faces_step1[test]

    # polygon = centroids_coarse_points_face
    # result = np.array([nbspatial.ray_tracing(point[0], point[1], polygon) for point in centroids_fine_faces])
    # fine_faces_in_coarse_face = fine_faces_step1[result]

    fine_faces_in_coarse_face = _get_fine_faces_in_triangle(
            centroids_coarse_points_face, 
            centroids_fine_faces, 
            fine_faces_step1
        )

    # if centroids_coarse_points_face.shape[0] < 4:
    #     fine_faces_in_coarse_face = _get_fine_faces_in_triangle(
    #         centroids_coarse_points_face, 
    #         centroids_fine_faces, 
    #         fine_faces_step1
    #     )
    # else:
    #     fine_faces_in_coarse_face = _get_fine_faces_in_polygon(
    #         centroids_coarse_points_face,
    #         centroids_fine_faces,
    #         fine_faces_step1,
    #         coarse_centroid
    #     )

    
    # poly = Path(centroids_coarse_points_face, closed=True)
    # test = poly.contains_points(centroids_fine_faces)
    # fine_faces_in_coarse_face = fine_faces_step1[test]



    # print(fine_faces_in_coarse_face)
    # print(fine_faces_in_coarse_face.shape[0])
    # import pdb; pdb.set_trace()

    return fine_faces_in_coarse_face

def verify_fine_adjacencies(fine_primal_ids, coarse_faces_ids, fine_adjacencies, fine_faces_of_faces, fine_faces_centroids, coarse_faces_centroids, coarse_faces_of_faces):
    fine_coarse_adjacencies = fine_primal_ids[fine_adjacencies]
    fine_coarse_adjacencies[fine_adjacencies[:, 1] == -1, 1] = -1
    
    for coarse_face in coarse_faces_ids:
        test1 = fine_coarse_adjacencies[:, 0] == coarse_face
        test2 = fine_coarse_adjacencies[:, 1] == coarse_face
        test3 = fine_coarse_adjacencies[:, 0] != coarse_face
        test4 = fine_coarse_adjacencies[:, 1] != coarse_face
        
        test5 = test1 & test4
        test6 = test2 & test3
        
        test7 = test6 | test5        
        test_fine_faces = fine_adjacencies[test7][fine_adjacencies[test7] != -1].flatten()
        test_fine_faces = np.unique(test_fine_faces[fine_primal_ids[test_fine_faces] == coarse_face])
        adjs = fine_faces_of_faces[test_fine_faces]
        
        for i, faces in enumerate(adjs):
            coarse_ids = fine_primal_ids[faces]
            test = coarse_ids == coarse_face
            if np.any(test):
                pass
            else:

                coarse_ids_candidatos = np.intersect1d(coarse_ids, coarse_faces_of_faces[coarse_face])         
                dists = np.linalg.norm(
                    fine_faces_centroids[test_fine_faces[i]] - coarse_faces_centroids[coarse_ids_candidatos], 
                    axis=1
                )
                
                coarse_id_selected = coarse_ids_candidatos[
                    dists <= dists.min()
                ]
                
                fine_primal_ids[test_fine_faces[i]] = coarse_id_selected[0]

def coarse_integrity_check(face_test, fine_adjacencies, fine_faces_centroids, fine_faces_id, primal_id, coarse_face, centroid_coarse_face, fine_edges_id, fine_edges_centroids):
    
    fine_faces_in_coarse_face = fine_faces_id[primal_id == coarse_face]
    coarse_centroid = np.mean(
        fine_faces_centroids[fine_faces_in_coarse_face],
        axis=0
    )
    # dists = np.linalg.norm(centroid_coarse_face - fine_faces_centroids[fine_faces_in_coarse_face], axis=1)
    dists = np.linalg.norm(coarse_centroid - fine_faces_centroids[fine_faces_in_coarse_face], axis=1)
    fine_face_in_coarse_centroid = fine_faces_in_coarse_face[dists <= dists.min()]
    
    edges_in_coarse_face = fine_edges_id[
        (primal_id[fine_adjacencies[:, 0]] == coarse_face) | 
        (primal_id[fine_adjacencies[:, 1]] == coarse_face)
    ]
    fine_adjacencies_local = fine_adjacencies[edges_in_coarse_face]
    fine_faces_for_local_graph = np.unique(fine_adjacencies_local)
    n_faces = fine_faces_for_local_graph.shape[0]

    local_map = np.arange(fine_faces_for_local_graph.shape[0])
    local_adjacencies = fine_adjacencies_local.copy()
    local_adjacencies[:, 0] = np.array([local_map[fine_faces_for_local_graph==i][0] for i in fine_adjacencies_local[:, 0]])
    local_adjacencies[:, 1] = np.array([local_map[fine_faces_for_local_graph==i][0] for i in fine_adjacencies_local[:, 1]])
    lr = local_adjacencies

    local_dists = np.zeros(local_adjacencies.shape)
    local_dists[:, 0] = np.linalg.norm(fine_edges_centroids[edges_in_coarse_face] - fine_faces_centroids[fine_adjacencies[edges_in_coarse_face, 0]], axis=1)
    local_dists[:, 1] = np.linalg.norm(fine_edges_centroids[edges_in_coarse_face] - fine_faces_centroids[fine_adjacencies[edges_in_coarse_face, 1]], axis=1)
    local_dists = local_dists.sum(axis=1)

    local_graph = sp.csr_matrix((local_dists, (lr[:, 0], lr[:, 1])), shape=(n_faces, n_faces))
    D, Pr = shortest_path(local_graph, directed=False, method='D', return_predecessors=True)
    from_node = local_map[fine_faces_for_local_graph==fine_face_in_coarse_centroid[0]][0]
    target_node = local_map[fine_faces_for_local_graph==face_test][0]
    local_path = utils_old.get_Path(Pr, from_node, target_node)
    path = np.array([fine_faces_for_local_graph[local_map==i][0] for i in local_path]).astype(np.int64)

    return path

def verify_fine_adjacencies_v2(fine_faces_id, primal_id, fine_adjacencies, fine_faces_centroids, coarse_faces_centroids, coarse_faces_id, fine_edges_id, fine_bool_boundary_edges, fine_edges_centroids, fine_faces_of_faces):

    fine_bool_internal_edges = ~fine_bool_boundary_edges
    
    for coarse_face in coarse_faces_id:
        t0 = time.time()
        centroid_coarse_face = coarse_faces_centroids[coarse_face]
        
        edges_in_boundary_coarse_face = fine_edges_id[
            (
                (primal_id[fine_adjacencies[:, 0]] == coarse_face) & 
                (primal_id[fine_adjacencies[:, 1]] != coarse_face)
            ) |
            (
                (primal_id[fine_adjacencies[:, 1]] == coarse_face) & 
                (primal_id[fine_adjacencies[:, 0]] != coarse_face)
            ) & 
            fine_bool_internal_edges
        ]

        faces_in_boundary_coarse_face = np.unique(fine_adjacencies[edges_in_boundary_coarse_face])
        faces_for_test = faces_in_boundary_coarse_face[
            primal_id[faces_in_boundary_coarse_face] == coarse_face
        ]

        # import pdb; pdb.set_trace()
        # faces_of_faces_for_test = fine_faces_of_faces[faces_for_test]
        # test = np.array([utils_old.allUnique(primal_id[i]) for i in faces_of_faces_for_test])
        # faces_for_test = faces_for_test[test]


        for face_test in faces_for_test:
            path = coarse_integrity_check(
                face_test,
                fine_adjacencies,
                fine_faces_centroids,
                fine_faces_id,
                primal_id,
                coarse_face,
                centroid_coarse_face,
                fine_edges_id,
                fine_edges_centroids
            )

            primal_id[path] = coarse_face
        t1 = time.time()

        print(f"Coarse face: {coarse_face} finished with {t1-t0} seconds \n")


def update_fine_faces(
        fine_faces_ids,
        coarse_centroids,
        fine_faces_centroids,
        faces_without_primal_id,
        coarse_faces_id,
        fine_adjacencies,
        fine_edges_centroids,
        fine_bool_boundary_edges,
        fine_edges_ids,
        primal_id
    ):

    fine_bool_internal_edges = ~fine_bool_boundary_edges
    fine_internal_edges = fine_edges_ids[fine_bool_internal_edges]

    internal_adjacencies = fine_adjacencies[fine_internal_edges]
    lr = internal_adjacencies
    n_faces = fine_faces_ids.shape[0]

    local_dists = np.zeros(internal_adjacencies.shape)
    local_dists[:, 0] = np.linalg.norm(fine_edges_centroids[fine_internal_edges] - fine_faces_centroids[internal_adjacencies[:, 0]], axis=1)
    local_dists[:, 1] = np.linalg.norm(fine_edges_centroids[fine_internal_edges] - fine_faces_centroids[internal_adjacencies[:, 1]], axis=1)
    local_dists = local_dists.sum(axis=1)
    local_graph = sp.csr_matrix((local_dists, (lr[:, 0], lr[:, 1])), shape=(n_faces, n_faces))
    D, Pr = shortest_path(local_graph, directed=False, method='D', return_predecessors=True)

    for face in faces_without_primal_id:
        face_centroid = fine_faces_centroids[face]
        dists = np.linalg.norm(
            coarse_centroids - face_centroid,
            axis=1
        )
        coarse_id_selected = coarse_faces_id[dists <= dists.min()]
        coarse_centroid = coarse_centroids[coarse_id_selected]
        dists = np.linalg.norm(
            fine_faces_centroids - coarse_centroid,
            axis=1
        )
        fine_face_in_coarse_centroid = fine_faces_ids[dists <= dists.min()][0]
        from_node = fine_face_in_coarse_centroid
        target_node = face
        path = utils_old.get_Path(Pr, from_node, target_node)
        primal_id[path] = coarse_id_selected


def create_coarse_volumes(faces_id_level0, faces_centroids_level0, faces_ids_level1, nodes_centroids_level1, nodes_of_faces_level1, adjacencies_level0, faces_of_faces_level0, faces_centroids_level1, faces_of_faces_level1, level:int, edges_ids_level0, bool_boundary_edges_level0, edges_centroids_level0):
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
            coarse_nodes_centroids=nodes_centroids_level1,
            coarse_centroid=faces_centroids_level1[coarse_face]
        )

        # fine_faces_in_coarse_face = _get_faces_in_coarse_volume(
        #     faces_id_level0,
        #     faces_centroids_level0,
        #     coarse_points_face=coarse_points_face,
        #     coarse_nodes_centroids=nodes_centroids_level1,
        # )

        fine_primal_ids[fine_faces_in_coarse_face] = coarse_face
    
    test = fine_primal_ids == -1

    if np.any(test):
        # TODO : criar funcao para verificar as faces que nao foram definidas em algum coarse volume
        # fine_faces_without_primal_id = faces_id_level0[test]
        # update_fine_faces(
        #     fine_faces_ids=faces_id_level0,
        #     coarse_centroids=faces_centroids_level1,
        #     fine_faces_centroids=faces_centroids_level0,
        #     faces_without_primal_id=fine_faces_without_primal_id,
        #     coarse_faces_id=faces_ids_level1,
        #     fine_adjacencies=adjacencies_level0,
        #     fine_edges_centroids=edges_centroids_level0,
        #     fine_bool_boundary_edges=bool_boundary_edges_level0,
        #     fine_edges_ids=edges_ids_level0,
        #     primal_id=fine_primal_ids
        # )
        raise NotImplementedError
    
    # verify_fine_adjacencies(
    #     fine_primal_ids=fine_primal_ids,
    #     coarse_faces_ids=faces_ids_level1,
    #     fine_adjacencies=adjacencies_level0,
    #     fine_faces_of_faces=faces_of_faces_level0,
    #     fine_faces_centroids=faces_centroids_level0,
    #     coarse_faces_centroids=faces_centroids_level1,
    #     coarse_faces_of_faces=faces_of_faces_level1
    # )

    verify_fine_adjacencies_v2(
        fine_faces_id=faces_id_level0,
        primal_id=fine_primal_ids,
        fine_adjacencies=adjacencies_level0,
        fine_faces_centroids=faces_centroids_level0,
        coarse_faces_centroids=faces_centroids_level1,
        coarse_faces_id=faces_ids_level1,
        fine_edges_id=edges_ids_level0,
        fine_bool_boundary_edges=bool_boundary_edges_level0,
        fine_edges_centroids=edges_centroids_level0,
        fine_faces_of_faces=faces_of_faces_level0
    )

    data = {
        defnames.get_primal_id_name_by_level(level): fine_primal_ids
    }
    
    return data
        