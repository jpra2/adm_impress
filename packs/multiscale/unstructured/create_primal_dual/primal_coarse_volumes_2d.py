# from packs.manager.meshmanager import MeshProperty
from packs.manager import MeshProperty, MeshData
from packs.utils.calculate_face_properties import sort_radial_sweep
from packs.manager.generic_data import PrimalCoarseData
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation

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

def get_coarse_structure(
        level: int,
        primal_id: np.ndarray,
        faces_id_level0: np.ndarray,
        adjacencies_level0: np.ndarray,
        edges_id_level0: np.ndarray,
        nodes_of_edges_level0: np.ndarray,
        bool_boundary_edges_level0: np.ndarray,
        bool_boundary_nodes_level0: np.ndarray,
        nodes_level0: np.ndarray,
        nodes_weight: np.ndarray,
        nodes_of_nodes_level0: np.ndarray,
        edges_of_nodes_level0: np.ndarray,
        faces_of_nodes_level0: np.ndarray,
        nodes_centroids_level0: np.ndarray,
        faces_centroids_level0: np.ndarray,
        permeability: np.ndarray,
        unitary_normal_edges: np.ndarray,
        dual_id_level0: np.ndarray,
        edges_dim_level0: np.ndarray,
        lsds: LsdsFluxCalculation
):
    
    
    cids = np.unique(primal_id)
    name = 'level_' + str(level) + '_c_'
    cadj_fine = primal_id[adjacencies_level0]
    cadj_fine[adjacencies_level0 == -1] = -1
    bedges = edges_id_level0[bool_boundary_edges_level0]
    bnodes = nodes_level0[bool_boundary_nodes_level0]

    resp = []

    for cid in cids:
        cname = name + str(cid)
        coarse_data = PrimalCoarseData(cname)
        faces_in = faces_id_level0[primal_id==cid]
        set_faces_in = set(faces_in)


        edges_in = edges_id_level0[
            (cadj_fine[:, 0] == cid) | (cadj_fine[:, 1] == cid)  
        ]

        nodes_in = np.unique(nodes_of_edges_level0[edges_in].flatten())

        local_edges = np.arange(edges_in.shape[0])
        local_faces = np.arange(faces_in.shape[0])
        local_nodes = np.arange(nodes_in.shape[0])
        local_bool_boundary_edges = np.isin(edges_in, bedges)
        local_bool_boundary_nodes = np.isin(nodes_in, bnodes)

        
        local_nodes_of_edges = nodes_of_edges_level0[edges_in]
        n1, n2 = local_nodes_of_edges.shape
        for i in range(n1):
            for j in range(n2):
                aux = local_nodes[nodes_in==local_nodes_of_edges[i,j]]
                local_nodes_of_edges[i,j] = aux

        nodes_of_nodes_aux = nodes_of_nodes_level0[nodes_in]
        local_nodes_of_nodes = []
        for aux in nodes_of_nodes_aux:
            test = np.isin(aux, nodes_in)
            nodes_aux = aux[test]
            nodes_aux = np.array([local_nodes[nodes_in==i] for i in nodes_aux])
            local_nodes_of_nodes.append(nodes_aux.flatten())
        local_nodes_of_nodes = np.array(local_nodes_of_nodes, dtype='O')

        edges_of_nodes_aux = edges_of_nodes_level0[nodes_in]
        local_edges_of_nodes = []
        for aux in edges_of_nodes_aux:
            test = np.isin(aux, edges_in)
            edges_aux = aux[test]
            edges_aux = np.array([local_edges[edges_in==i] for i in edges_aux])
            local_edges_of_nodes.append(edges_aux.flatten())
        local_edges_of_nodes = np.array(local_edges_of_nodes, dtype='O')

        faces_of_nodes_aux = faces_of_nodes_level0[nodes_in]
        local_faces_of_nodes = []
        for aux in faces_of_nodes_aux:
            test = np.isin(aux, faces_in)
            faces_aux = aux[test]
            faces_aux = np.array([local_faces[faces_in==i] for i in faces_aux])
            local_faces_of_nodes.append(faces_aux.flatten())
        local_faces_of_nodes = np.array(local_faces_of_nodes, dtype='O')
        
        intersect_edges = edges_id_level0[
            ((cadj_fine[:, 0] == cid) & (cadj_fine[:, 1] != cid)) |
            ((cadj_fine[:, 0] != cid) & (cadj_fine[:, 1] == cid))
        ]

        local_bool_intersect_edges = np.isin(edges_in, intersect_edges)
        intersect_nodes = np.unique(nodes_of_edges_level0[intersect_edges].flatten())
        local_bool_intersect_nodes = np.isin(nodes_in, intersect_nodes)
        local_bool_boundary_edges = local_bool_boundary_edges | local_bool_intersect_edges
        local_bool_boundary_nodes = local_bool_boundary_nodes | local_bool_intersect_nodes

        local_adjacencies = adjacencies_level0[edges_in].copy()

        for edge in local_edges[local_bool_intersect_edges]:       
            if set([local_adjacencies[edge, 1]]) & set_faces_in:
                local_adjacencies[edge] = local_adjacencies[edge, [1, 0]]
            local_adjacencies[edge, 1] = -1
        
        # test_bound = local_adjacencies == -1
        ni, nj = local_adjacencies.shape
        for i in range(ni):
            for j in range(nj):
                if local_adjacencies[i,j] == -1:
                    continue
                aux = local_faces[faces_in==local_adjacencies[i,j]]
                local_adjacencies[i,j] = aux

        test1 = np.isin(nodes_weight['node_id'], nodes_in)
        test2 =  np.isin(nodes_weight['face_id'], faces_in)
        test3 = test1 & test2

        local_nodes_weight = nodes_weight[test3]

        for node in nodes_in:
            local_nodes_weight['node_id'][local_nodes_weight['node_id']==node] = local_nodes[nodes_in==node]
        
        for face in faces_in:
            local_nodes_weight['face_id'][local_nodes_weight['face_id']==face] = local_faces[faces_in==face]
        
        test4 = np.isin(local_nodes_weight['node_id'], local_nodes[local_bool_boundary_nodes])
        test4 = ~test4
        local_nodes_weight = local_nodes_weight[test4]
        
        coarse_id = np.array([cid])

        nodes_to_calculate = local_nodes[local_bool_intersect_nodes | local_bool_boundary_nodes]

        local_xi_params = lsds.get_all_edges_flux_params(
            faces_centroids_level0[faces_in],
            local_bool_boundary_edges,
            nodes_centroids_level0[nodes_in],
            local_nodes_of_edges,
            local_adjacencies,
            local_faces,
            local_edges,
            unitary_normal_edges[edges_in],
            permeability[faces_in],
            edges_dim_level0[edges_in]
        )

        coarse_data.insert_or_update_data({
            coarse_data.my_data_names[0]: local_adjacencies,
            coarse_data.my_data_names[1]: local_nodes,
            coarse_data.my_data_names[2]: local_faces,
            coarse_data.my_data_names[3]: local_edges,
            coarse_data.my_data_names[4]: local_bool_boundary_edges,
            coarse_data.my_data_names[5]: local_bool_intersect_edges,
            coarse_data.my_data_names[6]: nodes_in,
            coarse_data.my_data_names[7]: faces_in,
            coarse_data.my_data_names[8]: edges_in,
            coarse_data.my_data_names[9]: local_nodes_weight,
            coarse_data.my_data_names[10]: local_bool_boundary_nodes,
            coarse_data.my_data_names[11]: local_bool_intersect_nodes,
            coarse_data.my_data_names[12]: coarse_id,
            coarse_data.my_data_names[13]: nodes_to_calculate,
            coarse_data.my_data_names[14]: local_nodes_of_nodes,
            coarse_data.my_data_names[15]: local_edges_of_nodes,
            coarse_data.my_data_names[16]: local_faces_of_nodes,
            coarse_data.my_data_names[17]: nodes_centroids_level0[nodes_in],
            coarse_data.my_data_names[18]: faces_centroids_level0[faces_in],
            coarse_data.my_data_names[19]: permeability[faces_in],
            coarse_data.my_data_names[20]: unitary_normal_edges[edges_in],
            coarse_data.my_data_names[21]: np.array([]),
            coarse_data.my_data_names[22]: np.array([]),
            coarse_data.my_data_names[23]: np.array([]),
            coarse_data.my_data_names[24]: np.array([]),
            coarse_data.my_data_names[25]: dual_id_level0[faces_in],
            coarse_data.my_data_names[26]: np.array([]),
            coarse_data.my_data_names[27]: np.array([]),
            coarse_data.my_data_names[28]: edges_dim_level0[edges_in],
            coarse_data.my_data_names[29]: local_xi_params['xi_params'],
            coarse_data.my_data_names[30]: local_xi_params['xi_params'],
            coarse_data.my_data_names[31]: test3,
            coarse_data.my_data_names[32]: local_nodes_of_edges,
        })

        coarse_data.export_data()   
        resp.append(coarse_data)

    return resp     

def load_coarse_structure(level, primal_ids):
    cids = np.unique(primal_ids)
    name = 'level_' + str(level) + '_c_'

    resp = []

    for cid in cids:
        cname = name + str(cid)
        coarse_data = PrimalCoarseData(cname)
        coarse_data.load_data()
        resp.append(coarse_data)
    
    return resp




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
        