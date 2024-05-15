from packs.manager.meshmanager import MeshProperty
from packs import defnames
import numpy as np
from packs.utils.utils_old import get_local_shortest_path_for_create_dual_edges

def create_dual_vertices(
        coarse_faces_id: np.ndarray, 
        coarse_faces_centroids: np.ndarray, 
        fine_faces_id: np.ndarray, 
        fine_faces_centroids: np.ndarray,
        primal_id: np.ndarray, 
        dual_id: np.ndarray
    ) -> None:

    dual_vertices = np.repeat(-1, coarse_faces_id.shape[0])
    for coarse_face in coarse_faces_id:
        coarse_centroid = coarse_faces_centroids[coarse_face]
        local_fine_faces = fine_faces_id[primal_id == coarse_face]
        dists = np.linalg.norm(
            fine_faces_centroids[local_fine_faces] - coarse_centroid,
            axis = 1
        )

        selected_vertice = local_fine_faces[dists <= dists.min()]
        dual_vertices[coarse_face] = selected_vertice
    
    test = dual_vertices == -1
    if np.any(test):
        raise NotImplementedError
    
    dual_id[dual_vertices] = defnames.dual_ids('vertice_id')

def create_dual_edges(
        coarse_faces_id: np.ndarray, 
        coarse_faces_centroids: np.ndarray, 
        fine_faces_id: np.ndarray, 
        fine_faces_centroids: np.ndarray, 
        primal_id: np.ndarray,
        fine_adjacencies: np.ndarray,
        coarse_adjacencies: np.ndarray,
        fine_edges: np.ndarray,
        coarse_edges: np.ndarray,
        fine_edges_centroids: np.ndarray,
        coarse_edges_centroids: np.ndarray,
        fine_faces_of_faces: np.ndarray,
        coarse_boundary_edges: np.ndarray,
        fine_boundary_edges: np.ndarray, 
        dual_id: np.ndarray
    ):

    ## first: loop in boundary edges
    for coarse_edge in coarse_boundary_edges:
        coarse_edge_centroid = coarse_edges_centroids[coarse_edge]
        coarse_face_adj = coarse_adjacencies[coarse_edge, 0]
        boundary_fine_edges_in_coarse_face = fine_boundary_edges[
            (primal_id[fine_adjacencies[fine_boundary_edges, 0]]==coarse_face_adj)
        ]
        dists = np.linalg.norm(
            fine_edges_centroids[boundary_fine_edges_in_coarse_face] - coarse_edge_centroid,
            axis=1
        )

        selected_fine_edge = boundary_fine_edges_in_coarse_face[dists <= dists.min()]
        selected_face_to_dual_edge = fine_adjacencies[selected_fine_edge, 0]

        dual_id[selected_face_to_dual_edge] = defnames.dual_ids('edge_id')
    
    # second: loop in internal_edges
    bool_coarse_boundary_edges = np.isin(coarse_edges, coarse_boundary_edges)
    bool_coarse_internal_edges = ~bool_coarse_boundary_edges

    for coarse_edge in coarse_edges[bool_coarse_internal_edges]:
        coarse_edge_centroid = coarse_edges_centroids[coarse_edge]
        coarse_faces_adj = coarse_adjacencies[coarse_edge]
        fine_internal_edges_between_coarse_faces = fine_edges[
            (
                (primal_id[fine_adjacencies[:, 0]] == coarse_faces_adj[0]) &
                (primal_id[fine_adjacencies[:, 1]] == coarse_faces_adj[1])
            ) |
            (
                (primal_id[fine_adjacencies[:, 1]] == coarse_faces_adj[0]) &
                (primal_id[fine_adjacencies[:, 0]] == coarse_faces_adj[1])
            )
        ]
        edges_centroids_between_coarses = fine_edges_centroids[
            fine_internal_edges_between_coarse_faces
        ]

        dists = np.linalg.norm(
            edges_centroids_between_coarses - coarse_edge_centroid,
            axis=1
        )

        selected_faces_to_dual_edge = fine_adjacencies[    
            fine_internal_edges_between_coarse_faces[
                dists <= dists.min()
            ]
        ]

        dual_id[selected_faces_to_dual_edge] = defnames.dual_ids('edge_id')

    ## three: loop for create paths in coarse faces
    # initial_dual_edges = fine_faces_id[
    #     dual_id==defnames.dual_ids('edge_id')
    # ]

    edges_paths = dict()

    dists = np.zeros(fine_adjacencies.shape, dtype=np.float64)
    dists[:, 0] = np.linalg.norm(
        fine_edges_centroids - fine_faces_centroids[fine_adjacencies[:, 0]],
        axis=1
    )
    dists[:, 1] = np.linalg.norm(
        fine_edges_centroids - fine_faces_centroids[fine_adjacencies[:, 1]],
        axis=1
    )
    dists[fine_adjacencies==-1] = np.inf

    for coarse_face in coarse_faces_id:
        dual_vertice_in_coarse = fine_faces_id[
            (primal_id==coarse_face) & (dual_id==defnames.dual_ids('vertice_id'))
        ]
        dual_edges_in_coarse = fine_faces_id[
            (primal_id==coarse_face) & (dual_id==defnames.dual_ids('edge_id')) 
        ]

        for dual_edge in dual_edges_in_coarse:
            path = get_local_shortest_path_for_create_dual_edges(
                fine_adjacencies,
                dists,
                coarse_face,
                primal_id[fine_adjacencies],
                dual_vertice_in_coarse,
                dual_edge
            )

            dual_id[path] = defnames.dual_ids('edge_id')
            edges_paths.update({
                (dual_vertice_in_coarse[0], dual_edge): path
            })

    # ## four: check viz of dual_vertices
    # for coarse_face in coarse_faces_id:
    #     dual_vertice = fine_faces_id[
    #         (primal_id==coarse_face) & (dual_id==defnames.dual_ids('vertice_id'))
    #     ]

    #     faces_of_dual_vertice = fine_faces_of_faces[dual_vertice[0]]
    #     test_edges_faces_of_dual_edges = dual_id[faces_of_dual_vertice] == defnames.dual_ids('edge_id')
    #     if test_edges_faces_of_dual_edges.sum() == faces_of_dual_vertice.shape[0]:
    #         pass
    #     else:
    #         faces_not_edges = faces_of_dual_vertice[
    #             dual_id[faces_of_dual_vertice] != defnames.dual_ids('edge_id')
    #         ]
    #         initial_dual_edges_coarse_face = initial_dual_edges[
    #             primal_id[initial_dual_edges] == coarse_face
    #         ]

    #         for face in faces_not_edges:
    #             dists_to_test = np.linalg.norm(
    #                 fine_faces_centroids[initial_dual_edges_coarse_face] - fine_faces_centroids[face],
    #                 axis=1
    #             )
                
    #             selected_dual_initial_face = initial_dual_edges_coarse_face[dists_to_test <= dists_to_test.min()]
    #             path_initial = edges_paths[(dual_vertice[0], selected_dual_initial_face[0])]
    #             dual_id[path_initial] = -1
    #             path = get_local_shortest_path_for_create_dual_edges(
    #                 fine_adjacencies,
    #                 dists,
    #                 coarse_face,
    #                 primal_id[fine_adjacencies],
    #                 [face],
    #                 selected_dual_initial_face[0],
    #                 node_to_remove=dual_vertice[0]
    #             )
    #             dual_id[face] = defnames.dual_ids('edge_id')
    #             dual_id[path] = defnames.dual_ids('edge_id')
    #             edges_paths.update(
    #                 {
    #                     (dual_vertice[0], selected_dual_initial_face[0]): np.concatenate([path, [face]])
    #                 }
    #             )
    
    # all_edges_paths = np.unique(np.concatenate(list(edges_paths.values())))
    # dual_id[all_edges_paths] = defnames.dual_ids('edge_id')
    dual_id[dual_id == -1] = defnames.dual_ids('face_id')

    return edges_paths

def create_dual_edges_v2(
        coarse_faces_id: np.ndarray, 
        coarse_faces_centroids: np.ndarray, 
        fine_faces_id: np.ndarray, 
        fine_faces_centroids: np.ndarray, 
        primal_id: np.ndarray,
        fine_adjacencies: np.ndarray,
        coarse_adjacencies: np.ndarray,
        fine_edges: np.ndarray,
        coarse_edges: np.ndarray,
        fine_edges_centroids: np.ndarray,
        coarse_edges_centroids: np.ndarray,
        fine_faces_of_faces: np.ndarray,
        coarse_boundary_edges: np.ndarray,
        fine_boundary_edges: np.ndarray, 
        dual_id: np.ndarray
    ):

    fine_dual_edge_to_coarse_edge = dict()

    ## first: loop in boundary edges
    for coarse_edge in coarse_boundary_edges:
        coarse_edge_centroid = coarse_edges_centroids[coarse_edge]
        coarse_face_adj = coarse_adjacencies[coarse_edge, 0]
        boundary_fine_edges_in_coarse_face = fine_boundary_edges[
            (primal_id[fine_adjacencies[fine_boundary_edges, 0]]==coarse_face_adj)
        ]
        dists = np.linalg.norm(
            fine_edges_centroids[boundary_fine_edges_in_coarse_face] - coarse_edge_centroid,
            axis=1
        )

        selected_fine_edge = boundary_fine_edges_in_coarse_face[dists <= dists.min()]
        selected_face_to_dual_edge = fine_adjacencies[selected_fine_edge, 0]

        dual_id[selected_face_to_dual_edge] = defnames.dual_ids('edge_id')
        fine_dual_edge_to_coarse_edge.update({selected_face_to_dual_edge[0]: coarse_edge})
    
    # second: loop in internal_edges
    bool_coarse_boundary_edges = np.isin(coarse_edges, coarse_boundary_edges)
    bool_coarse_internal_edges = ~bool_coarse_boundary_edges

    for coarse_edge in coarse_edges[bool_coarse_internal_edges]:
        coarse_edge_centroid = coarse_edges_centroids[coarse_edge]
        coarse_faces_adj = coarse_adjacencies[coarse_edge]
        fine_internal_edges_between_coarse_faces = fine_edges[
            (
                (primal_id[fine_adjacencies[:, 0]] == coarse_faces_adj[0]) &
                (primal_id[fine_adjacencies[:, 1]] == coarse_faces_adj[1])
            ) |
            (
                (primal_id[fine_adjacencies[:, 1]] == coarse_faces_adj[0]) &
                (primal_id[fine_adjacencies[:, 0]] == coarse_faces_adj[1])
            )
        ]
        edges_centroids_between_coarses = fine_edges_centroids[
            fine_internal_edges_between_coarse_faces
        ]

        dists = np.linalg.norm(
            edges_centroids_between_coarses - coarse_edge_centroid,
            axis=1
        )

        selected_faces_to_dual_edge = fine_adjacencies[    
            fine_internal_edges_between_coarse_faces[
                dists <= dists.min()
            ]
        ]

        dual_id[selected_faces_to_dual_edge] = defnames.dual_ids('edge_id')
        for fac in selected_faces_to_dual_edge[0]:
            fine_dual_edge_to_coarse_edge.update({fac: coarse_edge})

    ## three: loop for create paths in coarse faces
    # initial_dual_edges = fine_faces_id[
    #     dual_id==defnames.dual_ids('edge_id')
    # ]

    coarse_face_path = []
    coarse_edge_path = []
    paths = []

    dists = np.zeros(fine_adjacencies.shape, dtype=np.float64)
    dists[:, 0] = np.linalg.norm(
        fine_edges_centroids - fine_faces_centroids[fine_adjacencies[:, 0]],
        axis=1
    )
    dists[:, 1] = np.linalg.norm(
        fine_edges_centroids - fine_faces_centroids[fine_adjacencies[:, 1]],
        axis=1
    )
    dists[fine_adjacencies==-1] = np.inf

    for coarse_face in coarse_faces_id:
        dual_vertice_in_coarse = fine_faces_id[
            (primal_id==coarse_face) & (dual_id==defnames.dual_ids('vertice_id'))
        ]
        dual_edges_in_coarse = fine_faces_id[
            (primal_id==coarse_face) & (dual_id==defnames.dual_ids('edge_id')) 
        ]

        for dual_edge in dual_edges_in_coarse:
            path = get_local_shortest_path_for_create_dual_edges(
                fine_adjacencies,
                dists,
                coarse_face,
                primal_id[fine_adjacencies],
                dual_vertice_in_coarse,
                dual_edge
            )

            dual_id[path] = defnames.dual_ids('edge_id')
            path = np.append(path, [dual_vertice_in_coarse[0], dual_edge])
            coarse_face_path.append(coarse_face)
            coarse_edge_path.append(fine_dual_edge_to_coarse_edge[dual_edge])
            paths.append(path)
    
    coarse_face_path = np.array(coarse_face_path)
    coarse_edge_path = np.array(coarse_edge_path)
    paths = np.array(paths, dtype='O')

    # ## four: check viz of dual_vertices
    # for coarse_face in coarse_faces_id:
    #     dual_vertice = fine_faces_id[
    #         (primal_id==coarse_face) & (dual_id==defnames.dual_ids('vertice_id'))
    #     ]

    #     faces_of_dual_vertice = fine_faces_of_faces[dual_vertice[0]]
    #     test_edges_faces_of_dual_edges = dual_id[faces_of_dual_vertice] == defnames.dual_ids('edge_id')
    #     if test_edges_faces_of_dual_edges.sum() == faces_of_dual_vertice.shape[0]:
    #         pass
    #     else:
    #         faces_not_edges = faces_of_dual_vertice[
    #             dual_id[faces_of_dual_vertice] != defnames.dual_ids('edge_id')
    #         ]
    #         initial_dual_edges_coarse_face = initial_dual_edges[
    #             primal_id[initial_dual_edges] == coarse_face
    #         ]

    #         for face in faces_not_edges:
    #             dists_to_test = np.linalg.norm(
    #                 fine_faces_centroids[initial_dual_edges_coarse_face] - fine_faces_centroids[face],
    #                 axis=1
    #             )
                
    #             selected_dual_initial_face = initial_dual_edges_coarse_face[dists_to_test <= dists_to_test.min()]
    #             path_initial = edges_paths[(dual_vertice[0], selected_dual_initial_face[0])]
    #             dual_id[path_initial] = -1
    #             path = get_local_shortest_path_for_create_dual_edges(
    #                 fine_adjacencies,
    #                 dists,
    #                 coarse_face,
    #                 primal_id[fine_adjacencies],
    #                 [face],
    #                 selected_dual_initial_face[0],
    #                 node_to_remove=dual_vertice[0]
    #             )
    #             dual_id[face] = defnames.dual_ids('edge_id')
    #             dual_id[path] = defnames.dual_ids('edge_id')
    #             edges_paths.update(
    #                 {
    #                     (dual_vertice[0], selected_dual_initial_face[0]): np.concatenate([path, [face]])
    #                 }
    #             )
    
    # all_edges_paths = np.unique(np.concatenate(list(edges_paths.values())))
    # dual_id[all_edges_paths] = defnames.dual_ids('edge_id')
    dual_id[dual_id == -1] = defnames.dual_ids('face_id')
    

    return coarse_face_path, coarse_edge_path, paths



def get_dual_volumes_2d(dual_id: np.ndarray, fine_faces_id: np.ndarray, fine_faces_of_faces: np.ndarray, fine_faces_of_faces_by_faces: np.ndarray):

    dual_faces = fine_faces_id[dual_id == defnames.dual_ids('face_id')]
    dual_others = fine_faces_id[dual_id != defnames.dual_ids('face_id')]
    dual_volumes = []

    while dual_faces.shape[0] > 0:
        face0 = np.array([dual_faces[0]], dtype=np.int)
        test = np.array([True], dtype=bool)
        # dual_volume = [face0]
        while np.any(test):
            faces_of_face0 = np.unique(
                np.union1d(
                    np.concatenate(fine_faces_of_faces[face0]),
                    [face0]
                )
            )
            external_faces = np.setdiff1d(faces_of_face0, face0)
            test = dual_id[external_faces] == defnames.dual_ids('face_id')
            # dual_volume.append(
            #     np.setdiff1d(external_faces, np.concatenate(dual_volume))
            # )
            
            face0 = np.setdiff1d(faces_of_face0, dual_others)
        
        faces_of_face0_v2 = np.unique(
            np.concatenate(
                fine_faces_of_faces[faces_of_face0]
            )
        )
        vertices_to_get = faces_of_face0_v2[
            dual_id[faces_of_face0_v2] == defnames.dual_ids('vertice_id')
        ]
        faces_of_face0 = np.union1d(faces_of_face0, vertices_to_get)
        # dual_volume = np.unique(np.concatenate(dual_volume))
        dual_faces = np.setdiff1d(dual_faces, faces_of_face0)
        dual_volumes.append(faces_of_face0)
    
    dual_volumes = np.array(dual_volumes, dtype='O')
    
    return dual_volumes

def get_dual_volumes_2d_v2(
        dual_id: np.ndarray, 
        fine_faces_id: np.ndarray, 
        fine_faces_of_faces: np.ndarray, 
        fine_faces_of_faces_by_faces: np.ndarray, 
        coarse_nodes, 
        coarse_faces_of_nodes, 
        coarse_edges_of_nodes, 
        coarse_adjacencies, 
        fine_faces_centroids, 
        coarse_nodes_centroids, 
        fine_primal_id,
        coarse_face_path,
        coarse_edge_path,
        edge_paths
    ):

    dual_volumes = []

    for coarse_node in coarse_nodes:
        boundary = []
        coarse_edges_of_node = coarse_edges_of_nodes[coarse_node]

        boundary = np.unique(
            np.concatenate(
                edge_paths[
                    np.isin(coarse_edge_path, coarse_edges_of_node)
                ]
            )
        )

        coarse_node_centroid = coarse_nodes_centroids[coarse_node]
        all_primal_adj = np.unique(coarse_adjacencies[coarse_edges_of_node][
            coarse_adjacencies[coarse_edges_of_node] != -1
        ])
        fine_faces_local = fine_faces_id[np.isin(fine_primal_id, all_primal_adj)]
        dual_faces_local = fine_faces_local[dual_id[fine_faces_local] == defnames.dual_ids('face_id')]
        local_dists = np.linalg.norm(
            fine_faces_centroids[dual_faces_local] - coarse_node_centroid,
            axis=1
        )
        dual_face1 = dual_faces_local[local_dists <= local_dists.min()]

        test = np.array([True])
        while np.any(test):
            faces_of_local_duals = np.unique(np.concatenate(fine_faces_of_faces[dual_face1].flatten()))
            faces_of_local_duals = np.setdiff1d(faces_of_local_duals, dual_face1)
            test = dual_id[faces_of_local_duals] == defnames.dual_ids('face_id')
            dual_face1 = np.append(dual_face1, faces_of_local_duals[test])
        
        dual_volume = np.concatenate([boundary, dual_face1])
        dual_volumes.append(dual_volume)
    
    dual_volumes = np.array(dual_volumes, dtype='O')
    
    return dual_volumes

def get_dual_interaction_region(dual_volumes, coarse_faces_id, fine_faces_id, primal_id, dual_id):
    regions = []
    vertices_selected = np.repeat(-1, coarse_faces_id.shape[0])

    for i, coarse_id in enumerate(coarse_faces_id):
        region = []
        fine_vertice_id = fine_faces_id[
            (dual_id == defnames.dual_ids('vertice_id')) & (primal_id == coarse_id)
        ]
        vertices_selected[i] = fine_vertice_id

        for dual_volume in dual_volumes:
            if np.intersect1d(dual_volume, fine_vertice_id).shape[0] == 1:
                region.append(dual_volume)
           
        region = np.unique(np.concatenate(region))
        regions.append(region)
    
    regions = np.array(regions, dtype='O')
    
    return regions, vertices_selected

def get_dual_interaction_region_v2(
        dual_volumes, 
        coarse_faces_id, 
        fine_faces_id, 
        primal_id, 
        dual_id, 
        coarse_faces_path,
        coarse_edges_path,
        edge_paths,
        coarse_edges_of_faces
    ):
    regions = []
    vertices_selected = np.repeat(-1, coarse_faces_id.shape[0])
    boundarys = []
    internal_paths = []
    initial_ccs = []

    for i, coarse_id in enumerate(coarse_faces_id):
        region = []
        fine_vertice_id = fine_faces_id[
            (dual_id == defnames.dual_ids('vertice_id')) & (primal_id == coarse_id)
        ][0]
        vertices_selected[i] = fine_vertice_id
        

        for dual_volume in dual_volumes:
            if np.intersect1d(dual_volume, fine_vertice_id).shape[0] == 1:
                region.append(dual_volume)
           
        region = np.unique(np.concatenate(region))
        regions.append(region)

        others_vertices = region[dual_id[region] == defnames.dual_ids('vertice_id')]
        others_vertices = np.setdiff1d(others_vertices, fine_vertice_id)
        others_primal_ids = primal_id[others_vertices]

        paths_others_primal_ids = edge_paths[
            np.isin(coarse_faces_path, others_primal_ids)
        ]
        coarse_edges_others_primals = coarse_edges_path[
            np.isin(coarse_faces_path, others_primal_ids)
        ]

        edges_of_coarse_id = coarse_edges_of_faces[coarse_id]
        test = np.isin(coarse_edges_others_primals, edges_of_coarse_id)
        test = ~test
        paths_others_primal_ids = paths_others_primal_ids[test]
        coarse_edges_others_primals = coarse_edges_others_primals[test]
        boundary = np.intersect1d(
            np.unique(np.concatenate(paths_others_primal_ids)),
            region
        )
        boundarys.append(boundary)

        edges_path_of_vertice = np.unique(
            np.concatenate(edge_paths[
            np.isin(coarse_edges_path, edges_of_coarse_id)
        ]))

        internal_paths.append(edges_path_of_vertice)

        first_cc = np.intersect1d(boundary, edges_path_of_vertice)
        initial_ccs.append(first_cc)


    regions = np.array(regions, dtype='O')
    boundarys = np.array(boundarys, dtype='O')
    internal_paths = np.array(internal_paths, dtype='O')
    initial_ccs = np.array(initial_ccs, dtype='O')
    
    return regions, vertices_selected, boundarys, internal_paths, initial_ccs


def create_dual(fine_mesh_properties: MeshProperty, coarse_mesh_properties: MeshProperty, level:int):
    dual_id = np.repeat(-1, fine_mesh_properties['faces'].shape[0])

    create_dual_vertices(
        coarse_faces_id=coarse_mesh_properties['faces'],
        coarse_faces_centroids=coarse_mesh_properties['faces_centroids'],
        fine_faces_id=fine_mesh_properties['faces'],
        fine_faces_centroids=fine_mesh_properties['faces_centroids'],
        primal_id=fine_mesh_properties[defnames.get_primal_id_name_by_level(level)],
        dual_id=dual_id
    )

    coarse_face_path, coarse_edge_path, edge_paths = create_dual_edges_v2(
        coarse_faces_id=coarse_mesh_properties['faces'],
        coarse_faces_centroids=coarse_mesh_properties['faces_centroids'],
        fine_faces_id=fine_mesh_properties['faces'],
        fine_faces_centroids=fine_mesh_properties['faces_centroids'],
        primal_id=fine_mesh_properties[defnames.get_primal_id_name_by_level(level)],
        fine_adjacencies=fine_mesh_properties['adjacencies'],
        coarse_adjacencies=coarse_mesh_properties['adjacencies'],
        fine_edges=fine_mesh_properties['edges'],
        coarse_edges=coarse_mesh_properties['edges'],
        fine_edges_centroids=fine_mesh_properties.edges_centroids,
        coarse_edges_centroids=coarse_mesh_properties.edges_centroids,
        fine_faces_of_faces=fine_mesh_properties.faces_of_faces,
        coarse_boundary_edges=coarse_mesh_properties.boundary_edges,
        fine_boundary_edges=fine_mesh_properties.boundary_edges,
        dual_id=dual_id
    )

    test = dual_id == -1
    if np.any(test):
        raise NotImplementedError
    
    dual_volumes = get_dual_volumes_2d_v2(
        dual_id=dual_id,
        fine_faces_id=fine_mesh_properties['faces'],
        fine_faces_of_faces=fine_mesh_properties.faces_of_faces_by_nodes,
        fine_faces_of_faces_by_faces=fine_mesh_properties.faces_of_faces,
        coarse_nodes=coarse_mesh_properties['nodes'],
        coarse_faces_of_nodes=coarse_mesh_properties['faces_of_nodes'],
        coarse_edges_of_nodes=coarse_mesh_properties['edges_of_nodes'],
        coarse_adjacencies=coarse_mesh_properties['adjacencies'],
        fine_faces_centroids=fine_mesh_properties['faces_centroids'],
        coarse_nodes_centroids=coarse_mesh_properties['nodes_centroids'],
        fine_primal_id=fine_mesh_properties[defnames.get_primal_id_name_by_level(level)],
        coarse_face_path=coarse_face_path,
        coarse_edge_path=coarse_edge_path,
        edge_paths=edge_paths
    )

    regions, vertices_selected, boundarys, internal_paths, initial_ccs = get_dual_interaction_region_v2(
        dual_volumes,
        coarse_mesh_properties['faces'],
        fine_mesh_properties['faces'],
        fine_mesh_properties[defnames.get_primal_id_name_by_level(level)],
        dual_id,
        coarse_faces_path=coarse_face_path,
        coarse_edges_path=coarse_edge_path,
        edge_paths=edge_paths,
        coarse_edges_of_faces=coarse_mesh_properties['edges_of_faces']
    )

    level_str = defnames.level_str(level)
    
    data = {
        defnames.get_dual_id_name_by_level(level): dual_id,
        defnames.get_dual_volumes_name_by_level(level): dual_volumes,
        defnames.get_dual_interation_region_name_by_level(level): regions,
        defnames.vertices_selected + level_str: vertices_selected,
        defnames.boundary_dual_interaction + level_str: boundarys,
        defnames.internal_dual_path + level_str: internal_paths,
        defnames.dual_initial_ccs + level_str: initial_ccs
    }

    return data