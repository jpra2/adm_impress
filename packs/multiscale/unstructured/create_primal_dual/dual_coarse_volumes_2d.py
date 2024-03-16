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

def get_dual_interaction_region(dual_volumes, coarse_faces_id, fine_faces_id, primal_id, dual_id):
    regions = []
    for coarse_id in coarse_faces_id:
        region = []
        fine_vertice_id = fine_faces_id[
            (dual_id == defnames.dual_ids('vertice_id')) & (primal_id == coarse_id)
        ]

        for dual_volume in dual_volumes:
            if np.intersect1d(dual_volume, fine_vertice_id).shape[0] == 1:
                region.append(dual_volume)
        
        region = np.unique(np.concatenate(region))

        regions.append(region)
    
    regions = np.array(regions, dtype='O')
    return regions


def create_dual(fine_mesh_properties: MeshProperty, coarse_mesh_properties: MeshProperty, level:int):
    dual_id = np.repeat(-1, fine_mesh_properties['faces'].shape[0])

    create_dual_vertices(
        coarse_faces_id=coarse_mesh_properties['faces'],
        coarse_faces_centroids=coarse_mesh_properties['faces_centroids'],
        fine_faces_id=fine_mesh_properties['faces'],
        fine_faces_centroids=fine_mesh_properties['faces_centroids'],
        primal_id=fine_mesh_properties[defnames.fine_primal_id],
        dual_id=dual_id
    )

    create_dual_edges(
        coarse_faces_id=coarse_mesh_properties['faces'],
        coarse_faces_centroids=coarse_mesh_properties['faces_centroids'],
        fine_faces_id=fine_mesh_properties['faces'],
        fine_faces_centroids=fine_mesh_properties['faces_centroids'],
        primal_id=fine_mesh_properties[defnames.fine_primal_id],
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
    
    dual_volumes = get_dual_volumes_2d(
        dual_id=dual_id,
        fine_faces_id=fine_mesh_properties['faces'],
        fine_faces_of_faces=fine_mesh_properties.faces_of_faces_by_nodes,
        fine_faces_of_faces_by_faces=fine_mesh_properties.faces_of_faces
    )

    regions = get_dual_interaction_region(
        dual_volumes,
        coarse_mesh_properties['faces'],
        fine_mesh_properties['faces'],
        fine_mesh_properties[defnames.fine_primal_id],
        dual_id
    )
    
    data = {
        defnames.get_dual_id_name_by_level(level): dual_id,
        defnames.get_dual_volumes_name_by_level(level): dual_volumes,
        defnames.get_dual_interation_region_name_by_level(level): regions
    }

    return data