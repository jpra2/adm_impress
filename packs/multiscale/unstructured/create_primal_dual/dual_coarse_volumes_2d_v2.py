from packs.manager.meshmanager import MeshProperty
from packs import defnames
import numpy as np
from packs.utils.utils_old import get_local_shortest_path_for_create_dual_edges
from packs.multiscale.unstructured.create_primal_dual.dual_coarse_volumes_2d import get_dual_volumes_2d, get_dual_interaction_region

def create_dual_vertices(
        coarse_faces_id: np.ndarray, 
        coarse_faces_centroids: np.ndarray, 
        fine_faces_id: np.ndarray, 
        fine_faces_centroids: np.ndarray,
        primal_id: np.ndarray, 
        dual_id: np.ndarray,
        coarse_adjacencies: np.ndarray,
        coarse_edges_id: np.ndarray,
        coarse_edges_centroids: np.ndarray,
        coarse_boundary_edges: np.ndarray
    ) -> None:

    dual_vertices = np.repeat(-1, coarse_faces_id.shape[0])

    for coarse_edge in coarse_boundary_edges:
        coarse_face = coarse_adjacencies[coarse_edge, 0]
        coarse_edge_centroid = coarse_edges_centroids[coarse_edge]
        local_fine_faces = fine_faces_id[primal_id == coarse_face]
        dists = np.linalg.norm(
            fine_faces_centroids[local_fine_faces] - coarse_edge_centroid,
            axis=1
        )

        selected_vertice = local_fine_faces[dists <= dists.min()]
        dual_vertices[coarse_face] = selected_vertice

    for coarse_face in np.setdiff1d(coarse_faces_id, coarse_adjacencies[coarse_boundary_edges, 0]):
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

    fine_boundary_faces = fine_adjacencies[fine_boundary_edges, 0]

    ## first: loop in boundary faces
    # for coarse_edge in coarse_boundary_edges:
    #     coarse_edge_centroid = coarse_edges_centroids[coarse_edge]
    #     coarse_face_adj = coarse_adjacencies[coarse_edge, 0]
    #     boundary_fine_edges_in_coarse_face = fine_boundary_edges[
    #         (primal_id[fine_adjacencies[fine_boundary_edges, 0]]==coarse_face_adj)
    #     ]
    #     dists = np.linalg.norm(
    #         fine_edges_centroids[boundary_fine_edges_in_coarse_face] - coarse_edge_centroid,
    #         axis=1
    #     )

    #     selected_fine_edge = boundary_fine_edges_in_coarse_face[dists <= dists.min()]
    #     selected_face_to_dual_edge = fine_adjacencies[selected_fine_edge, 0]

    #     dual_id[selected_face_to_dual_edge] = defnames.dual_ids('edge_id')
    
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
    
    boundary_coarse_faces = coarse_adjacencies[coarse_adjacencies[:, 1]==-1, 0]

    for coarse_face in boundary_coarse_faces:
        dual_vertice_in_coarse = fine_faces_id[
            (primal_id==coarse_face) & (dual_id==defnames.dual_ids('vertice_id'))
        ]
        dual_edges_in_coarse = fine_faces_id[
            (primal_id==coarse_face) & (dual_id==defnames.dual_ids('edge_id')) 
        ]

        for dual_edge in dual_edges_in_coarse:
            faces_adj_vertice = fine_faces_of_faces[dual_vertice_in_coarse[0]]
            dists1 = np.linalg.norm(
                fine_faces_centroids[faces_adj_vertice] - fine_faces_centroids[dual_edge],
                axis=1
            )
            selected_edge1 = faces_adj_vertice[dists1 <= dists1.min()]

            faces_of_selected_edge1 = fine_faces_of_faces[selected_edge1[0]]
            faces_of_selected_edge1 = np.setdiff1d(faces_of_selected_edge1, fine_boundary_faces)
            dists2 = np.linalg.norm(
                fine_faces_centroids[dual_edge] - fine_faces_centroids[faces_of_selected_edge1],
                axis=1
            )
            selected_edge2 = faces_of_selected_edge1[dists2 <= dists2.min()]

            path = get_local_shortest_path_for_create_dual_edges(
                fine_adjacencies,
                dists,
                coarse_face,
                primal_id[fine_adjacencies],
                selected_edge2,
                dual_edge
            )

            dual_id[path] = defnames.dual_ids('edge_id')
            dual_id[selected_edge1] = defnames.dual_ids('edge_id')
            dual_id[selected_edge2] = defnames.dual_ids('edge_id')

    for coarse_face in np.setdiff1d(coarse_faces_id, boundary_coarse_faces):
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
    
    dual_id[dual_id == -1] = defnames.dual_ids('face_id')

def create_dual(fine_mesh_properties: MeshProperty, coarse_mesh_properties: MeshProperty, level:int):
    
    dual_id = np.repeat(-1, fine_mesh_properties['faces'].shape[0])

    create_dual_vertices(
        coarse_faces_id=coarse_mesh_properties['faces'],
        coarse_faces_centroids=coarse_mesh_properties['faces_centroids'],
        fine_faces_id=fine_mesh_properties['faces'],
        fine_faces_centroids=fine_mesh_properties['faces_centroids'],
        primal_id=fine_mesh_properties[defnames.get_primal_id_name_by_level(level)],
        dual_id=dual_id,
        coarse_adjacencies=coarse_mesh_properties['adjacencies'],
        coarse_edges_id=coarse_mesh_properties['edges'],
        coarse_edges_centroids=coarse_mesh_properties.edges_centroids,
        coarse_boundary_edges=coarse_mesh_properties.boundary_edges
    )

    create_dual_edges(
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
        fine_mesh_properties[defnames.get_primal_id_name_by_level(level)],
        dual_id
    )
    
    data = {
        defnames.get_dual_id_name_by_level(level): dual_id,
        defnames.get_dual_volumes_name_by_level(level): dual_volumes,
        defnames.get_dual_interation_region_name_by_level(level): regions
    }

    return data