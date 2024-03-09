import numpy as np
from packs.manager import MeshData, MeshProperty

def get_interfaces(
        fine_edges,
        primal_id,
        fine_adjacencies,
        coarse_ids_to_show
    ):

    bool_all_edges_in_intersection_and_boundary = primal_id[fine_adjacencies[:, 0]] != primal_id[fine_adjacencies[:, 1]]

    edges_to_print = []
    for coarse_id in coarse_ids_to_show:
        edges_intersection = fine_edges[
            (primal_id[fine_adjacencies[:, 0]] == coarse_id) |
            (primal_id[fine_adjacencies[:, 1]] == coarse_id) &
            bool_all_edges_in_intersection_and_boundary
        ]
        edges_to_print.append(edges_intersection)
    
    return np.unique(edges_to_print)











def print_fine_interfaces_coarse_mesh(
        fine_mesh_properties: MeshProperty,
        fine_mesh_path: str
    ):
    mesh_data = MeshData(fine_mesh_path)
