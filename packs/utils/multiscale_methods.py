import numpy as np
from packs.manager import MeshData, MeshProperty
from packs.defnames import get_primal_id_name_by_level

def get_interfaces_edges(
        fine_edges,
        primal_id,
        fine_adjacencies,
        coarse_ids_to_show
    ):

    bool_all_edges_in_intersection_and_boundary = primal_id[fine_adjacencies[:, 0]] != primal_id[fine_adjacencies[:, 1]]
    # import pdb; pdb.set_trace()

    edges_to_print = []
    for coarse_id in coarse_ids_to_show:
        edges_intersection = fine_edges[
            (
                (primal_id[fine_adjacencies[:, 0]] == coarse_id) |
                (primal_id[fine_adjacencies[:, 1]] == coarse_id)
            ) &
            bool_all_edges_in_intersection_and_boundary
        ]
        edges_to_print.append(edges_intersection)
    
    edges_to_print = np.unique(np.concatenate(edges_to_print))
    
    return edges_to_print











def print_fine_interfaces_coarse_mesh_2d(
        fine_mesh_properties: MeshProperty,
        fine_mesh_path: str,
        level: int,
        export_name: str
    ):

    edges_to_print = get_interfaces_edges(
        fine_mesh_properties['edges'],
        fine_mesh_properties[get_primal_id_name_by_level(level)],
        fine_mesh_properties['adjacencies'],
        np.unique(fine_mesh_properties[get_primal_id_name_by_level(level)])
    )

    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.export_only_the_elements(
        export_name,
        'edges',
        edges_to_print
    )



