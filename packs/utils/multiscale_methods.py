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
    
    if len(edges_to_print) > 0:
        pass
    else:
        return np.array(edges_to_print)
    
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

    edges_to_print = np.unique(np.concatenate([edges_to_print, fine_mesh_properties.boundary_edges]))

    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.export_only_the_elements(
        export_name,
        'edges',
        edges_to_print
    )

def print_adm_interfaces_2d(
        fine_mesh_properties: MeshProperty,
        fine_mesh_path: str,
        fine_levels: np.ndarray,
        export_name: str
):
    levels = np.setdiff1d(np.unique(fine_levels), [0])

    edges_to_print = get_interfaces_edges(
        fine_mesh_properties['edges'],
        fine_mesh_properties['faces'],
        fine_mesh_properties['adjacencies'],
        fine_mesh_properties['faces'][fine_levels == 0]
    )
    
    all_edges = [edges_to_print]
    for level in levels:
        edges_to_print = get_interfaces_edges(
            fine_mesh_properties['edges'],
            fine_mesh_properties[get_primal_id_name_by_level(level)],
            fine_mesh_properties['adjacencies'],
            np.unique(fine_mesh_properties[get_primal_id_name_by_level(level)][fine_levels==level])
        )
        all_edges.append(edges_to_print)
    
    all_edges = np.unique(np.concatenate(all_edges))
    all_edges = np.unique(np.concatenate([
        all_edges,
        fine_mesh_properties.boundary_edges
    ]))

    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.export_only_the_elements(export_name, element_type='edges', elements_array=all_edges)

