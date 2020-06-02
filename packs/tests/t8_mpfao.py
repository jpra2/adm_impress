from packs.preprocess.preprocess_for_mpfa import PreprocessForMpfa
from packs.running.initial_mesh_properties import initial_mesh
import pdb

M, elements_lv0, data_impress, wells = initial_mesh()

mpfa_data = PreprocessForMpfa()

midpoints_edges = mpfa_data.create_midpoints_edges(
    elements_lv0['edges'],
    elements_lv0['nodes'],
    data_impress['centroid_nodes'],
    elements_lv0['edges_node_nodes']
)

normais, points_subfaces, local_ids_subfaces = mpfa_data.create_subfaces(
    elements_lv0['faces'],
    elements_lv0['edges'],
    elements_lv0['faces_edge_edges'],
    elements_lv0['faces_node_nodes'],
    data_impress['centroid_faces'],
    data_impress['centroid_nodes'],
    midpoints_edges,
    set_ids=True
)

points_subfaces = mpfa_data.get_points_of_subfaces_from_complete_structured_array(
    points_subfaces,
    data_impress['centroid_faces'],
    data_impress['centroid_nodes'],
    midpoints_edges
)

gids_subvolumes, from_volumes_to_gids_subvolumes, structured_points_subvolumes = mpfa_data.create_all_subvolumes_mpfao(
    elements_lv0['volumes'],
    elements_lv0['faces'],
    elements_lv0['volumes_face_faces'],
    elements_lv0['nodes_node_faces'],
    elements_lv0['volumes_face_nodes']
)

pt = mpfa_data.get_points_from_st_subvolume_mpfao(
    structured_points_subvolumes[0],
    data_impress['centroid_volumes'],
    data_impress['centroid_faces']
)

pts = mpfa_data.get_points_from_st_subvolumes_mpfao(
    structured_points_subvolumes[[0,1,2]],
    data_impress['centroid_volumes'],
    data_impress['centroid_faces']
)






pdb.set_trace()

# mpfa_data.init_data(M, elements_lv0, data_impress)
# pdb.set_trace()
