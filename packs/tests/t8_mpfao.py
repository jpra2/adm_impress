from packs.preprocess.preprocess_for_mpfa import PreprocessForMpfa
from packs.running.initial_mesh_properties import initial_mesh
import pdb

M, elements_lv0, data_impress, wells = initial_mesh()
elements_lv0.create_adj_matrix_volumes_to_faces()
elements_lv0.create_adj_matrix_faces_to_edges()
elements_lv0.create_adj_matrix_edges_to_nodes()

f1 = elements_lv0.volumes_to_nodes([0,1,2,3])

# fs0 = elements_lv0.get_adj_volumes_to_faces([0, 1])

# edges_v012 = elements_lv0.get_adj_edges_to_volumes([0,1,2])


mpfa_data = PreprocessForMpfa()

pdb.set_trace()
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

gids_subvolumes, from_volumes_to_gids_subvolumes, str_points_subvolumes = mpfa_data.create_all_subvolumes_mpfao(
    elements_lv0['volumes'],
    elements_lv0['faces'],
    elements_lv0['volumes_face_faces'],
    elements_lv0['nodes_node_faces'],
    elements_lv0['volumes_face_nodes']
)



pdb.set_trace()

# mpfa_data.init_data(M, elements_lv0, data_impress)
# pdb.set_trace()
