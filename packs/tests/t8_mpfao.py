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
    midpoints_edges
)

pdb.set_trace()

# mpfa_data.init_data(M, elements_lv0, data_impress)
# pdb.set_trace()
