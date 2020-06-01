from packs.preprocess.preprocess_for_mpfa import PreprocessForMpfa
from packs.running.initial_mesh_properties import initial_mesh
import pdb

M, elements_lv0, data_impress, wells = initial_mesh()
mpfa_data = PreprocessForMpfa(data_impress)
mpfa_data.init_data(M, elements_lv0, data_impress)
pdb.set_trace()
