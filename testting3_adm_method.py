from packs.running.initial_mesh_properties import initial_mesh
from packs.adm.adm_method import AdmMethod
from packs.directories import data_loaded
import scipy.sparse as sp
import numpy as np
import time

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']
n = data_loaded['n_test']

M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)
adm_method = AdmMethod(wells['all_wells'], 2, M, data_impress)
adm_method.set_level_wells()
adm_method.set_adm_mesh()
import pdb; pdb.set_trace()

data_impress.update_variables_to_mesh()
M.core.print(folder='results',file='test'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
