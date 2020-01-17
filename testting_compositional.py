import pdb
from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional import CompositionalIMPEC
from packs.directories import data_loaded
import scipy.sparse as sp
import numpy as np
import time

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']

M, elements_lv0, data_impress, wells, fluid_properties = initial_mesh(load=load, convert=convert)
tpfa_solver = CompositionalIMPEC.update_transmissibility(M, data_impress, elements_lv0, fluid_properties)
#data_impress.update_variables_to_mesh()
#M.core.print(file='test'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
