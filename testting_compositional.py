import pdb
from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional.CompositionalIMPEC import CompositionalTPFA
from packs.directories import data_loaded
from packs.compositional.stability_check import StabilityCheck
import scipy.sparse as sp
import numpy as np
import time

import update_inputs

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']
w, Bin, R, Tc, Pc, Vc, T, P, Mw, C7, z = update_inputs.inputs_components_properties(data_loaded['compositional_data'])

M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)
fluid_properties = StabilityCheck(w, Bin, R, Tc, Pc, Vc, T, P, Mw, C7)
fluid_properties.run(z)
t = 0
tfinal = 1
while t < tfinal:
    CompositionalTPFA( M, data_impress, wells, fluid_properties, elements_lv0, load)
    t = tfinal
    #check stability and perform flash calculation

#data_impress.update_variables_to_mesh()
#M.core.print(file='test'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
