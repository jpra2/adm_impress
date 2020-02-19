import pdb
from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional.compositionalIMPEC import IMPEC
from packs.directories import data_loaded
from packs.compositional.stability_check import StabilityCheck
from packs.compositional.properties_calculation import PropertiesCalc
import scipy.sparse as sp
import numpy as np
import time

import update_inputs

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']


M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)

n_volumes = data_impress.len_entities['volumes']

w, Bin, R, Tc, Pc, Vc, T, P, Mw, C7, z = update_inputs.inputs_components_properties(data_loaded, n_volumes)

fluid_properties = StabilityCheck(w, Bin, R, Tc, Pc, Vc, T, P[0], Mw, C7, z)
update_inputs.inputs_water_properties(data_loaded, fluid_properties)
PropertiesCalc(M, data_impress, wells, fluid_properties, elements_lv0, load)

fluid_properties.P = P
fluid_properties.component_phase_mole_numbers = fluid_properties.component_molar_fractions * fluid_properties.phase_mole_numbers
fluid_properties.component_mole_numbers = np.sum(fluid_properties.component_phase_mole_numbers, axis = 1)

t = 0
tfinal = 1

while t < tfinal:
    fluid_properties.P = P
    IMPEC( M, data_impress, wells, fluid_properties, elements_lv0, load)
    PropertiesCalc(M, data_impress, wells, fluid_properties, elements_lv0, load)
    t = tfinal
    #check stability and perform flash calculation

#data_impress.update_variables_to_mesh()
#M.core.print(file='test'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
