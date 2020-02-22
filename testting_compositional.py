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
from update_inputs import FluidProperties, ComponentProperties

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']


M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)

n_volumes = data_impress.len_entities['volumes']

kprop = ComponentProperties(data_loaded)
z, P, T, R, Nc = update_inputs.inputs_overall_properties(data_loaded)
fp = StabilityCheck(z, P, T, R, Nc, kprop)
fluid_properties = FluidProperties(fp, data_loaded, n_volumes)
PropertiesCalc(M, data_impress, wells, fluid_properties, load)
fluid_properties.inputs_missing_properties()

t = 0
tfinal = 1

while t < tfinal:
    # fluid_properties.P = P
    IMPEC( M, data_impress, wells, fluid_properties, fp, kprop, load)
    PropertiesCalc(M, data_impress, wells, fluid_properties, load)

    for i in range(n_volumes):
        P = fluid_properties.P[i]
        z = fluid_properties.z[:,i]
        fp = StabilityCheck(z, P, T, R, Nc, kprop)
        fluid_properties.update_all_volumes(fp, i)

    t = tfinal
    #check stability and perform flash calculation

#data_impress.update_variables_to_mesh()
#M.core.print(file='test'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
