import pdb
from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional.compositionalIMPEC import CompositionalFVM
from packs.directories import data_loaded
from packs.compositional.stability_check import StabilityCheck
from packs.compositional.properties_calculation import PropertiesCalc
import scipy.sparse as sp
import numpy as np
import time

from packs.compositional.update_time import delta_time
import update_inputs

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']


M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)
n_volumes = data_impress.len_entities['volumes']
fprop, fprop_block, kprop = update_inputs.update(M, data_impress, wells, load, data_loaded, n_volumes)

t = 0
tfinal = 1
deltaT = 0.02
while t < tfinal:
    t_obj = delta_time(fprop) #get wanted properties in t=n
    CompositionalFVM(M, data_impress, wells, fprop, fprop_block, kprop, load, deltaT)
    prop = PropertiesCalc(M, data_impress, wells, fprop, load)
    prop.run_inside_loop(data_impress, wells, fprop)
    for i in range(1,n_volumes):
        P = fprop.P[i]
        z = fprop.z[0:fprop.Nc,i] #água não entra
        fprop_block = StabilityCheck(z, P, fprop.T, fprop.R, fprop.Nc, kprop)
        fprop.update_all_volumes(fprop_block, i)
    import pdb; pdb.set_trace()
    deltaT = t_obj.update_deltaT(deltaT, fprop)#get deltaT with properties in t=n and t=n+1
    t = t + deltaT
    #check stability and perform flash calculation

#data_impress.update_variables_to_mesh()
#M.core.print(file='test'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
