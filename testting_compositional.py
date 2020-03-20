import pdb
from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional.compositionalIMPEC import CompositionalFVM
from packs.directories import data_loaded
from packs.compositional.simulation import run_simulation
import scipy.sparse as sp
import numpy as np
import time
from packs.compositional.update_time import delta_time
import update_inputs_compositional

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']


M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)
n_volumes = data_impress.len_entities['volumes']
fprop, fprop_block, kprop = update_inputs_compositional.update(M, data_impress, wells, load, data_loaded, n_volumes)

sim = run_simulation(fprop, data_impress)
loop_max = 50
t = 0
loop = 0
tfinal = 1

while t < tfinal and loop < loop_max:

    sim.run(M, data_impress, wells, fprop, fprop_block, kprop, load, n_volumes)
    t = sim.t
    loop = sim.loop

'''
calculate the time of this test to compare results:
 t = 0.157 * porosity * viscosity * ct * L**2 / K #I figure ct stands for rock compressibility or total compressibility
'''
#data_impress.update_variables_to_mesh()
#M.core.print(file='test'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
