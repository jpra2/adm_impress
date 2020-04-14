import pdb
from packs.running.compositional_initial_mesh_properties import initial_mesh
from packs.compositional.compositionalIMPEC import CompositionalFVM
from packs.directories import data_loaded
from run_compositional import run_simulation
import run_compositional
import scipy.sparse as sp
import numpy as np
from packs.compositional.update_time import delta_time
import matplotlib.pyplot as plt
import update_inputs_compositional

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']


#M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)
#n_volumes = data_impress.len_entities['volumes']
#fprop, fprop_block, kprop = update_inputs_compositional.update(M, data_impress, wells, load, data_loaded, n_volumes)

delta_t_initial = 100
loop_max = 20
t = 0
loop = 0
tmax = 0.01*86400 #422693.9470089848 #seg #0.01* 86400

M, data_impress, wells, fprop, fprop_block, kprop, load, n_volumes = run_compositional.initialize(load, convert)

sim = run_simulation(delta_t_initial, data_impress, fprop)

while t < tmax: # and loop < loop_max:

    sim.run(M, data_impress, wells, fprop, fprop_block, kprop, load, n_volumes)
    t = sim.t
    loop = sim.loop
    if (t + sim.delta_t) > tmax:
        sim.delta_t = tmax - t;
    #print(loop)
    print(sim.t)

sim.save_infos(data_impress, M)

print(fprop.P)

'''
calculate the time of this test to compare results:
 t = 0.157 * porosity * viscosity * ct * L**2 / K #I figure ct stands for rock compressibility or total compressibility
 t = 0.157 * 0.2 * 2.498e(-4) * 7.25e-8 * 609.6**2/5e-13
'''
#data_impress.update_variables_to_mesh()
#M.core.print(file='test'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
