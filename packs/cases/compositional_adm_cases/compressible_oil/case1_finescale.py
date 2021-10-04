import pdb
from packs.directories import data_loaded
from run_compositional import run_simulation
import time
from packs.data_class.compositional_data import CompositionalData
from packs.data_class.compositional_cumulative_datamanager import CumulativeCompositionalDataManager
""" ---------------- LOAD STOP CRITERIA AND MESH DATA ---------------------- """
import numpy as np
from packs.cases.compositional_adm_cases.compressible_oil import all_functions


# description = 'case3_finescale_3k'
description = 'case17_finescale_6k_5000_'
compositional_data = CompositionalData(description=description)
cumulative_compositional_datamanager = CumulativeCompositionalDataManager(description=description)
cumulative_compositional_datamanager.create()
# # datas_comp = cumulative_compositional_datamanager.load_all_datas()
# cumulative_compositional_datamanager.delete_all_datas()
# compositional_data.delete()
# import pdb; pdb.set_trace()

loop_array = all_functions.get_empty_loop_array()


name_current = 'current_compositional_results_'
name_all = data_loaded['name_save_file'] + '_'
mesh = 'mesh/' + data_loaded['mesh_name']

if data_loaded['use_vpi']:
    stop_criteria = max(data_loaded['compositional_data']['vpis_para_gravar_vtk'])
else: stop_criteria = data_loaded['compositional_data']['maximum_time']

loop_max = 1000
run_criteria = 0
loop = 0
""" ----------------------------- RUN CODE --------------------------------- """

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']

t = time.time()
sim = run_simulation(name_current, name_all)
M, data_impress, wells, fprop, load = sim.initialize(load, convert, mesh)

# #################
# ## update data before init
# from packs.cases.compositional_adm_cases.compressible_oil import update_variables_before_init as update_var

# update_var.update_variables_for_initial_run_finescale(fprop, sim, compositional_data)
# ########################

# import pdb; pdb.set_trace()

while run_criteria < stop_criteria:# and loop < loop_max:
    # import pdb; pdb.set_trace()
    
    t0 = time.time()
    sim.run(M, wells, fprop, load)
    simulation_time = time.time() - t0
    
    if data_loaded['use_vpi']:
        'If using time-step unit as vpi'
        run_criteria = sim.vpi
    else:
        'If using time-step unit as second'
        run_criteria = sim.t
        if sim.time_save[-1] == 0.0:
            'if time_save = [0.0] only, it means that all time-steps are going \
            to be saved'
            t_next = sim.t + sim.delta_t
        else:
            'This is made so time to save is equal to the simulation time - this \
            only works (for now) if time to save is in seconds and not vpi'
            t_next = sim.time_save[sim.time_save > sim.t]
            if len(t_next)>=1: t_next = t_next[0]

        'If the current simulation time plus the computed time-step is bigger \
        than the final simulation time, correct the time-step so the current \
        simulation time plus delta_t is equal to the final time'
        if sim.t + sim.delta_t > t_next:
            sim.delta_t = t_next - sim.t
    loop = sim.loop
    print(sim.t)
    
    loop_array['total_simulation_time'][0] += simulation_time
    loop_array['n_total_loops'][0] += 1
    
    loop_array['loop'][0] = loop
    loop_array['t'][0] = sim.t
    loop_array['vpi'][0] = sim.vpi
    loop_array['simulation_time'][0] = simulation_time
    loop_array['oil_production'][0] = sim.oil_production
    loop_array['gas_production'][0] = sim.gas_production
    compositional_data.update({
        'pressure': fprop.P,
        'Sg': fprop.Sg,
        'Sw': fprop.Sw,
        'So': fprop.So,
        'q': fprop.q,
        'global_composition': fprop.z,
        'mols': fprop.Nk,
        'xkj': fprop.xkj,
        'Vp': fprop.Vp,
        'loop_array': loop_array
    })
    cumulative_compositional_datamanager.insert_data(compositional_data._data)
    
    if loop % 500 == 0:
        compositional_data.export_to_npz()
        cumulative_compositional_datamanager.export()
        # import pdb; pdb.set_trace()
    
    # import pdb; pdb.set_trace()
        
    
# loop_array['loop'][0] = loop
# loop_array['t'][0] = sim.t
# loop_array['vpi'][0] = sim.vpi
# loop_array['simulation_time'][0] = simulation_time
# loop_array['oil_production'][0] = sim.oil_production
# loop_array['gas_production'][0] = sim.gas_production
# compositional_data.update({
#     'pressure': fprop.P,
#     'Sg': fprop.Sg,
#     'Sw': fprop.Sw,
#     'So': fprop.So,
#     'global_composition': fprop.z,
#     'mols': fprop.Nk,
#     'xkj': fprop.xkj,
#     'Vp': fprop.Vp,
#     'loop_array': loop_array
# })
# cumulative_compositional_datamanager.insert_data(compositional_data._data)
compositional_data.export_to_npz()
cumulative_compositional_datamanager.export()

tf = time.time()
print('Total computational time: ', tf-t) #total simulation time
import pdb; pdb.set_trace()
sim.save_infos(data_impress, M) #Save data to file
