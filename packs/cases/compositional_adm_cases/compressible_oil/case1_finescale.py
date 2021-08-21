import pdb
from packs.directories import data_loaded
from run_compositional import run_simulation
import time
from packs.data_class.compositional_data import CompositionalData
from packs.data_class.compositional_cumulative_datamanager import CumulativeCompositionalDataManager
""" ---------------- LOAD STOP CRITERIA AND MESH DATA ---------------------- """
import numpy as np

description = 'case1_finescale_'
compositional_data = CompositionalData(description=description)
cumulative_compositional_datamanager = CumulativeCompositionalDataManager(description=description)
cumulative_compositional_datamanager.create()
# # datas_comp = cumulative_compositional_datamanager.load_all_datas()
# cumulative_compositional_datamanager.delete_all_datas()
# compositional_data.delete()
# import pdb; pdb.set_trace()

loop_array = np.zeros(1, dtype=[('loop', int), ('t', float)])


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

# import pdb; pdb.set_trace()

while run_criteria < stop_criteria:# and loop < loop_max:
    
    sim.run(M, wells, fprop, load)
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
    
    loop_array['loop'][0] = loop
    loop_array['t'][0] = sim.t
    compositional_data.update({
        'pressure': fprop.P,
        'Sg': fprop.Sg,
        'Sw': fprop.Sw,
        'So': fprop.So,
        'composition': fprop.Csi_j,
        'loop_array': loop_array
    })
    cumulative_compositional_datamanager.insert_data(compositional_data._data)
    
    if loop % 500 == 0:
        compositional_data.export_to_npz()
        cumulative_compositional_datamanager.export()
        import pdb; pdb.set_trace()
        
    


tf = time.time()
print('Total computational time: ', tf-t) #total simulation time
import pdb; pdb.set_trace()
sim.save_infos(data_impress, M) #Save data to file
