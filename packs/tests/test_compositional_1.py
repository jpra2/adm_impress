from packs.multiscale.preprocess.dual_primal.create_dual_and_primal_mesh import MultilevelData

import pdb
from packs.directories import data_loaded
from run_compositional import run_simulation
import time

""" ---------------- LOAD STOP CRITERIA AND MESH DATA ---------------------- """

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

multilevel_structure = MultilevelData(data_impress, M)

import pdb; pdb.set_trace()