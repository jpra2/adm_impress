from packs.running.initial_mesh_properties import initial_mesh
from packs.pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from packs.directories import data_loaded

load = data_loaded['load_mesh']
convert = data_loaded['convert_english_to_SI']
n = data_loaded['n_test']

M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)
tpfa_solver = FineScaleTpfaPressureSolver(data_impress, elements_lv0, wells)
tpfa_solver.run()
data_impress.update_variables_to_mesh()
M.core.print(file='test'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
