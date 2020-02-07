from packs.running.initial_mesh_properties import initial_mesh
from packs.pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from packs.directories import data_loaded

M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)
adm_method = AdmMethod(wells['all_wells'], 2, M, data_impress, elements_lv0)
adm_method.so_nv1 = True

adm_method.restart_levels()
adm_method.set_level_wells()
adm_method.set_adm_mesh()

import pdb; pdb.set_trace()
