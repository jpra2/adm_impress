from packs.running.initial_mesh_properties import initial_mesh
from packs.adm.adm_method import AdmMethod
from packs.biphasic.biphasic_tpfa import BiphasicTpfa
from packs.directories import data_loaded
from testting1_monophasic_multilevel import(M, elements_lv0, wells, data_impress,
                                            multilevel_operators as mlo, tpfa_solver)
import scipy.sparse as sp
import numpy as np
import time

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']
n = data_loaded['n_test']
biphasic = data_loaded['biphasic']

# M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)

if biphasic:
    b1 = BiphasicTpfa(M, data_impress, elements_lv0, wells)

adm_method = AdmMethod(wells['all_wells'], 2, M, data_impress, elements_lv0)
adm_method.restart_levels()
adm_method.set_level_wells()
adm_method.set_adm_mesh()
adm_method.organize_ops_adm(mlo['prolongation_level_1'],
                            mlo['restriction_level_1'],
                            1)

adm_method.organize_ops_adm(mlo['prolongation_level_2'],
                            mlo['restriction_level_2'],
                            2)

T, b = tpfa_solver.run()
adm_method.solve_multiscale_pressure(T, b)
adm_method.set_pms_flux_intersect_faces()
t1 = time.time()
# adm_method.set_pcorr()
adm_method.set_paralel_pcorr()
t2 = time.time()
print(t2-t1)
import pdb; pdb.set_trace()
# b1.run_2()
# adm_method.set_pms_flux_volumes()

# import pdb; pdb.set_trace()
# p2 = adm_method.solver.direct_solver(T, b)
# data_impress['pressure'] = p2
# tpfa_solver.get_flux_faces_and_volumes()

data_impress.update_variables_to_mesh()
if biphasic:
    n=1
import pdb; pdb.set_trace()
M.core.print(folder='results', file='test'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
