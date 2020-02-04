from packs.running.initial_mesh_properties import initial_mesh
from packs.adm.adm_method import AdmMethod
from packs.biphasic.biphasic_tpfa import BiphasicTpfa
from packs.directories import data_loaded
from testting1_monophasic_multilevel import(M, elements_lv0, wells, data_impress,
                                            multilevel_operators as mlo, tpfa_solver)
import scipy.sparse as sp
import numpy as np
import time


adm_method = AdmMethod(wells['all_wells'], 2, M, data_impress, elements_lv0)
T, b = tpfa_solver.run()
adm_method.set_initial_mesh(mlo, T, b)

adm_method.organize_ops_adm(mlo['prolongation_level_1'],
                            mlo['restriction_level_1'],
                            1)

adm_method.organize_ops_adm(mlo['prolongation_level_2'],
                            mlo['restriction_level_2'],
                            2)

adm_method.solve_multiscale_pressure(T, b)
adm_method.set_pcorr()

data_impress.update_variables_to_mesh()
M.core.print(folder='results', file='test_'+ str(0), extension='.vtk', config_input='input_cards/print_settings0.yml')
