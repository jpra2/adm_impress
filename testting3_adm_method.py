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
# adm_method.restart_levels()
# adm_method.set_level_wells()
# adm_method.set_adm_mesh()
T, b = b1.get_T_and_b()
adm_method.set_initial_mesh(mlo, T, b)

import pdb; pdb.set_trace()
verif = True


while verif:

    adm_method.organize_ops_adm(mlo['prolongation_level_1'],
                                mlo['restriction_level_1'],
                                1)

    adm_method.organize_ops_adm(mlo['prolongation_level_2'],
                                mlo['restriction_level_2'],
                                2)
    import pdb; pdb.set_trace()

    adm_method.solve_multiscale_pressure(T, b)
    adm_method.set_pcorr()
    # adm_method.set_paralel_pcorr()
    # p2 = adm_method.solver.direct_solver(T, b)
    # data_impress['pressure'] = p2
    # data_impress['erro'] = np.absolute((p2-data_impress['pms']))

    b1.run_2()

    adm_method.restart_levels_2()
    # adm_method.set_level_wells()
    adm_method.set_saturation_level()
    adm_method.set_adm_mesh()

    # n=0
    # data_impress.update_variables_to_mesh()
    # M.core.print(folder='results', file='test_'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')

    T, b = b1.get_T_and_b()

    import pdb; pdb.set_trace()



M.core.print(folder='results', file='test_'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')
import pdb; pdb.set_trace()
