from packs.running.initial_mesh_properties import initial_mesh
from packs.adm.adm_method import AdmMethod
from packs.biphasic.biphasic_ms.biphasic_multiscale import BiphasicTpfaMultiscale
from packs.directories import data_loaded
from testting1_monophasic_multilevel import(M, elements_lv0, wells, data_impress,
                                            multilevel_operators as mlo)
import scipy.sparse as sp
import numpy as np
import time
import pdb

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']
n = data_loaded['n_test']
biphasic = data_loaded['biphasic']
n_levels = data_loaded['n_levels']

# M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)

if biphasic:
    b1 = BiphasicTpfaMultiscale(M, data_impress, elements_lv0, wells)

adm_method = AdmMethod(wells['all_wells'], 2, M, data_impress, elements_lv0)
# adm_method.restart_levels()
# adm_method.set_level_wells()
# adm_method.set_adm_mesh()
T, b = b1.get_T_and_b()
# import pdb; pdb.set_trace()
# adm_method.set_initial_mesh(mlo, T, b)
####################
adm_method.restart_levels()
adm_method.set_level_wells()
adm_method.set_adm_mesh()
#########################

# import pdb; pdb.set_trace()
verif = True


while verif:

    for level in range(1, n_levels):

        adm_method.organize_ops_adm(mlo['prolongation_level_'+str(level)],
                                    mlo['restriction_level_'+str(level)],
                                    level)

    adm_method.solve_multiscale_pressure(T, b)
    adm_method.set_pcorr()
    pdb.set_trace()
    b1.get_velocity_faces()
    b1.get_flux_volumes()
    pdb.set_trace()
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
