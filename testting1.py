# from packs.solvers.solvers_scipy.solver_sp import SolverSp
# from packs.running.run_simulation import RunSimulation
# # from packs.direct_solution.monophasic.monophasic1 import Monophasic
# from packs.direct_solution.biphasic.biphasic1 import Biphasic
# from packs.preprocess.preprocess1 import Preprocess1
import pdb
#
# rodar = RunSimulation(state=5)
# M = rodar.M
#
#
# prep1 = Preprocess1()
# prep1.set_saturation_regions(M)
#
# # m1 = Monophasic(M)
# m1 = Biphasic(M)
# m1.get_transmissibility_matrix_without_contours()
# m1.get_transmissibility_matrix()
# m1.get_RHS_term()
#
#
# solver = SolverSp()
# x = solver.direct_solver(m1.datas['T'], m1.datas['b'])
# m1.get_solution(x)
# m1.get_flux_faces_and_volumes()
#
# pdb.set_trace()
# M.data.update_variables_to_mesh()

#################################
##test Monophasic
# from packs.type_simulation.monophasic_tpfa import monophasicTpfa
# m1 = monophasicTpfa(M)
# m1.run()
# import pdb; pdb.set_trace()
#################################

#############################
# # test biphasic
# from packs.type_simulation.biphasic_simulation.biphasic_tpfa import biphasicTpfa
# import time
# b1 = biphasicTpfa(M, load=False)
# b1.run()
# b1.run()
# b1.update_flux_w_and_o_volumes()
# b1.update_delta_t()
# b1.update_saturation()
# b1.update_relative_permeability()
# b1.update_mobilities()
# b1.update_transmissibility()

# def r1():
#     b1.run()
#
# def r2():
#     r1()
#     b1.mesh.data.update_variables_to_mesh()
#     b1.mesh.core.print(file='test', extension='.vtk', config_input="input_cards/print_settings0.yml")
#
#
# def mostrar():
#     b1.mesh.data.update_variables_to_mesh()
#     b1.mesh.core.print(file='test', extension='.vtk', config_input="input_cards/print_settings0.yml")
#
# verif = True
# contador = 1
# while verif:
#     if contador % 2 == 0:
#         contador = 1
#         import pdb; pdb.set_trace()
#     r1()
#     contador += 1





from packs.load.preprocessor0 import M
from packs.data_class.data_impress import Data
from packs.data_class.elements_lv0 import ElementsLv0
from packs.contours.wells import Wells
from packs.preprocess.prep0_0 import Preprocess0

# import pdb; pdb.set_trace()
elements_lv0 = ElementsLv0(M)
data_impress = Data(M, elements_lv0)
Preprocess0(M, elements_lv0)
wells = Wells(M)
wells.export_all_datas_to_npz()
data_impress.update_variables_to_mesh()
import pdb; pdb.set_trace()










import pdb; pdb.set_trace()







# import pickle
#
# from tcc.load_save_initialize.load_infos import LoadInfos
# from tcc.dual_mesh.create_dual_mesh import DualMesh1
# import numpy as np
# from . import directories
# import pdb; pdb.set_trace()
#
# # file_name = os.path.join(path_flying, 'mesh_obj.txt')
# # with open(file_name, 'wb') as handle:
# #     pickle.dump(M, handle)
#
#
#
# LoadInfos(M)
# DualMesh1(M)
#
# import pdb; pdb.set_trace()
