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






# from packs.data_class.data_impress import Data
# from packs.data_class.elements_lv0 import ElementsLv0
# from packs.contours.wells import Wells
# from packs.convert_unit.conversion import Conversion

# import pdb; pdb.set_trace()

# def initial_mesh(load=False, convert=False):
#     from packs.load.preprocessor0 import M
#
#     elements_lv0 = ElementsLv0(M, load=load)
#     data_impress = Data(M, elements_lv0, load=load)
#     wells = Wells(M, load=load)
#     if convert:
#         conversion = Conversion(wells, data_impress)
#         conversion.convert_English_to_SI()
#
#     if not load:
#
#         wells.update_values_to_mesh()
#         wells.export_to_npz()
#         data_impress.update_variables_to_mesh()
#         data_impress.export_to_npz()
#
#     return M, elements_lv0, data_impress, wells

from packs.running.initial_mesh_properties import initial_mesh
from packs.pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.directories import data_loaded
import scipy.sparse as sp
import numpy as np
import time

def get_gids_and_primal_id(gids, primal_ids):

    gids2 = np.unique(gids)
    primal_ids2 = []
    for i in gids2:
        primal_id = np.unique(primal_ids[gids==i])
        if len(primal_id) > 1:
            raise ValueError('erro get_gids_and_primal_id')
        primal_ids2.append(primal_id[0])

    primal_ids2 = np.array(primal_ids2)

    return gids2, primal_ids2

def mostrar(i, data_impress, M, op1, rest1):
    l0 = np.concatenate(op1[:,i].toarray())
    el0 = np.concatenate(rest1[i].toarray())
    data_impress['verif_po'] = l0
    data_impress['verif_rest'] = el0
    data_impress.update_variables_to_mesh(['verif_po', 'verif_rest'])
    M.core.print(file='results/test_'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')

def mostrar_2(i, data_impress, M, op, rest, gid0, gid_coarse1, gid_coarse2):
    l0 = np.concatenate(op[:,i].toarray())
    el0 = np.concatenate(rest[i].toarray())
    el2 = np.zeros(len(gid0))
    l2 = el2.copy()
    for fid in el0:
        if fid == 0:
            continue

        el2[gid_coarse2==fid] = np.ones(len(el2[gid_coarse2==fid]))

    for fid, val in enumerate(l0):
        if val == 0:
            continue
        n = len(gid_coarse1[gid_coarse1==fid])

        l2[gid_coarse1==fid] = np.repeat(val, n)

    data_impress['verif_po'] = l2
    data_impress['verif_rest'] = el2
    data_impress.update_variables_to_mesh(['verif_po', 'verif_rest'])
    M.core.print(file='results/test_'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')

def dados_unitarios(data_impress):
    data_impress['hs'] = np.ones(len(data_impress['hs'])*3).reshape([len(data_impress['hs']), 3])
    data_impress['volume'] = np.ones(len(data_impress['volume']))
    data_impress['area'] = np.ones(len(data_impress['area']))
    data_impress['permeability'] = np.ones(data_impress['permeability'].shape)
    data_impress['k_harm'] = np.ones(len(data_impress['k_harm']))
    data_impress['dist_cent'] = np.ones(len(data_impress['dist_cent']))
    data_impress['transmissibility'] = np.ones(len(data_impress['transmissibility']))
    data_impress['pretransmissibility'] = data_impress['transmissibility'].copy()
    data_impress.export_to_npz()

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']
n = data_loaded['n_test']
load_operators = data_loaded['load_operators']

M, elements_lv0, data_impress, wells = initial_mesh(load=load, convert=convert)
# dados_unitarios(data_impress)

#######################
tpfa_solver = FineScaleTpfaPressureSolver(data_impress, elements_lv0, wells)
tpfa_solver.get_transmissibility_matrix_without_boundary_conditions()
# tpfa_solver.get_RHS_term()
# tpfa_solver.get_transmissibility_matrix()
multilevel_operators = MultilevelOperators(2, data_impress, M.multilevel_data, load=load_operators)

if load_operators:
    pass
else:
    multilevel_operators.run(tpfa_solver['Tini'])
