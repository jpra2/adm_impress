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
    rr = set(np.where(l0>0)[0])
    rr2 = set(np.where(el0>0)[0])
    if rr & rr2 != rr2:
        import pdb; pdb.set_trace()

def mostrar_2(i, data_impress, M, op, rest, gid0, gid_coarse1, gid_coarse2):
    l0 = np.concatenate(op[:,i].toarray())
    el0 = np.concatenate(rest[i].toarray())
    el2 = np.zeros(len(gid0))
    l2 = el2.copy()
    cont = 0
    for fid, val in enumerate(el0):
        if val == 0:
            cont += 1
            continue
        else:
            el2[gid_coarse1==fid] = np.ones(len(el2[gid_coarse1==fid]))

    for fid, val in enumerate(l0):
        if val == 0:
            continue
        n = len(gid_coarse1[gid_coarse1==fid])

        l2[gid_coarse1==fid] = np.repeat(val, n)

    data_impress['verif_po'] = l2
    data_impress['verif_rest'] = el2
    data_impress.update_variables_to_mesh(['verif_po', 'verif_rest'])
    M.core.print(file='results/test_'+ str(0), extension='.vtk', config_input='input_cards/print_settings0.yml')
    import pdb; pdb.set_trace()

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
get_correction_term = data_loaded['get_correction_term']
n_levels = int(data_loaded['n_levels'])
_debug = data_loaded['_debug']

M, elements_lv0, data_impress, wells = initial_mesh()

######################
tpfa_solver = FineScaleTpfaPressureSolver(data_impress, elements_lv0, wells)
# tpfa_solver.get_transmissibility_matrix_without_boundary_conditions()
T, b = tpfa_solver.run()
# tpfa_solver.get_RHS_term()
# tpfa_solver.get_transmissibility_matrix()

multilevel_operators = MultilevelOperators(n_levels, data_impress, elements_lv0, M.multilevel_data, load=load_operators, get_correction_term=get_correction_term)

if load_operators:
    pass
else:
    if get_correction_term:
        multilevel_operators.run(tpfa_solver['Tini'], total_source_term=b, q_grav=data_impress['flux_grav_volumes'])
        # multilevel_operators.run_paralel_ant0(tpfa_solver['Tini'], total_source_term=b, q_grav=data_impress['flux_grav_volumes'])
    else:
        # multilevel_operators.run(tpfa_solver['Tini'])
        multilevel_operators.run_paralel(tpfa_solver['Tini'], M.multilevel_data['dual_structure_level_1'], 0, False)
        # multilevel_operators.run_paralel_ant0(tpfa_solver['Tini'])

# import pdb; pdb.set_trace()
