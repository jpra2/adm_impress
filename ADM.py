from packs.running.initial_mesh_properties import initial_mesh
from packs.pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.directories import data_loaded
import scipy.sparse as sp
import numpy as np
import time
# from packs.adm.adm_method import AdmMethod
from packs.adm.non_uniform.adm_method_non_nested import AdmNonNested
'''
def get_gids_and_primal_id(gids, primal_ids):
    gids2 = np.unique(gids)
    primal_ids2 = []
    for i in gids2:
        primal_id = np.unique(primal_ids[gids==i])
        if len(primal_id) > 1:
            raise ValueError('erro get_gids_and_primal_id')
        primal_ids2.append(primal_id[0])
    primal_ids2 = np.array(primal_ids2)
    return gids2, primal_ids2'''

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
    # M.core.print(file='test_'+ str(0), extension='.vtk', config_input='input_cards/print_settings0.yml')
    import pdb; pdb.set_trace()

# def dados_unitarios(data_impress):
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
#
if load_operators:
    pass
else:
    # multilevel_operators.run(tpfa_solver['Tini'])
    multilevel_operators.run_paralel(tpfa_solver['Tini'])

mlo=multilevel_operators

# adm_method = AdmMethod(wells['all_wells'], n_levels, M, data_impress, elements_lv0)
adm_method = AdmNonNested(wells['all_wells'], n_levels, M, data_impress, elements_lv0)
T, b = tpfa_solver.run()

adm_method.restart_levels()
# adm_method.set_level_wells()
# adm_method.verificate_levels()
# adm_method.set_adm_mesh()
gids_0 = data_impress['GID_0']

# adm_method.set_adm_mesh_non_nested(gids_0[data_impress['LEVEL']==0])
# adm_method.set_initial_mesh(mlo, T, b)

adm_method.organize_ops_adm(mlo['prolongation_level_1'],
                            mlo['restriction_level_1'],
                            1)
# adm_method.plot_operator(adm_method[adm_method.adm_op_n+'1'], mlo['prolongation_level_1'], 0)

if n_levels > 2:
    adm_method.organize_ops_adm(mlo['prolongation_level_2'],
                                mlo['restriction_level_2'],
                                2)


adm_method.set_initial_mesh(mlo, T, b)
############################teste#################################
OP_AMS=mlo['prolongation_level_1']
OR_AMS=mlo['restriction_level_1']
Tc=OR_AMS*T*OP_AMS
bc=OR_AMS*b
from scipy.sparse import linalg
pc=linalg.spsolve(Tc,bc)
pf=linalg.spsolve(T,b)
pms=OP_AMS*pc
OP_ADM = adm_method['adm_prolongation_level_1']
OR_ADM = adm_method['adm_restriction_level_1']
Tcadm=OR_ADM*T*OP_ADM
bcadm = OR_ADM*b
pcadm=linalg.spsolve(Tcadm,bcadm)
padm=OP_ADM*pcadm

eadm=np.linalg.norm(abs(padm-pf))/np.linalg.norm(pf)
eams=np.linalg.norm(abs(pms-pf))/np.linalg.norm(pf)
print("erro_adm: {}, erro_ams: {}".format(eadm,eams))
# import pdb; pdb.set_trace()
# adm_method.organize_ops_adm_level_1( OP_AMS, OR_AMS, level, _pcorr=None)
#########################################################################
Tc2=Tc.copy()
Tc2.setdiag(0)
DTc=np.array(Tc[range(Tc.shape[0]),range(Tc.shape[0])])[0]
MTc=Tc2.min(axis=1).toarray().T[0]
netasams=abs(MTc/DTc)

Tcadm2=Tcadm.copy()
Tcadm2.setdiag(0)
DTcadm=np.array(Tcadm[range(Tcadm.shape[0]),range(Tcadm.shape[0])])[0]
MTcadm=Tcadm2.min(axis=1).toarray().T[0]
netasadm=abs(MTcadm/DTcadm)
# np.set_printoptions(precision=0)
# print(Tc.toarray())
# M.multilevel_data._data['adm_prolongation_level_1']
print("netamax: adm: {}, ams: {}".format(netasadm.max(), netasams.max()))
import pdb; pdb.set_trace()
adm_method.solve_multiscale_pressure(T, b)
adm_method.set_pcorr()
data_impress['pcorr'][data_impress['LEVEL']==0] = data_impress['pms'][data_impress['LEVEL']==0]

data_impress['pressure'] = adm_method.solver.direct_solver(T, b)
data_impress['erro'] = np.absolute((data_impress['pressure'] - data_impress['pms'])/data_impress['pms'])
data_impress['erro_pcorr_pdm'] = np.absolute(data_impress['pcorr'] - data_impress['pms'])

data_impress.update_variables_to_mesh()
# M.core.print(folder='results' file='test_'+ str(0), extension='.vtk', config_input='input_cards/print_settings0.yml')
