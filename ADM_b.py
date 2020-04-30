from packs.running.initial_mesh_properties import initial_mesh
from packs.pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.directories import data_loaded
from packs.multiscale.operators.prolongation.AMS.Paralell.group_dual_volumes import group_dual_volumes_and_get_OP
import scipy.sparse as sp
import numpy as np
import time
import matplotlib.pyplot as plt
# from packs.adm.adm_method import AdmMethod
from packs.adm.non_uniform.adm_method_non_nested import AdmNonNested
from packs.biphasic.biphasic_tpfa import BiphasicTpfa
plt.close('all')
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
biphasic = data_loaded['biphasic']

M, elements_lv0, data_impress, wells = initial_mesh()

if biphasic:
    b1 = BiphasicTpfa(M, data_impress, elements_lv0, wells)

multilevel_operators = MultilevelOperators(n_levels, data_impress, elements_lv0, M.multilevel_data, load=load_operators, get_correction_term=get_correction_term)
mlo = multilevel_operators
T, b = b1.get_T_and_b()


if load_operators:
    pass
else:
    # multilevel_operators.run(tpfa_solver['Tini'])
    multilevel_operators.run_paralel(b1['Tini'], M.multilevel_data['dual_structure_level_1'], 0, False)
PP2=mlo['prolongation_level_'+str(1)]
mlo=multilevel_operators
tpfa_solver = FineScaleTpfaPressureSolver(data_impress, elements_lv0, wells)
tpfa_solver.run()

OP=group_dual_volumes_and_get_OP(mlo, T, M, data_impress, tpfa_solver, neta_lim=1.0)
# adm_method = AdmMethod(wells['all_wells'], n_levels, M, data_impress, elements_lv0)
adm_method = AdmNonNested(wells['all_wells'], n_levels, M, data_impress, elements_lv0)

adm_method.restart_levels()
# adm_method.set_level_wells()
# adm_method.set_level_wells_2()
adm_method.set_level_wells_3()
adm_method.equalize_levels()

# adm_method.verificate_levels()
# adm_method.set_adm_mesh()
gids_0 = data_impress['GID_0']

adm_method.set_adm_mesh_non_nested(gids_0[data_impress['LEVEL']==0])
# adm_method.set_initial_mesh(mlo, T, b)

meshset_volumes = M.core.mb.create_meshset()
meshset_faces = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)
M.core.mb.add_entities(meshset_faces, M.core.all_faces)

OP_AMS=mlo['prolongation_level_1']
OR_AMS=mlo['restriction_level_1']
Tc=OR_AMS*T*OP_AMS
bc=OR_AMS*b
from scipy.sparse import linalg
pc=linalg.spsolve(Tc,bc)
adm_method.organize_ops_adm(mlo['prolongation_level_'+str(1)],
                            mlo['restriction_level_'+str(1)],
                            1)
pms=OP_AMS*pc
OP_ADM = adm_method['adm_prolongation_level_1']

OR_ADM = adm_method['adm_restriction_level_1']
Tcadm=OR_ADM*T*OP_ADM
bcadm = OR_ADM*b
pcadm=linalg.spsolve(Tcadm,bcadm)
padm=OP_ADM*pcadm


pf=linalg.spsolve(T,b)
eadm=np.linalg.norm(abs(padm-pf))/np.linalg.norm(pf)
eams=np.linalg.norm(abs(pms-pf))/np.linalg.norm(pf)
print("erro_adm: {}, erro_ams: {}".format(eadm,eams))
data_impress['pressure'] = padm
data_impress['tpfa_pressure'] = adm_method.solver.direct_solver(T, b)

data_impress.update_variables_to_mesh()
M.core.mb.write_file('results/testt_00'+'.vtk', [meshset_volumes])

nn = 50
cont = 1

verif = True
pare = False
np.save('results/jac_iterarion.npy',np.array([0]))
while verif:
    # multilevel_operators = MultilevelOperators(n_levels, data_impress, elements_lv0, M.multilevel_data, load=load_operators, get_correction_term=get_correction_term)
    # multilevel_operators.run_paralel(b1['Tini'], M.multilevel_data['dual_structure_level_1'], 0, False)
    # mlo = multilevel_operators
    for level in range(1, n_levels):

        adm_method.organize_ops_adm(mlo['prolongation_level_'+str(level)],
                                    mlo['restriction_level_'+str(level)],
                                    level)

    # op1 = adm_method['adm_prolongation_level_1']
    # or1 = adm_method['adm_restriction_level_1']
    # import pdb; pdb.set_trace()
    adm_method.restart_levels()
    # adm_method.set_level_wells_2()
    adm_method.set_level_wells_3()
    adm_method.set_saturation_level()
    adm_method.equalize_levels()

    adm_method.solve_multiscale_pressure(T, b)
    adm_method.set_pms_flux(b1['Tini'], wells, pare=pare)
    b1.get_velocity_faces()
    b1.get_flux_volumes()
    # b1.print_test()
    # p2 = adm_method.solver.direct_solver(T, b)
    # data_impress['pressure'] = p2
    # data_impress['erro'] = np.absolute((p2-data_impress['pms']))
    b1.run_2()

    # M.core.mb.write_file('results/testt_f'+str(cont)+'.vtk', [meshset_faces])

    if cont % nn == 0:
        import pdb; pdb.set_trace()

    # adm_method.restart_levels_2()
    # adm_method.set_level_wells()
    #### adm_method.restart_levels()
    #### # adm_method.set_level_wells_2()
    #### adm_method.set_level_wells_3()
    #### adm_method.set_saturation_level()
    # adm_method.set_saturation_level_imposed_joined_coarse()
    # adm_method.set_saturation_level_imposed_bound_level_continuity()
    # adm_method.set_saturation_level_original()
    adm_method.equalize_levels()
    gid_0 = data_impress['GID_0'][data_impress['LEVEL']==0]
    gid_1 = data_impress['GID_0'][data_impress['LEVEL']==1]

    adm_method.set_adm_mesh_non_nested(v0=gid_0, v1=gid_1, pare=True)
    # b1.print_test()

    # n=0
    data_impress.update_variables_to_mesh()
    M.core.print(folder='results', file='test_'+ str(n), extension='.vtk', config_input='input_cards/print_settings0.yml')

    T, b = b1.get_T_and_b()

    cont += 1
