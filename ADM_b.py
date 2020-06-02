from packs.running.initial_mesh_properties import initial_mesh
from packs.pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from packs.multiscale.preprocess.prep_neumann import NeumannSubdomains
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.directories import data_loaded
from packs.multiscale.operators.prolongation.AMS.Paralell.group_dual_volumes import group_dual_volumes_and_get_OP
from packs.multiscale.tpfalize_operator import tpfalize
import scipy.sparse as sp
from pymoab import types
import numpy as np
import time
import matplotlib.pyplot as plt
# from packs.adm.adm_method import AdmMethod
from packs.adm.non_uniform.adm_method_non_nested import AdmNonNested
# from packs.biphasic.biphasic_tpfa import BiphasicTpfa
from packs.biphasic.biphasic_ms.biphasic_multiscale import BiphasicTpfaMultiscale
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
def plot_operator(OP_AMS, primals):
    for i in range(len(primals)):
        tag_ams=M.core.mb.tag_get_handle("OP_AMS_"+str(primals[i]), 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        fb_ams = OP_AMS[:, primals[i]].toarray()
        M.core.mb.tag_set_data(tag_ams, M.core.all_volumes, fb_ams)

def write_file_with_tag_range(tag,lims):
    mc=M.core.mb.create_meshset()
    tag_ams=M.core.mb.tag_get_handle(tag, 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
    values=M.core.mb.tag_get_data(tag_ams,M.core.all_volumes,flat=True)
    entities=np.array(M.core.all_volumes)[(values>lims[0]) & (values<lims[1])]
    M.core.mb.add_entities(mc,entities)
    M.core.mb.write_file('results/'+tag+'.vtk', [mc])

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
    # M.core.print(file='results/test_'+ str(0), extension='.vtk', config_input='input_cards/print_settings0.yml')
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
    b1 = BiphasicTpfaMultiscale(M, data_impress, elements_lv0, wells)

multilevel_operators = MultilevelOperators(n_levels, data_impress, elements_lv0, M.multilevel_data, load=load_operators, get_correction_term=get_correction_term)
mlo = multilevel_operators

T, b = b1.get_T_and_b()
try:
    perms=np.load("flying/permeability.npy")
    perms_xx=perms[:,0]
    data_impress["perm_x"]=perms_xx
except:
    pass

if load_operators:
    pass
else:
    multilevel_operators.run_paralel(b1['Tini'], M.multilevel_data['dual_structure_level_1'], 0, False)


PP2=mlo['prolongation_level_'+str(1)]
mlo=multilevel_operators



tpfa_solver = FineScaleTpfaPressureSolver(data_impress, elements_lv0, wells)
tpfa_solver.run()
neta_lim=1.0
elements_lv0['neta_lim']=neta_lim
OP=group_dual_volumes_and_get_OP(mlo, T, M, data_impress, tpfa_solver, neta_lim=neta_lim)

# import pdb; pdb.set_trace()
# mlo=tpfalize(M,mlo,data_impress)
mlo['prolongation_lcd_level_1']=sp.find(mlo['prolongation_level_1'])
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
#
# meshset_volumes = M.core.mb.create_meshset()
# M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)

meshset_volumes = M.core.mb.create_meshset()
meshset_faces = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)
M.core.mb.add_entities(meshset_faces, M.core.all_faces)

OP_AMS=mlo['prolongation_level_1']
OR_AMS=mlo['restriction_level_1']



# plot_operator(OP_AMS,np.arange(OP_AMS.shape[1]))
# write_file_with_tag_range('OP_AMS_63',[0,np.inf])

Tc=OR_AMS*T*OP_AMS
bc=OR_AMS*b
from scipy.sparse import linalg
pc=linalg.spsolve(Tc,bc)

adm_method.organize_ops_adm(mlo, 1)
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
# M.core.mb.write_file('results/testt_00'+'.vtk', [meshset_volumes])
########## for plotting faces adm resolution in a faces file
int_faces=M.faces.internal
adjs=M.faces.bridge_adjacencies(int_faces,2,3)
ad0=adjs[:,0]
ad1=adjs[:,1]
# #######################################

neumann_subds=NeumannSubdomains(elements_lv0, adm_method.ml_data, data_impress, wells)
nn = 7000
pp = 100
cont = 1

verif = True
pare = False
np.save('results/jac_iterarion.npy',np.array([0]))
while verif:

    for level in range(1, n_levels):
        adm_method.organize_ops_adm(mlo, level)

    adm_method.set_level_wells_3()
    # gid1=data_impress["GID_1"]
    # data_impress["LEVEL"][ad0[gid1[ad0]!=gid1[ad1]]]=0
    # data_impress["LEVEL"][ad1[gid1[ad0]!=gid1[ad1]]]=0
    # data_impress["LEVEL"][ad1[ad0!=ad1]]=0
    # data_impress["LEVEL"][data_impress['saturation']>0.2]=0
    adm_method.set_saturation_level_simple()

    adm_method.solve_multiscale_pressure(T, b)

    adm_method.set_pms_flux(wells, neumann_subds)

    ####################################
    op_adm = adm_method._data['adm_prolongation_level_1']
    or_adm = adm_method._data['adm_restriction_level_1']
    Tc=or_adm*T*op_adm
    lcd=sp.find(Tc)
    lc=lcd[0]
    cc=lcd[1]
    dc=lcd[2]
    diags=lc==cc
    off_diag=lc!=cc
    val_diag=dc[diags]
    lcc_diag=lc[diags]
    val_off_diag=dc[off_diag]

    lc_off_diag=lc[off_diag]
    val_diag_off_diag=val_diag[lcc_diag][lc_off_diag]
    netas=val_off_diag/val_diag_off_diag
    netasp=netas>0
    netas=netas[netasp]
    I=lc[off_diag][netasp]
    J=cc[off_diag][netasp]
    if len(netas)>0:
        print(netas.max())
        if netas.max()>0:
            # import pdb; pdb.set_trace()
            afff=1

    b1.get_velocity_faces()

    b1.get_flux_volumes()

    b1.run_2()

    if cont % nn == 0:
        import pdb; pdb.set_trace()

    gid_0 = data_impress['GID_0'][data_impress['LEVEL']==0]
    gid_1 = data_impress['GID_0'][data_impress['LEVEL']==1]
    adm_method.set_adm_mesh_non_nested(v0=gid_0, v1=gid_1, pare=True)

    if cont % pp == 0:
        print("Creating_file")
        meshset_plot_faces=M.core.mb.create_meshset()
        lv=data_impress['LEVEL']
        gid_coarse=data_impress['GID_1']
        bounds_coarse=int_faces[gid_coarse[ad0]!=gid_coarse[ad1]]
        lvs0=int_faces[(lv[ad0]==0) | (lv[ad1]==0)]
        facs_plot=np.concatenate([bounds_coarse,lvs0])
        M.core.mb.add_entities(meshset_plot_faces,np.array(M.core.all_faces)[facs_plot])
        data_impress.update_variables_to_mesh()
        M.core.mb.write_file('results/testt_'+str(cont)+'.vtk', [meshset_volumes])
        M.core.mb.write_file('results/faces_'+str(cont)+'.vtk', [meshset_plot_faces])
        print("File created at time-step: ",cont)

    T, b = b1.get_T_and_b()
    print('timestep: ',cont)
    cont += 1
