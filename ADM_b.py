from packs.running.initial_mesh_properties import initial_mesh
from packs.pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from packs.multiscale.preprocess.prep_neumann import NeumannSubdomains
from packs.adm.non_uniform import monotonize_adm
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.directories import data_loaded
from packs.multiscale.operators.prolongation.AMS.Paralell.group_dual_volumes import group_dual_volumes_and_get_OP
from packs.multiscale.operators.prolongation.AMS.ams_tpfa_new1 import AMSTpfa
from packs.multiscale.tpfalize_operator import tpfalize
from packs.adm.non_uniform import monotonic_adm_subds
# from packs.multiscale.create_local_submatrices_LHS import get_local_matrices_and_global_ids
from packs.multiscale.create_local_submatrices_LHS import LocalLU
from packs.multiscale.correction_function.CF import get_correction_function
import scipy.sparse as sp
from pymoab import types
import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
folder=np.load('flying/folder.npy')[0]

# from packs.adm.adm_method import AdmMethod
from packs.adm.non_uniform.adm_method_non_nested import AdmNonNested
# from packs.biphasic.biphasic_tpfa import BiphasicTpfa
from packs.biphasic.biphasic_ms.biphasic_multiscale import BiphasicTpfaMultiscale
from packs.tpfa.biphasic.load_or_preprocess_biphasic import preprocessar, carregar
from packs.tpfa.biphasic import TpfaBiphasicCons
from packs.tpfa.monophasic import TpfaMonophasic
from packs.multiscale.test_conservation import ConservationTest

def tpfalize_matrix(Mat, ids_coarse):
    lcd=sp.find(mat)
    import pdb; pdb.set_trace()

def save_multilevel_results():
    t_comp.append(t1-t0)
    vals_n1_adm.append(adm_method.n1_adm)
    vals_vpi.append(b1.vpi)
    vals_delta_t.append(b1.delta_t)
    vals_wor.append(b1.wor)

def plot_field(di, field):

    x=di['centroid_volumes'][:,1]
    y=di['centroid_volumes'][:,0]
    z=di['pressure']
    zd=di['DUAL_1']
    zc=di['coupled_flag']
    nx=len(np.unique(x))
    ny=len(np.unique(y))
    X=x.reshape(ny,nx)
    Y=y.reshape(ny,nx)
    Z=z.reshape(ny,nx)

    plt.close('all')
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)
    ax.zaxis.set_tick_params(labelsize=15)
    ax.view_init(30, 120)
    z[z<0.0]=np.nan
    z[z>1.0]=np.nan
    ax.plot_surface(X, Y, Z, cmap='jet',alpha=0.8, vmin=0.0,vmax=1.0)
    fig.colorbar(plt.cm.ScalarMappable(cmap='jet'), ax=ax,fraction=0.03)
    plt.savefig('results/'+folder+'/ms/'+ms_case+'/'+field+'.svg',transparent=True)



def export_multilevel_results(vals_n1_adm,vals_vpi,vals_delta_t,vals_wor, t_comp,
    el2, elinf, es_L2, es_Linf,vpis_for_save, ep_haji_L2,ep_haji_Linf, er_L2, er_Linf,
    ev_L2,ev_Linf):
    vals_n1_adm=100*np.array(vals_n1_adm+[len(data_impress['GID_0'])])/len(data_impress['pressure'])
    vals_vpi=np.array(vals_vpi)
    vals_delta_t=np.array(vals_delta_t)
    vals_wor=np.array(vals_wor)
    t_comp=np.array(t_comp)
    el2=np.array(el2)*100
    elinf=np.array(elinf)*100
    es_L2=np.array(es_L2)*100
    es_Linf=np.array(es_Linf)*100
    ev_L2=np.array(ev_L2)*100
    ev_Linf=np.array(ev_Linf)*100

    ep_haji_L2=np.array(ep_haji_L2)
    ep_haji_Linf=np.array(ep_haji_Linf)
    er_L2=np.array(er_L2)
    er_Linf=np.array(er_Linf)
    vars=[vals_vpi,vals_n1_adm,vals_delta_t,vals_wor, t_comp, el2, elinf, es_L2,
        es_Linf, vpis_for_save, ep_haji_L2,ep_haji_Linf, er_L2, er_Linf, ev_L2, ev_Linf]
    names=['vpi','n1_adm', 'delta_t', 'wor', 't_comp', 'el2', 'elinf', 'es_L2',
        'es_Linf', 'vpis_for_save','ep_haji_L2','ep_haji_Linf', 'er_L2', 'er_Linf', 'ev_L2', 'ev_Linf']
    for i in range(len(vars)):
        np.save('results/'+folder+'/ms/'+ms_case+names[i]+'.npy',vars[i])

def plot_operator(T,OP_AMS, primals, marker='',file_prefix=None):
    OP_AMS_c=T*OP_AMS
    for i in range(len(primals)):
        tag_ams=M.core.mb.tag_get_handle("OP_AMS_"+marker+str(primals[i]), 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        fb_ams = OP_AMS[:, primals[i]].toarray()
        tag_tp=M.core.mb.tag_get_handle("TP_"+marker+str(primals[i]), 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        tp_ams = OP_AMS_c[:, primals[i]].toarray()
        M.core.mb.tag_set_data(tag_ams, M.core.all_volumes, fb_ams)
        M.core.mb.tag_set_data(tag_tp, M.core.all_volumes, tp_ams)
        if file_prefix:
            ms=M.core.mb.create_meshset()
            M.core.mb.add_entities(ms,M.core.all_volumes)
            M.core.mb.write_file(file_prefix+str(i),[ms])

def plot_net_flux(OR_AMS,OP_AMS):
    Tc=OR_AMS*T*OP_AMS
    diag=np.array(Tc[range(Tc.shape[0]),range(Tc.shape[0])])[0]
    diag=np.array(T[range(T.shape[0]),range(T.shape[0])])[0]

    diag_orig=data_impress['initial_diag']
    diag_orig[diag_orig==0]=1
    diagf=np.array(T[range(T.shape[0]),range(T.shape[0])])[0]
    gids1=data_impress['GID_1']
    gids0=data_impress['GID_0']
    OP_AMS_c=OP_AMS.copy()

    phi_diag=OP_AMS_c[gids0,gids1].toarray()[0]
    OP_AMS_c[gids0,gids1]=0
    phi_off=OP_AMS_c.tocsc().max(axis=1).T.toarray()[0]


    data_impress['raz_diag'] = ((1-phi_diag)/phi_diag)*abs((diagf-diag_orig)/diag_orig)

    gid1=data_impress['GID_1']
    diags_f=diag[gid1]
    mm=T*OP_AMS
    mm2=mm.copy()
    # mm[mm>0]=0
    mm2[mm<0]=0

    # mms=np.array(mm.sum(axis=1)).T[0]
    mms=-(mm.max(axis=1)-mm.min(axis=1)).T.toarray()[0]
    mms2=np.array(mm.sum(axis=1)).T[0]
    data_impress['raz_flux_tag']=mms/diags_f
    # if (mms/diags_f).max() > 1:

    data_impress['raz_pos']=mms2/diags_f

def get_finescale_netas():
    t0=time.time()
    val_lim=1
    tc=OR_AMS*T*OP_AMS
    nc=Tc.shape[0]
    dc=np.array(tc[range(nc),range(nc)])[0]
    tp=T*OP_AMS
    # np.set_printoptions(1)
    gids0=data_impress['GID_0']
    gids1=data_impress['GID_1']
    print(time.time()-t0,"initie")
    tpn=tp.copy()
    # tpn[tpn>0]=0
    lines=gids0
    df=dc[gids1]
    t0=time.time()
    dd=sp.csc_matrix((1/df,(lines,lines)),shape=T.shape)
    nf=dd*tpn

    netasp=nf.max(axis=1).T.toarray()[0]

    lcd=mlo['prolongation_lcd_level_1']
    l=lcd[0]
    c=lcd[1]
    d=lcd[2]
    teste=True
    print(time.time()-t0,'dldldl')
    while teste:
        # netasp[wells['all_wells']]=0
        ls=np.arange(len(netasp))[netasp>val_lim]
        d2=d.copy()
        for li in ls:
            d2[l==li]=0
        OP_AMS_c=sp.csc_matrix((d2, (l,c)),shape=OP_AMS.shape)
        # OP_AMS_c=OP_AMS.copy()
        # OP_AMS_c[netasp>val_lim,:]=0
        t2p=(T*OP_AMS_c).tolil()
        # t2p[gids0,gids1]=0
        n2f=dd*t2p
        netas2p=n2f.max(axis=1).T.toarray()[0]
        netas2p[wells['all_wells']]=0
        if (netas2p>netasp).sum()==0:
            teste=False
        netasp[netas2p>netasp]=netas2p[netas2p>netasp]

    data_impress['nfp']=netasp
    ff=OP_AMS[gids0,gids1].toarray()[0]

    data_impress['nfn']=(1-ff)/ff

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
meshset_volumes = M.core.mb.create_meshset()
meshset_faces = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)
M.core.mb.add_entities(meshset_faces, M.core.all_faces)
wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate, phisical_properties = preprocessar(M, data_impress, wells)
# wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate, phisical_properties = carregar()
#########total_gravity_velocity,
monophasic = TpfaMonophasic()
biphasic = TpfaBiphasicCons()
biphasic_data['krw'], biphasic_data['kro'] = biphasic.get_krw_and_kro(biphasic_data['saturation'])
mob_w, mob_o = biphasic.get_mobilities_w_o(biphasic_data['krw'], biphasic_data['kro'])
biphasic_data['upwind_w'], biphasic_data['upwind_o'] = biphasic.set_initial_upwind_internal_faces(
    biphasic_data['saturation'],
    elements.get('volumes_adj_internal_faces'),
    wells['ws_inj'],
    elements.internal_faces,
    elements.volumes_to_faces(elements.volumes),
    elements.boundary_faces,
    elements.get('map_internal_faces')
)
biphasic_data['mob_w_internal_faces'] = mob_w[elements.get('volumes_adj_internal_faces')[biphasic_data['upwind_w']]]
biphasic_data['mob_o_internal_faces'] = mob_o[elements.get('volumes_adj_internal_faces')[biphasic_data['upwind_o']]]

biphasic_data['transmissibility_faces'] = biphasic.get_transmissibility_faces(
    geom['areas'],
    elements.internal_faces,
    elements.boundary_faces,
    elements.get('volumes_adj_internal_faces'),
    elements.get('volumes_adj_boundary_faces'),
    biphasic_data['mob_w_internal_faces'],
    biphasic_data['mob_o_internal_faces'],
    mob_w,
    mob_o,
    rock_data['keq_faces']
)
##############################
biphasic_data['transmissibility_faces'][:] = 1.0
delta_s = biphasic_data['saturation'][elements.get('volumes_adj_internal_faces')]
delta_s = np.absolute(delta_s[:,1] - delta_s[:,0])
delta_s = delta_s > 0
biphasic_data['transmissibility_faces'][elements.internal_faces[delta_s]] = 2.0

# ff = biphasic_data['transmissibility_faces'] > 1
# print(ff.sum())
#
# pdb.set_trace()
#############################


tt = biphasic_data['transmissibility_faces'].copy()

biphasic_data['g_velocity_w_internal_faces'], biphasic_data['g_velocity_o_internal_faces'] = biphasic.get_g_velocity_w_o_internal_faces(
    phisical_properties.gravity_vector,
    biphasic_data['mob_w_internal_faces'],
    biphasic_data['mob_o_internal_faces'],
    biphasic.properties.rho_w,
    biphasic.properties.rho_o,
    geom['hi'],
    rock_data['keq_faces'][elements.internal_faces]
)
g_total_velocity_internal_faces = biphasic_data['g_velocity_w_internal_faces'] + biphasic_data['g_velocity_o_internal_faces']

biphasic_data['g_source_w_internal_faces'], biphasic_data['g_source_o_internal_faces'] = biphasic.get_g_source_w_o_internal_faces(
    geom['areas'][elements.internal_faces],
    geom['u_direction_internal_faces'],
    biphasic_data['g_velocity_w_internal_faces'],
    biphasic_data['g_velocity_o_internal_faces']
)

wells2.add_gravity_2(
    elements.volumes,
    phisical_properties.gravity_vector,
    geom['centroid_volumes'],
    elements.get('volumes_adj_volumes_by_faces'),
    geom['centroid_nodes'],
    biphasic_data['saturation'],
    biphasic.properties.rho_w,
    biphasic.properties.rho_o
)

g_source_total_internal_faces = biphasic_data['g_source_w_internal_faces'] + biphasic_data['g_source_o_internal_faces']
g_source_total_volumes = phisical_properties.get_total_g_source_volumes(
    elements.volumes,
    elements.get('volumes_adj_internal_faces'),
    g_source_total_internal_faces
)

# g_source_total_volumes = np.zeros_like(g_source_total_volumes)
# data_impress['g_total_source_term'] = g_source_total_volumes
data_impress['g_total_source_term'] = g_source_total_volumes
T = monophasic.mount_transmissibility_matrix(
    biphasic_data['transmissibility_faces'][elements.internal_faces],
    elements.internal_faces,
    elements.get('volumes_adj_internal_faces'),
    elements.volumes
)

T_with_boundary, b = monophasic.get_linear_problem(
    wells2['ws_p'],
    wells2['ws_q'],
    wells2['values_p'],
    wells2['values_q'],
    g_source_total_volumes,
    T
)

data_impress['transmissibility'] = biphasic_data['transmissibility_faces'].copy()

# data_impress.update_variables_to_mesh()
# M.core.mb.write_file('results/testt_00'+'.vtk', [meshset_volumes])
# pdb.set_trace()

#########
# pdb.set_trace()

# separated_dual_structure = DualStructure(load=True)['dual_structure']

# local_matrices_and_global_ids=get_local_matrices_and_global_ids(separated_dual_structure, data_impress['transmissibility'])

# if biphasic:
#     b1 = BiphasicTpfaMultiscale(M, data_impress, elements_lv0, wells)

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
    # multilevel_operators.run_paralel(b1['Tini'], M.multilevel_data['dual_structure_level_1'], 0, False)
    multilevel_operators.run_paralel(T, M.multilevel_data['dual_structure_level_1'], 0, False)


# PP2=mlo['prolongation_level_'+str(1)]

mlo=multilevel_operators

# tpfa_solver = FineScaleTpfaPressureSolver(data_impress, elements_lv0, wells)
# tpfa_solver.run()
# neta_lim=np.load('flying/neta_lim_dual.npy')[0]
neta_lim=10
elements_lv0['neta_lim']=neta_lim

OP_AMS=mlo['prolongation_level_1']
phiks=np.array(OP_AMS[data_impress['GID_0'],data_impress['GID_1']])[0]
r_phiks=(1-phiks)/phiks

maxs=sp.csc_matrix((r_phiks,(data_impress['GID_0'],data_impress['GID_1'])),shape=OP_AMS.shape).max(axis=0).toarray()[0]

M.maxs=maxs

###########
dual_flags=data_impress['DUAL_1']
gids0=data_impress['GID_0']
gids1=data_impress['GID_1']
adjs=M.faces.bridge_adjacencies(M.faces.internal,2,3)
ads=adjs[(gids1[adjs][:,0]!=gids1[adjs][:,1]) & (dual_flags[adjs[:,0]]<=2) & (dual_flags[adjs[:,0]]<=2)]
vals=r_phiks[ads]
maxs_ds=vals.max(axis=1)


aa=data_impress['GID_1']
mexs=np.zeros(aa.max()+1)
np.maximum.at(mexs,aa,r_phiks)
##############
# import pdb; pdb.set_trace()
# data_impress['nfn'][:]=maxs[data_impress['GID_1']]
data_impress['nfn'][ads[:,0]]=maxs_ds
data_impress['nfn'][ads[:,1]]=maxs_ds
M.edges_and_vals=[maxs_ds,gids1[ads]]
OP=group_dual_volumes_and_get_OP(mlo, T, M, data_impress, tpfa_solver, neta_lim=neta_lim)
# OP=group_dual_volumes_and_get_OP(mlo, T, M, data_impress, tpfa_solver, neta_lim=20)

OP=group_dual_volumes_and_get_OP(mlo, T_with_boundary, M, data_impress, T, neta_lim=neta_lim)

##########################
ams_tpfa = AMSTpfa(data_impress['GID_0'], data_impress['GID_1'], data_impress['DUAL_1'])
Twire = ams_tpfa.get_Twire(T)
As = ams_tpfa.get_as_off_diagonal(Twire)
separated_dual_structures = DualStructure(load=True)['dual_structure']
local_lu_matrices = LocalLU()
# pdb.set_trace()
transmissibility = data_impress['transmissibility']
# transmissibility = np.ones(len(data_impress['transmissibility']))
local_lu_matrices.update_lu_objects(separated_dual_structures, transmissibility)

ml_data = M.multilevel_data
volumes_without_grav_level_0 = ml_data['volumes_without_grav_level_0']
data_impress['verif_rest'][:] = 0.0
data_impress['verif_rest'][volumes_without_grav_level_0] = 1.0


t0=time.time()
# b2 = g_source_total_volumes.copy()
b2 = b.copy()
b2[wells2['ws_p']] = 0.0
# b2[wells[]]
# b2[volumes_without_grav_level_0] = 0

# b2[data_impress['DUAL_1'] == 2] = 0.0

# b2=np.ones_like(b2)10000
# b2[wells['ws_q']] = -b[wells['ws_q']]
# b2[wells['ws_p']] = 0
# cfs = get_correction_function(local_lu_matrices.local_lu_and_global_ids, As, np.ones_like(b2))
cfs = get_correction_function(local_lu_matrices.local_lu_and_global_ids, As, b2)
# cfs = get_correction_function(local_lu_matrices.local_lu_and_global_ids, As, np.zeros_like(b2))
# cfs[wells['ws_p']] = 0
# cfs[:] = 0
data_impress['gama']=cfs
print('cf {}'.format(time.time()-t0))
# import pdb; pdb.set_trace()
#############################

# mlo=tpfalize(M,mlo,data_impress)
mlo['prolongation_lcd_level_1']=sp.find(mlo['prolongation_level_1'])
# adm_method = AdmMethod(wells['all_wells'], n_levels, M, data_impress, elements_lv0)
adm_method = AdmNonNested(wells['all_wells'], n_levels, M, data_impress, elements_lv0)

adm_method.restart_levels()
# adm_method.set_level_wells()
# adm_method.set_level_wells_2()
# adm_method.set_level_wells_3()
adm_method.set_level_wells_only()
adm_method.equalize_levels()

# adm_method.verificate_levels()
# adm_method.set_adm_mesh()
gids_0 = data_impress['GID_0']

adm_method.set_adm_mesh_non_nested(gids_0[data_impress['LEVEL']==0])
# adm_method.set_initial_mesh(mlo, T, b)
#
# meshset_volumes = M.core.mb.create_meshset()
# M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)

# meshset_volumes = M.core.mb.create_meshset()
# meshset_faces = M.core.mb.create_meshset()
# M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)
# M.core.mb.add_entities(meshset_faces, M.core.all_faces)

OP_AMS=mlo['prolongation_level_1']
OR_AMS=mlo['restriction_level_1']




# write_file_with_tag_range('OP_AMS_63',[0,np.inf])
from scipy.sparse import csc_matrix
import matplotlib as mpl
mpl.rc('font', **{'size'   : 22})

def save_matrices_as_png(matrices, names):

    for matrix,name in zip(matrices,names):
        plt.close('all')
        colors_pos = plt.cm.Blues(np.linspace(0, 1, 256))
        colors_neg = plt.cm.Reds(np.linspace(0, 1, 256))[::-1]
        if name=='OP_AMS' or name=='OP_ADM':
            all_colors=colors_pos
            v_center=1e-10
        else:
            all_colors = np.vstack((colors_neg, colors_pos))
            v_center=0

        mat_map = mpl.colors.LinearSegmentedColormap.from_list('mat_map', all_colors)
        divnorm = mpl.colors.DivergingNorm(vmin=matrix.toarray().min(), vcenter=v_center, vmax=matrix.toarray().max())

        plt.matshow(matrix.toarray(),cmap=mat_map,norm=divnorm,rasterized=True)
        plt.gca().set_yticks(range(matrix.shape[0]))
        plt.gca().set_xticks(range(matrix.shape[1]))
        plt.gcf().set_size_inches(matrix.shape[1],matrix.shape[0])
        plt.gca().set_xticks([x - 0.5 for x in plt.gca().get_xticks()][1:], minor='true')
        plt.gca().set_yticks([y - 0.5 for y in plt.gca().get_yticks()][1:], minor='true')
        plt.grid(which='minor')

        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                c=matrix[i,j]
                if abs(c)>1e-15:
                    c='{:.2f}'.format(c)
                    plt.text(j, i, str(c), va='center', ha='center',fontsize=15, color = 'lime',weight='bold')

        plt.savefig('results/'+name+'.png')

def save_matrices_as_png_with_highlighted_lines(matrices, names, Lines0):
    colors=['lime','cyan','y','m']

    cmap = mpl.colors.ListedColormap(colors)
    cmap.set_over('red')
    cmap.set_under('blue')

    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    mat_=np.array([[1,3],[0,2]])
    plt.matshow(mat_,cmap=cmap)
    plt.gca().set_yticks(range(mat_.shape[0]))
    plt.gca().set_xticks(range(mat_.shape[1]))
    plt.gcf().set_size_inches(mat_.shape[1]*2,mat_.shape[0]*2)
    ticks=[0.0,1.0]
    ticks_M=[-0.5,0.5,1.5]

    plt.gca().set_xticks(ticks, minor='true')
    plt.gca().set_yticks(ticks, minor='true')
    plt.gca().set_yticks(ticks_M, minor='false')
    plt.gca().set_xticks(ticks_M, minor='false')
    plt.gca().tick_params(which='major',width=3,color='black')
    plt.grid(which='minor',color='black',linewidth=2)
    plt.grid(which='major',color='black')
    plt.savefig('results/mesh.png')

    for matrix, name in zip(matrices,names):



        if matrix.shape[0]>4:
            Lines=[data_impress['GID_0'][data_impress['GID_1']==l] for l in Lines0]
        else:
            Lines=[[l] for l in Lines0]
        plt.close('all')
        plt.matshow(np.zeros_like(matrix.toarray()),cmap=plt.cm.Blues)
        for (lines, color) in zip(Lines,colors):
            for l in lines:
                plt.axhspan(l-0.5,l+0.5, color=color, alpha=0.7)
        plt.gca().set_yticks(range(matrix.shape[0]))
        plt.gca().set_xticks(range(matrix.shape[1]))
        plt.gcf().set_size_inches(matrix.shape[1],matrix.shape[0])
        plt.gca().set_xticks([x - 0.5 for x in plt.gca().get_xticks()][1:], minor='true')
        plt.gca().set_yticks([y - 0.5 for y in plt.gca().get_yticks()][1:], minor='true')
        plt.grid(which='minor')
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                c=matrix[i,j]
                if abs(c)>1e-15:
                    c='{:.2f}'.format(c)
                    plt.text(j, i, str(c), va='center', ha='center',fontsize=15)
        plt.savefig('results/'+name+'.png')

# Tc=OR_AMS*T*OP_AMg_source_total_volumes = np.zeros_like(g_source_total_volumes)
Tc=OR_AMS*T_with_boundary*OP_AMS
bc=OR_AMS*(b - T_with_boundary*cfs)
from scipy.sparse import linalg
pc=linalg.spsolve(Tc,bc)

adm_method.organize_ops_adm(mlo, 1)
pms=OP_AMS*pc + cfs
OP_ADM = adm_method['adm_prolongation_level_1']

OR_ADM = adm_method['adm_restriction_level_1']
# Tcadm=OR_ADM*T*OP_ADM
Tcadm=OR_ADM*T_with_boundary*OP_ADM
bcadm = OR_ADM*(b - T_with_boundary*cfs)
# bcadm = OR_ADM*(b - T*cfs)
pcadm=linalg.spsolve(Tcadm,bcadm)
padm=OP_ADM*pcadm + cfs
# padm=OP_ADM*pcadm
data_impress['pms_adm_pressure'] = padm
data_impress.update_variables_to_mesh()

# data_impress.update_variables_to_mesh()
# M.core.mb.write_file('results/trash.vtk',[meshset_volumes])
# import pdb; pdb.set_trace()

'''
matrices=[T,Tc,OP_AMS,T*OP_AMS, OP_ADM, Tcadm,T*OP_ADM]
names=['T','Tc','OP_AMS', 'TP', 'OP_ADM', 'Tcadm','TP_ADM']
save_matrices_as_png(matrices,names)

diags=np.array(Tc[range(Tc.shape[0]),range(Tc.shape[0])])[0]

gids1=data_impress['GID_1']
gids0=data_impress['GID_0']
dia_vec=csc_matrix(np.array([diags[gids1]]).T)
vals=1/diags[gids1]
mult=csc_matrix((vals,(gids0,gids0)),shape=(len(gids0),len(gids0)))
tp_mod=T*OP_AMS
tp_mod[gids0,gids1]=0
mat_neta=mult*tp_mod
mat_neta[mat_neta<0]=0
max_neta=mat_neta.max(axis=1).tocsc()

matricesh=[T,OP_AMS,T*OP_AMS,Tc, mat_neta,dia_vec, max_neta]
namesh=['T_hi','OP_AMS_hi','TP_hi','Tc_hi','neta','dia_vec','max_neta']
lines=[0,1,2,3]
save_matrices_as_png_with_highlighted_lines(matricesh,namesh, lines)

'''

# pf=linalg.spsolve(T,b)
pf=linalg.spsolve(T_with_boundary,b)
eadm=np.linalg.norm(abs(padm-pf))/np.linalg.norm(pf)
eams=np.linalg.norm(abs(pms-pf))/np.linalg.norm(pf)
print("erro_adm: {}, erro_ams: {}".format(eadm,eams))
data_impress['pressure'] = padm
# data_impress['tpfa_pressure'] = adm_method.solver.direct_solver(T, b)
data_impress['tpfa_pressure'] = pf.copy()
data_impress['erro_p'] = padm-pf

conservation_test = ConservationTest()

# total_flux_internal_faces, velocity_internal_faces = conservation_test.conservation_with_gravity(
#     elements.volumes,
#     data_impress['GID_1'],
#     g_source_total_internal_faces,
#     pf,
#     T,
#     ml_data['coarse_faces_level_1'],
#     ml_data['coarse_intersect_faces_level_1'],
#     geom['areas'],
#     ml_data['coarse_primal_id_level_1'],
#     elements.get('volumes_adj_internal_faces'),
#     ml_data['coarse_internal_faces_level_1'],
#     elements.get('map_internal_faces'),
#     geom['abs_u_normal_faces'],
#     biphasic_data['mob_w_internal_faces'],
#     biphasic_data['mob_o_internal_faces'],
#     phisical_properties.gravity_vector,
#     rock_data['keq_faces'],
#     biphasic.properties.rho_w,
#     biphasic.properties.rho_o,
#     geom['hi'],
#     g_source_total_volumes,
#     ml_data['vertex_level_1'],
#     g_source_total_internal_faces,
#     g_total_velocity
# )

t0 = time.time()
flux_coarse_volumes, test_vector, local_pressure = conservation_test.conservation_with_gravity(
    elements.volumes,
    data_impress['GID_1'],
    padm,
    T,
    ml_data['coarse_faces_level_1'],
    ml_data['coarse_intersect_faces_level_1'],
    geom['areas'],
    ml_data['coarse_primal_id_level_1'],
    elements.get('volumes_adj_internal_faces'),
    ml_data['coarse_internal_faces_level_1'],
    elements.get('map_internal_faces'),
    geom['abs_u_normal_faces'],
    phisical_properties.gravity_vector,
    rock_data['keq_faces'],
    biphasic.properties.rho_w,
    biphasic.properties.rho_o,
    geom['hi'],
    g_source_total_volumes,
    ml_data['vertex_level_1'],
    g_source_total_internal_faces,
    g_total_velocity_internal_faces,
    wells2['ws_p'],
    wells2['values_p'],
    wells2['ws_q'],
    wells2['values_q']
)
t1=time.time()
print(t1-t0)
data_impress['verif_po'][:] = 0.0
# data_impress['verif_rest'][:] = 0.0
# data_impress['verif_rest'][:] = local_pressure


for primal_id, flux in enumerate(flux_coarse_volumes):
    data_impress['verif_po'][elements.volumes[data_impress['GID_1'] == primal_id]] = flux
    # data_impress['verif_rest'][elements.volumes[data_impress['GID_1'] == primal_id]] = test_vector[primal_id]

data_impress.update_variables_to_mesh()
# M.core.mb.write_file('results/testt_00'+'.vtk', [meshset_volumes])
M.core.mb.write_file('results/trash'+'.vtk', [meshset_volumes])
pdb.set_trace()
########## for plotting faces adm resolution in a faces file
int_faces=M.faces.internal
adjs=M.faces.bridge_adjacencies(int_faces,2,3)
ad0=adjs[:,0]
ad1=adjs[:,1]
# #######################################

diagf=np.array(T[range(T.shape[0]),range(T.shape[0])])[0]

# data_impress['initial_diag']=diagf

neumann_subds=NeumannSubdomains(elements_lv0, adm_method.ml_data, data_impress, wells)
nn = 1000
pp = 1000
cont = 1
vpis_for_save=np.load('flying/vpis_for_save.npy')
vpis_for_vtk=np.load('flying/vpis_for_vtk.npy')
count_save=0
count_vtk=0

verif = True
pare = False
np.save('results/jac_iterarion.npy',np.array([0]))
phis_k=OP_AMS[data_impress['GID_0'],data_impress['GID_1']].toarray()[0]
data_impress['raz_pos']=(1-phis_k)/phis_k

phiK_raz_lim=np.load('flying/phiK_raz_lim.npy')[0]

preprocessed_primal_objects, critical_groups=monotonic_adm_subds.get_preprossed_monotonic_primal_objects(data_impress, elements_lv0, OP_AMS, neumann_subds.neumann_subds, phiK_raz_lim=phiK_raz_lim)
try:
    l_groups=np.concatenate([np.repeat(i,len(critical_groups[i])) for i in range(len(critical_groups))])
    groups_c=np.concatenate(critical_groups)
except:
    l_groups=np.array([0])
    groups_c=np.array([0])
    # l_groups=np.repeat(0,len(critical_groups[0]))
    # groups_c=critical_groups[0]
# else:
#     l_groups=np.array(critical_groups)
#     groups_c=critical_groups
ms_case=np.load("flying/ms_case.npy")[0]
data_impress['coupled_flag'][data_impress['DUAL_1']>=2]=0


coupl=100*(data_impress['coupled_flag']==1).sum()/len(data_impress['coupled_flag'])
np.save('results/'+folder+'/ms/'+ms_case+'/coupl'+'.npy',np.array([coupl]))

# adm_method.set_level_wells_3()
adm_method.set_level_wells_only()
vals_n1_adm=[]
vals_vpi=[]
vals_delta_t=[]
vals_wor=[]
t_comp=[]
el2=[]
elinf=[]
es_L2=[]
es_Linf=[]
ev_L2=[]
ev_Linf=[]

er_L2=[]
er_Linf=[]
ep_haji_L2=[]
ep_haji_Linf=[]


neta_lim_finescale=np.load('flying/neta_lim_finescale.npy')[0]
type_of_refinement=np.load('flying/type_of_refinement.npy')[0]
delta_sat_max=np.load('flying/delta_sat_max.npy')[0]

while verif:

    t00=time.time()
    transmissibility=data_impress['transmissibility']
    if type_of_refinement=='uni':
        volumes, netasp_array=monotonic_adm_subds.get_monotonizing_volumes(preprocessed_primal_objects, transmissibility)
        maxs=np.zeros(len(np.unique(volumes)))
        np.maximum.at(maxs,volumes,netasp_array)

        data_impress['nfp'][np.unique(volumes)]=maxs
        # netasp_array=np.maximum(netasp_array,netasp_array*data_impress['raz_phi'][volumes])
        # vols_orig=volumes[netasp_array>neta_lim_finescale]

        #try:
            # import pdb; pdb.set_trace()
            # netasp_array[data_impress['LEVEL'][volumes]==0]=neta_lim_finescale+1
        vols_orig=monotonic_adm_subds.get_monotonizing_level(l_groups, groups_c, critical_groups,data_impress,elements_lv0, volumes,netasp_array, neta_lim_finescale)
        # except:
        #     print('treta')
        #     import pdb; pdb.set_trace()
        #     vols_orig=data_impress['GID_0'][(maxs>neta_lim_finescale) | (data_impress['LEVEL']==0)]

            # vols_orig=data_impress['GID_0'][data_impress['LEVEL']==0]
        #
        # vols_orig=monotonic_adm_subds.get_monotonizing_level(l_groups, groups_c, critical_groups,data_impress,volumes,netasp_array, neta_lim_finescale)
        # vols_orig=np.array([])

    gid_0 = data_impress['GID_0'][data_impress['LEVEL']==0]
    gid_1 = data_impress['GID_0'][data_impress['LEVEL']==1]

    adm_method.set_adm_mesh_non_nested(v0=gid_0, v1=gid_1, pare=True)

    t0=time.time()
    for level in range(1, n_levels):
        adm_method.organize_ops_adm(mlo, level)
    '''# or_adm=adm_method._data['adm_restriction_level_1']
    # op_adm=adm_method._data['adm_prolongation_level_1']

    # idl1=data_impress['LEVEL_ID_1']
    # map1=np.zeros(idl1.max()+1)
    # map1[idl1]=data_impress['GID_1']
    # monotonize_adm.verify_monotonize_adm(or_adm, T, op_adm, neta_lim,map1)
    # np.concatenate(np.array(critical_groups)'''
    # adm_method.set_level_wells_3()
    adm_method.set_level_wells_only()

    if type_of_refinement=='uni':
        if len(vols_orig)>0:
            adm_method.set_monotonizing_level(vols_orig)
        adm_method.set_saturation_level_simple(delta_sat_max)
    else:
        adm_method.set_saturation_level_uniform(delta_sat_max)

    t0=time.time()
    # from implicit_impress.jacobian.symbolic_jacobian import symbolic_J as s_J
    from implicit_impress.jacobian.impress_assembly import assembly

    # M.pressure[:]=np.array([OP_ADM*linalg.spsolve(OR_ADM*T*OP_ADM,OR_ADM*b)]).T
    # JJ=assembly(M,0.00001)
    # J=JJ.J
    # nvols=int(J.shape[0]/2)
    # Jsp=J[nvols:,nvols:]
    # Jpp=J[0:nvols,0:nvols]
    # Jss=J[nvols:,nvols:]
    # Fs=JJ.q[nvols:]
    # S=(OR_ADM*Fs-OR_ADM*Jsp*OR_ADM.T*linalg.spsolve(OR_ADM*T*OP_ADM,OR_ADM*b))/(OR_ADM*OR_ADM.T)[np.arange(OR_ADM.shape[0]),np.arange(OR_ADM.shape[0])]
    # import pdb; pdb.set_trace()

    adm_method.solve_multiscale_pressure(T, b)

    adm_method.set_pms_flux(wells, neumann_subds) #

    b1.get_velocity_faces()

    b1.get_flux_volumes()

    b1.run_2()
    t1=time.time()

    pms=data_impress['pressure']
    if True:
        # np.save('flying/original_ms_solution.npy',pms)
        po=pms.copy()
        vo=data_impress['velocity_faces']
        # np.save('flying/velocity_faces_AMS.npy',vo)
    else:

        po=np.load('flying/original_ms_solution.npy')
        vo=np.load('flying/velocity_faces_AMS.npy')

    er_L2.append(np.linalg.norm(T*pms-b)/np.linalg.norm(T*po-b))

    er_Linf.append((T*pms-b).max()/(T*po-b).max())

    ep_haji_L2.append(np.linalg.norm(abs(pms-pf))/np.linalg.norm(pf-po))
    ep_haji_Linf.append(abs(pms-pf).max()/(pf-po).max())

    pf=linalg.spsolve(T,b)

    vf=np.load('flying/velocity_faces_finescale.npy')
    vadm=data_impress['velocity_faces']

    ev_L2.append(np.linalg.norm(abs(vf-vadm).max(axis=1))/np.linalg.norm(abs(vf).max(axis=1)))
    ev_Linf.append(abs(vf-vadm).max()/abs(vf).max())


    eadm_2=np.linalg.norm(abs(pms-pf))/np.linalg.norm(pf)
    eadm_inf=abs(pms-pf).max()/pf.max()
    el2.append(eadm_2)
    elinf.append(eadm_inf)

    p_wells=pms[wells['all_wells']]
    # pms=(pms-p_wells.min())/(p_wells.max()-p_wells.min())
    data_impress['pressure']=pms
    data_impress['DUAL_PRESSURE'][data_impress['DUAL_1']>=2]=pms[data_impress['DUAL_1']>=2]

    p_wells=pf[wells['all_wells']]
    # pf=(pf-p_wells.min())/(p_wells.max()-p_wells.min())
    data_impress['tpfa_pressure']=pf
    # plot_field(data_impress,'pressure')
    # plot_field(data_impress,'tpfa_pressure')
    save_multilevel_results()
    if vpis_for_save[count_save]<b1.vpi:
        sat_adm=data_impress['saturation']
        sat_f=np.load('flying/saturation_'+str(vpis_for_save[count_save])+'.npy')
        es_L2.append(np.linalg.norm(sat_f-sat_adm)/np.linalg.norm(sat_f))
        es_Linf.append(abs(sat_f-sat_adm).max()/sat_f.max())
        count_save+=1

    if len(vpis_for_vtk)>0 and vpis_for_vtk[count_vtk]<b1.vpi:
        refinement=100*((data_impress['LEVEL']==0).sum()-len(p_wells))/len(data_impress['LEVEL'])
        np.save('results/'+folder+'/ms/'+ms_case+'/refinement'+'.npy',np.array([refinement]))

        print(b1.vpi,'vpi')
        print("Creating_file")
        meshset_plot_faces=M.core.mb.create_meshset()
        lv=data_impress['LEVEL']
        gid_coarse=data_impress['GID_1']
        bounds_coarse=int_faces[gid_coarse[ad0]!=gid_coarse[ad1]]
        lvs0=int_faces[(lv[ad0]==0) | (lv[ad1]==0)]
        facs_plot=np.concatenate([bounds_coarse,lvs0])
        M.core.mb.add_entities(meshset_plot_faces,np.array(M.core.all_faces)[facs_plot])

        data_impress['raz_pos'][:]=0
        data_impress['raz_pos'][data_impress['DUAL_1']==2]=data_impress['raz_phi'][data_impress['DUAL_1']==2]
        data_impress.update_variables_to_mesh()
        file_count=str(int(100*vpis_for_vtk[count_vtk]))

        M.core.mb.write_file('results/'+folder+'/ms/'+ms_case+'vtks/volumes_'+file_count+'.vtk', [meshset_volumes])
        M.core.mb.write_file('results/'+folder+'/ms/'+ms_case+'vtks/faces_'+file_count+'.vtk', [meshset_plot_faces])

        if vpis_for_vtk[count_vtk]==vpis_for_vtk.max():
            export_multilevel_results(vals_n1_adm,vals_vpi,vals_delta_t,vals_wor,
            t_comp, el2, elinf, es_L2, es_Linf,vpis_for_vtk[:count_vtk+1],
            ep_haji_L2,ep_haji_Linf, er_L2, er_Linf,ev_L2,ev_Linf)
            verif=False
        print("File created at time-step: ",cont)
        count_vtk+=1


    # if cont % nn == 0 or count_save==len(vpis_for_save)-1:
    #     export_multilevel_results(vals_n1_adm,vals_vpi,vals_delta_t,vals_wor, t_comp, el2, elinf, es_L2, es_Linf,vpis_for_save[:count_save])
    #     verif=False


    T, b = b1.get_T_and_b()
    print(b1.wor, b1.vpi, adm_method.n1_adm, time.time()-t00, eadm_2)
    cont += 1
