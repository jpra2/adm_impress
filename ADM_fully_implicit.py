from implicit_impress.jacobian.symbolic_jacobian import symbolic_J as s_J
from packs.running.initial_mesh_properties import initial_mesh
from packs.directories import data_loaded
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.multiscale.operators.prolongation.AMS.Paralell.group_dual_volumes import group_dual_volumes_and_get_OP
from packs.biphasic.biphasic_ms.biphasic_multiscale import BiphasicTpfaMultiscale
from packs.pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from packs.adm.non_uniform.adm_method_non_nested import AdmNonNested
from implicit_impress.jacobian.impress_assembly import assembly
import matplotlib.pyplot as plt
from run_test_cases_mono import format_plot
import scipy.sparse as sp
import numpy as np
import time

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']
n = data_loaded['n_test']
load_operators = data_loaded['load_operators']
get_correction_term = data_loaded['get_correction_term']
n_levels = int(data_loaded['n_levels'])
_debug = data_loaded['_debug']
biphasic = data_loaded['biphasic']

neta_lim=1.0

M, elements_lv0, data_impress, wells = initial_mesh()
b1 = BiphasicTpfaMultiscale(M, data_impress, elements_lv0, wells)
multilevel_operators = MultilevelOperators(n_levels, data_impress,
 elements_lv0, M.multilevel_data, load=load_operators,
 get_correction_term=get_correction_term)

T, b = b1.get_T_and_b()

multilevel_operators.run_paralel(b1['Tini'], M.multilevel_data['dual_structure_level_1'], 0, False)
mlo=multilevel_operators
tpfa_solver = FineScaleTpfaPressureSolver(data_impress, elements_lv0, wells)
tpfa_solver.run()
elements_lv0['neta_lim']=neta_lim
mlo['prolongation_lcd_level_1']=sp.find(mlo['prolongation_level_1'])
adm_method = AdmNonNested(wells['all_wells'], n_levels, M, data_impress, elements_lv0)
adm_method.restart_levels()
# adm_method.set_level_wells_3()
adm_method.set_level_wells()
# adm_method.equalize_levels()
gids_0 = data_impress['GID_0']
adm_method.set_adm_mesh_non_nested(gids_0[data_impress['LEVEL']==0])
# OP_AMS=mlo['prolongation_level_1']
# OR_AMS=mlo['restriction_level_1']
adm_method.organize_ops_adm(mlo, 1)
# OP_ADM = adm_method['adm_prolongation_level_1']
# OR_ADM = adm_method['adm_restriction_level_1']
#
adjs=M.faces.bridge_adjacencies(M.faces.internal_elements[:],2,3)
Ts=M.k_harm[M.faces.internal].T[0]
F_Jacobian = s_J()
vec=np.zeros_like(M.faces.all)
vec[M.faces.internal]=1

def newton_iteration_ADM(data_impress, time_step, wells, rel_tol=1e-10):
    print(data_impress['swns'].sum(),'admsum')
    converged=False
    count=0
    dt=time_step
    data_impress['swn1s']=data_impress['swns'].copy()
    p=data_impress['pressure'].copy()
    s=data_impress['swns'].copy()
    while not converged:
        adm_method.restart_levels()
        adm_method.set_level_wells_only()
        adm_method.set_saturation_level_simple(delta_sat_max)
        adm_method.set_adm_mesh_non_nested(gids_0[adm_method.data_impress['LEVEL']==0])
        adm_method.organize_ops_adm(mlo, 1)
        OP_ADM = adm_method._data['adm_prolongation_level_1']
        OR_ADM = adm_method._data['adm_restriction_level_1']
        R, P = get_R_and_P(OR_ADM, OP_ADM)
        FIM=assembly(adjs, Ts, data_impress, dt, wells, F_Jacobian)
        J=FIM.J
        q=FIM.q
        sol=-P*ADM_solver(J, q, R, P)
        n=int(len(q)/2)

        p+=sol[0:n]
        s+=sol[n:]
        print(s.sum(),'sums')
        data_impress['saturation']=data_impress['swns'].copy()

        converged=max(abs(sol[n:]))<rel_tol
        count+=1
        if count>10:
            pass
        print(count,max(abs(sol[n:])))
    sats=s.copy()
    sats=get_sat_averager(sats,data_impress)
    # sats=OR_ADM.T*(OR_ADM*sats/np.array(OR_ADM.sum(axis=1)).T[0])
    data_impress['swns']=sats.copy()
    data_impress['saturation']=sats.copy()
    return data_impress['pressure'], data_impress['swns']

def get_sat_averager(sats, data_impress):
    gids0=data_impress['GID_0']
    gids1=data_impress['GID_1']
    level=data_impress['LEVEL']
    gids_adm_c=data_impress['LEVEL_ID_1']
    for i in np.unique(data_impress['LEVEL_ID_1'][level==1]):
        vols=gids0[gids_adm_c==i]
        t00=time.time()
        faces=np.unique(np.concatenate(elements_lv0['volumes_face_faces'][vols]))
        faces=faces[vec[faces]==1]
        t0=time.time()
        adjs=np.vstack(elements_lv0['faces_face_volumes'][faces])
        adjs_int=adjs[(gids_adm_c[adjs[:,0]]==i) & (gids_adm_c[adjs[:,1]]==i)]
        if len(adjs_int)>0:
            mapv=np.zeros(vols.max()+1)
            mapv[vols]=np.arange(len(vols))
            la=mapv[adjs_int]
            lines=np.concatenate([la[:,0], la[:,1]])
            cols=np.concatenate([la[:,1], la[:,0]])
            nvols=len(vols)
            gr=sp.csc_matrix((np.ones_like(lines), (lines, cols)),shape=(nvols, nvols))
            n_l,labels=sp.csgraph.connected_components(gr,connection='strong')
            for i in range(n_l):
                gv=vols[labels==i]
                sats[gv] = sats[gv].sum()/len(gv)
    return sats

def newton_iteration_finescale(data_impress, time_step, wells, rel_tol=1e-5):
    print(data_impress['swns'].sum(),'fssum')
    converged=False
    count=0
    dt=time_step
    data_impress['swn1s']=data_impress['swns'].copy()
    while not converged:
        FIM=assembly(adjs, Ts, data_impress, dt, wells, F_Jacobian)
        J=FIM.J
        q=FIM.q
        sol=-sp.linalg.spsolve(J, q)
        n=int(len(q)/2)
        data_impress['pressure']+=sol[0:n]
        data_impress['swns']+=sol[n:]
        converged=max(abs(sol[n:]))<rel_tol
        count+=1
        print(count,max(abs(sol[n:])))
    return data_impress['pressure'], data_impress['swns']

def newton_iteration_fs(M, time_step, rel_tol=1e-5):
    converged=False
    count=0

    while not converged:
        FIM=assembly(M, time_step)
        J=FIM.J
        q=FIM.q
        sol=-sp.linalg.spsolve(J, q)

        n=int(len(q)/2)
        # M.pressure[:]+=np.array([sol[0:n]]).T
        # M.swns[:]+=np.array([sol[n:]]).T
        # M.swns[:][M.swns[:]>1.0]=1.0
        # M.swns[:][M.swns[:]<0.0]=0.0

        data_impress['pressure']+=sol[0:n]
        data_impress['swns'][data_impress['swns']>1]=1.0
        data_impress['swns'][data_impress['swns']<0]=0.0
        # data_impress['saturation']=M.swns[:].T[0]

        converged=max(abs(sol[n:]))<rel_tol
        count+=1
        print(count,max(abs(sol[n:])))

    print(count,max(abs(sol[n:])))

def get_R_and_P(OR_ADM, OP_ADM):
    lp, cp, dp = sp.find(OP_ADM)
    lr, cr, dr = sp.find(OR_ADM)
    n_f, n_ADM=OP_ADM.shape
    lP=np.concatenate([lp, cr+n_f])
    cP=np.concatenate([cp, lr+n_ADM])
    dP=np.concatenate([dp, dr])

    lR=np.concatenate([lr, lr+n_ADM])
    cR=np.concatenate([cr, cr+n_f])
    dR=np.concatenate([dr, dr])

    R=sp.csc_matrix((dR, (lR, cR)), shape=(2*n_ADM, 2*n_f))
    P=sp.csc_matrix((dP, (lP, cP)), shape=(2*n_f, 2*n_ADM))
    return R, P
@profile
def print_results(data_impress,name):
    M.pressure[:]=data_impress['pressure']
    M.swns[:]=data_impress['swns']
    M.swn1s[:]=data_impress['swn1s']
    meshset_volumes=M.core.mb.create_meshset()
    M.core.mb.add_entities(meshset_volumes,np.array(M.core.all_volumes))
    data_impress.update_variables_to_mesh()
    M.core.mb.write_file('results/biphasic_FIM/'+name+str(i)+'.vtk', [meshset_volumes])
    int_faces=M.faces.internal
    # adjs=M.faces.bridge_adjacencies(int_faces,2,3)
    ad0=adjs[:,0]
    ad1=adjs[:,1]
    meshset_plot_faces=M.core.mb.create_meshset()
    lv=data_impress['LEVEL']
    gid_coarse=data_impress['GID_1']
    bounds_coarse=int_faces[gid_coarse[ad0]!=gid_coarse[ad1]]
    lvs0=int_faces[(lv[ad0]==0) | (lv[ad1]==0)]
    facs_plot=np.concatenate([bounds_coarse,lvs0])
    M.core.mb.add_entities(meshset_plot_faces,np.array(M.core.all_faces)[facs_plot])
    data_impress.update_variables_to_mesh()
    M.core.mb.write_file('results/biphasic_FIM/'+name+'f_'+str(i)+'.vtk', [meshset_plot_faces])
    print('File saved at time-step', i)

def ADM_solver(J, q, R, P):
    sol=sp.linalg.spsolve(R*J*P,R*q)
    return sol

delta_sat_max=0.1
# time_step=0.001*(0.3*len(data_impress['pressure']))
time_step=0.0001


continue_old_simulation=False
if continue_old_simulation:
    try:
        vect = np.load('flying/saturations.npy')
        M.swns[:] = vect[:-1]
        i = int(vect[-1])
        M.pressure[:] = np.load('flying/pressures.npy')
        print('Saturation loaded from previous simulation!')
    except:
        i=0
        M.swns[:] = 0
        M.pressure[:] = 0
        print('Tried to continue but started a new simulation!')

else:
    i=0
    print('Started a new simulation as ordened!')

i0=i
data_impress['swn1s'][wells['ws_inj']]=1.0
data_impress['swns'][wells['ws_inj']]=1.0
data_impress['pressure'][wells['ws_p']]=wells['values_p']
data_impress_ADM=data_impress._data.copy()
data_impress_fs=data_impress_ADM.copy()
del data_impress
ep_l2=[0]
ep_linf=[0]
es_l2=[0]
es_linf=[0]
vpi=[0]
pva=[0]
max_vpi=0.001

while max(vpi)<max_vpi:
    print('Time-step {}, vpi/vpi_max {}/{}'.format(i,max(vpi),max_vpi) )
    if (i==2) or ((i==i0) and (i>2)):
        time_step*=5

    # np.save('flying/saturations.npy', np.append(M.swn1s[:],i))
    # np.save('flying/pressures.npy', M.pressure[:].T[0])
    print(data_impress_fs['swns'].sum(),'soksk')
    data_impress_ADM['pressure'], data_impress_ADM['swns'] = newton_iteration_ADM(data_impress_ADM, time_step, wells)
    print(data_impress_fs['swns'].sum(),'soksk2')
    # SADM=data_impress_ADM['swns']
    # Sfs=data_impress_fs['swns']
    # print(SADM.sum(),Sfs.sum(),'sum_sats')
    data_impress_fs['pressure'], data_impress_fs['swns'] = newton_iteration_finescale(data_impress_fs, time_step, wells)

    # PADM=data_impress_ADM['pressure']
    # Pfs=data_impress_fs['pressure']
    SADM=data_impress_ADM['swns']
    Sfs=data_impress_fs['swns']
    print(SADM.sum(),Sfs.sum(),'sum_sats')
    # vpi.append(data_impress_ADM['swns'].sum()*0.3/len(data_impress_ADM['swns']))
    # ep_l2.append(np.linalg.norm(PADM-Pfs)/np.linalg.norm(Pfs))
    # ep_linf.append(abs(PADM-Pfs).max()/abs(Pfs).max())
    # es_l2.append(np.linalg.norm(SADM-Sfs)/np.linalg.norm(Sfs))
    # es_linf.append(abs(SADM-Sfs).max()/abs(Sfs).max())
    # pva.append((data_impress_ADM['LEVEL_ID_1'].max()+1)/len(PADM))

    # if i//500 ==i/500:
    #     print_results(data_impress, 'FIM_ADM_')
    #     data_impress['pressure']=data_impress2['pressure'].copy()
    #     data_impress['swns']=data_impress2['swns'].copy()
    #     data_impress['swn1s']=data_impress2['swn1s'].copy()
    #     print_results(data_impress, 'FIM_fs_')
    # i+=1

# ab=vpi
#
# plot_vars=[[        ab,        ab,        ab,         ab,         ab],     # Abcissas
#            [     ep_l2,   ep_linf,     es_l2,    es_linf,        pva],     # Ordenadas
#            [ 'lin_lin', 'lin_lin', 'lin_lin',  'lin_lin',  'lin_lin'],
#            [   'ep_l2', 'ep_linf',   'es_l2',  'es_linf',      'pva']]     # Escalas dos Eixos
# for i in range(len(plot_vars[0])):
#     try:
#         plt.close('all')
#         plt.plot(plot_vars[0][i],plot_vars[1][i])
#         format_plot(plot_vars[2][i],plot_vars[0][i],plot_vars[1][i])
#         plt.savefig('results/biphasic_FIM/'+plot_vars[3][i]+'.svg', bbox_inches='tight')
#     except:
        # import pdb; pdb.set_trace()
import pdb; pdb.set_trace()
