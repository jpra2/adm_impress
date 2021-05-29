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

# @profile
def newton_iteration_ADM(data_impress, time_step, wells, rel_tol=1e-3):
    converged=False
    count=0
    dt=time_step
    data_impress['swn1s']=data_impress['swns'].copy()
    adm_method.data_impress['saturation']=data_impress['swns'].copy()
    all_ids=data_impress['GID_0']
    not_prod=np.setdiff1d(all_ids,wells['all_wells'])
    while not converged:
        data_impress['swns'][wells['ws_inj']]=1
        adm_method.restart_levels()
        # adm_method.set_level_wells_only()
        adm_method.set_level_wells_3()
        adm_method.set_saturation_level_homogeneo(delta_sat_max)
        adm_method.set_adm_mesh_non_nested(gids_0[adm_method.data_impress['LEVEL']==0])
        adm_method.organize_ops_adm(mlo, 1)
        OP_ADM = adm_method._data['adm_prolongation_level_1']
        OR_ADM = adm_method._data['adm_restriction_level_1']
        R, P = get_R_and_P(OR_ADM, OP_ADM)

        FIM=assembly(adjs, Ts, data_impress, dt, wells, F_Jacobian)
        J=FIM.J
        q=FIM.q
        sol=ADM_solver(J, q, R, P)
        if count==0 and data_impress['swns'].sum()==1:
            sol[data_impress['LEVEL_ID_1'][wells['ws_p']]]=-wells['values_p'].copy()
        sol=P*sol
        if count==0 and data_impress['swns'].sum()==1:
            sol[wells['ws_p']]=0

        n=int(len(q)/2)

        data_impress['pressure']-=sol[0:n]
        data_impress['swns']-=sol[n:]
        data_impress['swns'][wells['ws_inj']]=1
        adm_method.data_impress['saturation']=data_impress['swns'].copy()

        converged=max(abs(sol[n:][not_prod]))<rel_tol
        print(max(abs(sol[n:][not_prod])),max(abs(sol)),'ADM')

        count+=1
        if count>20:
            print('excedded maximum number of iterations ADM')
            return False, count
            break

    data_impress['swns'][wells['ws_prod']]=data_impress['swns'][wells['viz_prod']].sum()/len(wells['viz_prod'])
    sats=data_impress['swns'].copy()
    sats=get_sat_averager(sats,data_impress)
    # sats=OR_ADM.T*(OR_ADM*sats/np.array(OR_ADM.sum(axis=1)).T[0])
    data_impress['swns']=sats.copy()
    data_impress['swns'][sats>1]=1.0
    data_impress['swns'][sats<0]=0.0
    return True, count

def newton_iteration_finescale(data_impress, time_step, wells, rel_tol=1e-3):
    converged=False
    count=0
    dt=time_step
    data_impress['swn1s']=data_impress['swns'].copy()
    all_ids=data_impress['GID_0']
    not_prod=np.setdiff1d(all_ids,wells['all_wells'])
    while not converged:
        data_impress['swns'][wells['ws_inj']]=1

        FIM=assembly(adjs, Ts, data_impress, dt, wells, F_Jacobian)

        J=FIM.J
        q=FIM.q
        sol=-sp.linalg.spsolve(J, q)
        n=int(len(q)/2)
        data_impress['pressure']+=sol[0:n]

        data_impress['swns']+=sol[n:]
        data_impress['swns'][wells['ws_inj']]=1
        converged=max(abs(sol[n:][not_prod]))<rel_tol
        print(max(abs(sol[n:][not_prod])),max(abs(sol)),'fs')
        count+=1
        if count>20:
            print('excedded maximum number of iterations finescale')
            return False, count


    data_impress['swns'][wells['ws_prod']]=data_impress['swns'][wells['viz_prod']].sum()/len(wells['viz_prod'])
    return True, count

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
# @profile
def print_results(data_impress,name, i):
    print('saving file at timestep: {}'.format(i))
    M.pressure[:]=data_impress['pressure'].copy()
    M.swns[:]=data_impress['swns'].copy()
    M.swn1s[:]=data_impress['swn1s'].copy()
    M.LEVEL[:]=data_impress['LEVEL'].copy()
    meshset_volumes=M.core.mb.create_meshset()
    M.core.mb.add_entities(meshset_volumes,np.array(M.core.all_volumes))
    # data_impress.update_variables_to_mesh()
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
    # data_impress.update_variables_to_mesh()
    M.core.mb.write_file('results/biphasic_FIM/'+name+'f_'+str(i)+'.vtk', [meshset_plot_faces])
    print('File saved at time-step', i)

def plot_matrix(Mat, name):
    plt.close('all')
    if Mat.shape[0]>100:
        print('are you sure to plot matrix with shape: {}'.format(Mat.shape))
        import pdb; pdb.set_trace()
    Tc=Mat.toarray()
    Tc[Tc==0]=np.nan
    plt.matshow(Tc)
    data=sp.find(Mat)
    for i, j, z in zip(data[0],data[1],data[2]):
        plt.text(j, i, '{:0.1f}'.format(z), ha='center', va='center', color='white', size=3)
    plt.xticks(np.unique(data[1]),size=5)
    plt.yticks(np.unique(data[0]),size=5)
    plt.savefig('results/biphasic_FIM/'+name+'.svg')

def ADM_solver(J, q, R, P):
    # plot_matrix(R*J*P,'coarse')
    # plot_matrix(J,'fine')
    # plot_matrix(P,'prolongation')
    # plot_matrix(R,'restriction')
    # import pdb; pdb.set_trace()
    sol=sp.linalg.spsolve(R*J*P,R*q)
    return sol
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
def restart_simulation(data_impress):
    data_impress['swns']=np.zeros_like(data_impress['swns'])
    data_impress['pressure']=np.zeros_like(data_impress['swns'])
    data_impress['swns'][wells['ws_inj']]=1.0
    data_impress['pressure'][wells['ws_p']]=wells['values_p']

data_impress_fs=data_impress._data.copy()
data_impress_ADM=data_impress._data.copy()
restart_simulation(data_impress_fs)
restart_simulation(data_impress)
vpi=[]
pva=[]
p=[]
s=[]
i=1
# time_step=0.08
time_step=0.0005
# time_step=1.0
count_save=0
wells['viz_prod']=np.concatenate(elements_lv0['volumes_face_volumes'][wells['ws_prod']])
delta_sat_max=0.03
max_vpi=0.01
vpis_for_save=np.arange(0,0.91,0.01)

ep_l2=[]
ep_linf=[]
es_l2=[]
es_linf=[]
vpi=[]
pva_ADM=[]
vpi_n=0
while (i==1) or (max(vpi)<max_vpi):
    wells['count']=i-1
    if i==5:
        time_step*=2
    if i==10:
        time_step*=2.5
    if i==20:
        time_step*=2

    converged=False
    while not converged:
        p_adm_ant=data_impress['pressure'].copy()
        s_adm_ant=data_impress['swns'].copy()
        p_fs_ant=data_impress_fs['pressure'].copy()
        s_fs_ant=data_impress_fs['swns'].copy()
        adm_conv, fs_iters=newton_iteration_finescale(data_impress_fs, time_step, wells)
        fs_conv, adm_iters=newton_iteration_ADM(data_impress, time_step, wells)
        converged= adm_conv & fs_conv
        if not converged:
            print('reducing time-step from {} to {}'.format(time_step,time_step/1.5))
            time_step/=1.5
            data_impress['pressure']=p_adm_ant.copy()
            data_impress['swns']=s_adm_ant.copy()
            data_impress['pressure']=p_fs_ant.copy()
            data_impress['swns']=s_fs_ant.copy()
        elif max(fs_iters,adm_iters)<4:
            print('increasing time-step from {} to {}'.format(time_step,1.5*time_step))
            time_step*=1.5

    # vpi_n-=wells['values_q'].sum()*time_step/(len(data_impress['pressure'])*0.3)
    vpi_n=(data_impress['swns'].sum()+1)/(len(data_impress['pressure']))
    padm=data_impress['pressure']
    pf=data_impress_fs['pressure']
    sadm=data_impress['swns']
    sf=data_impress_fs['swns']

    vpi.append(vpi_n)
    ep_l2.append(np.linalg.norm(padm-pf)/np.linalg.norm(pf))
    ep_linf.append(abs(padm-pf).max()/abs(pf).max())
    es_l2.append(np.linalg.norm(sadm-sf)/np.linalg.norm(sf))
    es_linf.append(abs(sadm-sf).max()/abs(sf).max())
    pva_ADM.append((data_impress['LEVEL_ID_1'].max()+1)/len(data_impress_ADM['swns']))
    # print('simulation: {}, time-step: {}, vpi: {}'.format(simulation, i, vpi_n))
    print(vpi_n,i,'volume poroso injetado e iteração')
    if (count_save<len(vpis_for_save)) and (vpi_n>vpis_for_save[count_save]):
        print_results(data_impress,'FIM_'+'ADM'+'_',count_save)
        print_results(data_impress_fs,'FIM_'+'finescale'+'_',count_save)
        count_save+=1
    i+=1

ab=np.array(vpi)


plot_vars=[[        ab,        ab,        ab,         ab,         ab],     # Abcissas
           [     ep_l2,   ep_linf,     es_l2,    es_linf,    pva_ADM],     # Ordenadas
           [ 'lin_lin', 'lin_lin', 'lin_lin',  'lin_lin',  'lin_lin'],
           [   'ep_l2', 'ep_linf',   'es_l2',  'es_linf',      'pva']]     # Escalas dos Eixos
for i in range(len(plot_vars[0])):
    try:
        plt.close('all')
        plt.plot(plot_vars[0][i],plot_vars[1][i])
        format_plot(plot_vars[2][i],plot_vars[0][i],plot_vars[1][i])
        plt.savefig('results/biphasic_FIM/'+plot_vars[3][i]+'.svg', bbox_inches='tight')
    except:
        import pdb; pdb.set_trace()
import pdb; pdb.set_trace()
