from implicit_impress.jacobian.symbolic_jacobian import symbolic_J as s_J
from packs.running.initial_mesh_properties import initial_mesh
from packs.directories import data_loaded
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.multiscale.operators.prolongation.AMS.Paralell.group_dual_volumes import group_dual_volumes_and_get_OP
from packs.biphasic.biphasic_ms.biphasic_multiscale import BiphasicTpfaMultiscale
from packs.pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from packs.adm.non_uniform.adm_method_non_nested import AdmNonNested

import matplotlib.pyplot as plt
from run_test_cases_mono import format_plot
from packs.adm import FIM_ADM
import scipy.sparse as sp
import numpy as np
import time

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



i=0
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
vec=np.zeros_like(M.faces.all)
vec[M.faces.internal]=1
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
        adm_conv, fs_iters=FIM_ADM.newton_iteration_finescale(F_Jacobian, Ts, adjs, data_impress_fs, time_step, wells)
        fs_conv, adm_iters=FIM_ADM.newton_iteration_ADM(vec, elements_lv0['faces_face_volumes'], elements_lv0['volumes_face_faces'], mlo, adm_method, F_Jacobian, Ts, adjs, data_impress, time_step, wells)
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
