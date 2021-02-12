from packs.running.initial_mesh_properties import initial_mesh
from packs.directories import data_loaded
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.multiscale.operators.prolongation.AMS.Paralell.group_dual_volumes import group_dual_volumes_and_get_OP
from packs.biphasic.biphasic_ms.biphasic_multiscale import BiphasicTpfaMultiscale
from packs.pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from packs.adm.non_uniform.adm_method_non_nested import AdmNonNested
from implicit_impress.jacobian.impress_assembly import assembly
import scipy.sparse as sp
import numpy as np

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
OP_AMS=mlo['prolongation_level_1']
OR_AMS=mlo['restriction_level_1']
adm_method.organize_ops_adm(mlo, 1)
OP_ADM = adm_method['adm_prolongation_level_1']
OR_ADM = adm_method['adm_restriction_level_1']

def newton_iteration(M, time_step, rel_tol=1e-8):
    converged=False
    count=0
    while not converged:
        FIM=assembly(M,time_step)
        J=FIM.J
        q=FIM.q
        sol=-ADM_solver(J, q, R, P)
        n=int(len(q)/2)
        M.pressure[:]+=np.array([sol[0:n]]).T
        M.swns[:]+=np.array([sol[n:]]).T
        M.swns[:][M.swns[:]<0.0]=0.0
        converged=(abs(np.linalg.norm(sol)/np.linalg.norm(np.concatenate([M.pressure[:],M.swns[:]]).T[0])))<rel_tol
        count+=1
    print(count,(abs(np.linalg.norm(sol)/np.linalg.norm(np.concatenate([M.pressure[:],M.swns[:]]).T[0]))))
def ADM_solver(J, q, R, P):
    sol=P*sp.linalg.spsolve(R*J*P,R*q)
    return sol

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

time_step=0.00001
for i in range(50):
    M.swn1s[:]=M.swns[:]
    M.swn1s[-1]=1.0
    newton_iteration(M, time_step)
    data_impress['pressure']=M.pressure[:]
    data_impress['swns']=M.swns[:]
    meshset_volumes=M.core.mb.create_meshset()
    M.core.mb.add_entities(meshset_volumes,np.array(M.core.all_volumes))
    data_impress.update_variables_to_mesh()
    M.core.mb.write_file('results/biphasic/FIM_'+str(i)+'.vtk', [meshset_volumes])
import pdb; pdb.set_trace()
