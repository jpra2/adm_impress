from ..data_class.data_manager import DataManager
import numpy as np
import scipy.sparse as sp
from ..solvers.solvers_scipy.solver_sp import SolverSp
from ..flux_calculation.flux_tpfa import TpfaFlux2
import multiprocessing as mp
from .local_solution import LocalSolution
from. obj_infos import InfosForProcess
from ..directories import data_loaded
from ..errors.err import ConservativeVolumeError, PmsFluxFacesError
from scipy.sparse import linalg, find
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix
from ..directories import file_adm_mesh_def
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve_triangular
from packs.adm.smoothers.test_smoothers import compare_smoothers
# from packs.multiscale.correction_function.CF import CorrectionFunction
import time
import pdb


def get_levelantids_levelids(level_ids_ant, level_ids):

    gids2 = np.unique(level_ids_ant)
    level_ids2 = []
    for i in gids2:
        test = level_ids_ant==i
        level_id = np.unique(level_ids[test])
        if len(level_id) > 1:
            raise ValueError('erro get_level_id')
        level_ids2.append(level_id[0])

    level_ids2 = np.array(level_ids2)

    return gids2, level_ids2

def solve_local_local_problem(solver, neigh_intern_faces, transmissibility, volumes,indices_p, values_p, indices_q=[], values_q=[]):

    # t0 = transmissibility
    v02 = neigh_intern_faces
    indices_p2 = indices_p
    t0 = transmissibility
    n = len(volumes)
    # local_ids = np.arange(n)
    # map_volumes = dict(zip(volumes, local_ids))
    # v02 = np.zeros(v0.shape)
    # v02[:,0] = np.array([map_volumes[k] for k in v0[:, 0]])
    # v02[:,1] = np.array([map_volumes[k] for k in v0[:, 1]])
    # indices_p2 = np.array(map_volumes[k] for k in indices_p)
    # indices_q2 = np.array(map_volumes[k] for k in indices_p)

    lines = np.concatenate([v02[:, 0], v02[:, 1], v02[:, 0], v02[:, 1]])
    cols = np.concatenate([v02[:, 1], v02[:, 0], v02[:, 0], v02[:, 1]])
    data = np.concatenate([t0, t0, -t0, -t0])
    T = sp.csc_matrix((data, (lines, cols)), shape=(n, n))
    T = T.tolil()

    T[indices_p2] = 0
    T[indices_p2, indices_p2] = np.ones(len(indices_p2))

    b = np.zeros(n)
    b[indices_p] = values_p
    b[indices_q] += values_q
    x = solver(T.tocsc(), b)
    return x

def set_boundary_conditions(T: 'transmissibility matrix',
                            b: 'source term',
                            indices_q: 'indices flux prescription',
                            values_q: 'flux prescription',
                            indices_p: 'indices pressure prescription',
                            values_p: 'values pressure prescription'):
    n = T.shape[0]

    T = T.tolil()
    np = len(indices_p)
    nq = len(indices_q)

    T[indices_p] = sp.lil_matrix((np, n))
    T[indices_p, indices_p] = np.ones(np)

    b[indices_q] += values_q

    return T.tocsc(), b

def Jacobi(xini, T, b):
    b = b.reshape([len(b), 1])

    nf = len(xini)
    jacobi_options = file_adm_mesh_def['jacobi_options']
    n = jacobi_options['n_verif']
    _dev = file_adm_mesh_def['_dev']

    titer = time.time()

    ran=range(nf)
    D=T.diagonal()
    l_inv=range(nf)
    data_inv=1/D
    D_inv=sp.csc_matrix((data_inv,(l_inv,l_inv)),shape=(nf, nf))
    D=sp.csc_matrix((D,(l_inv,l_inv)),shape=(nf, nf))
    R=T-D
    x0=sp.csc_matrix(xini).transpose()
    cont=0
    for i in range(n):x0=D_inv*(b-R*x0)
    delta_ant=abs((D_inv*(b-R*x0)-x0)).max()
    cont+=n
    for i in range(n):x0=D_inv*(b-R*x0)
    delta=abs((D_inv*(b-R*x0)-x0)).max()
    cont+=n
    while  delta<0.6*delta_ant:
        delta_ant=delta
        for i in range(n):x0=D_inv*(b-R*x0)
        delta=abs((D_inv*(b-R*x0)-x0)).max()
        cont+=n
    x0=np.array(x0).T[0]

    if _dev:
        print(time.time()-titer,n, "iterou ")
    return x0

class AdmMethod(DataManager, TpfaFlux2):

    def __init__(self, all_wells_ids, n_levels, M, data_impress, elements_lv0, load=False):
        data_name = 'AdmMethod.npz'
        super().__init__(data_name=data_name, load=load)
        self.mesh = M
        self.elements_lv0 = elements_lv0
        self.ml_data = M.multilevel_data
        self.all_wells_ids = all_wells_ids
        self.n_levels = n_levels
        # self.delta_sat_max = 0.05
        # self.delta_sat_max = 2.0
        # self.delta_sat_max = np.load('flying/delta_sat_max.npy')[0]
        self.delta_sat_max = 0.1
        self.data_impress = data_impress
        self.number_vols_in_levels = np.zeros(self.n_levels+1, dtype=int)
        gids_0 = self.data_impress['GID_0']
        self.data_impress['LEVEL_ID_0'] = gids_0.copy()
        self.solver = SolverSp()

        from ..directories import data_loaded
        self.get_correction_term = data_loaded['get_correction_term']

        self.adm_op_n = 'adm_prolongation_level_'
        self.adm_rest_n = 'adm_restriction_level_'
        self.pcorr_n = 'pcorr_level_'

        self.n_cpu = mp.cpu_count()
        self.n_workers = self.n_cpu
        self._so_nv1 = False
        if self.n_levels == 2:
            self.so_nv1 = True

    def set_level_wells(self):
        self.data_impress['LEVEL'][self.all_wells_ids] = np.zeros(len(self.all_wells_ids))

        so_nv1 = self.so_nv1

        if so_nv1:
            self.data_impress['LEVEL'] = np.ones(len(self.data_impress['GID_0']), dtype=int)
            self.data_impress['LEVEL'][self.all_wells_ids] = np.zeros(len(self.all_wells_ids), dtype=int)

    def set_level_wells_2(self):
        self.data_impress['LEVEL'][self.all_wells_ids] = np.zeros(len(self.all_wells_ids))

    def set_adm_mesh(self):

        levels = self.data_impress['LEVEL']
        gids_0 = self.data_impress['GID_0']
        gids_1 = self.data_impress['GID_1']
        gids_2 = self.data_impress['GID_2']
        vvv2 = range(len(np.unique(gids_2)))
        n0 = len(levels)

        list_L1_ID = np.repeat(-1, n0)
        list_L2_ID = np.repeat(-1, n0)

        n1=0
        n2=0
        n_vols = 0
        # meshset_by_L2 = mb.get_child_meshsets(self.L2_meshset)
        print('\n')
        print("INICIOU GERACAO DA MALHA ADM")
        print('\n')

        for v2 in vvv2:
            #1
            # n_vols_l3 = 0
            nivel3 = True
            nivel2 = False
            nivel1 = False
            vols2 = gids_0[gids_2==v2]
            # gids_1_1 = gids_1[gids_2==v2]
            gids_1_1 = gids_1[vols2]
            vvv1 = np.unique(gids_1_1)
            n_vols_2 = len(vols2)
            conj_vols_1 = set()

            for v1 in vvv1:
                #2
                # elem_by_L1 = mb.get_entities_by_handle(m1)
                vols1 = vols2[gids_1_1==v1]
                nn1 = len(vols1)
                # n_vols += nn1
                # n_vols_l3 += nn1
                levels_vols_1 = levels[vols1]
                set_verif = set(levels_vols_1)

                if set([0]) & set_verif: # se houver volumes no nivel 0
                    #3
                    # volumes.append(elem_by_L1)
                    # meshsets_nv1.add(m1)
                    conj_vols_1.add(v1)
                    nivel3 = False
                    nivel1 = True
                    level = 0
                    list_L1_ID[vols1] = np.arange(n1, n1+nn1)
                    # list_L1_ID.append(np.arange(n1, n1+nn1))
                    list_L2_ID[vols1] = np.arange(n2, n2+nn1)
                    # list_L2_ID.append(np.arange(n2, n2+nn1))
                    levels[vols1] = np.repeat(level, nn1)
                    # list_L3_ID.append(np.repeat(level, nn1))
                    n1 += nn1
                    n2 += nn1
                #2
                elif set([1]) & set_verif: # se houver volumes no nivel 1
                    #3
                    # volumes.append(elem_by_L1)
                    # meshsets_nv2.add(m1)
                    conj_vols_1.add(v1)
                    nivel3 = False
                    nivel2 = True
                    level = 1
                    list_L1_ID[vols1] = np.repeat(n1, nn1)
                    list_L2_ID[vols1] = np.repeat(n2, nn1)
                    levels[vols1] = np.repeat(level, nn1)

                    n1 += 1
                    n2 += 1
            #1
            if nivel3:
                #2
                level = 2
                for v1 in vvv1:
                    #3
                    vols1 = vols2[gids_1_1==v1]
                    nn1 = len(vols1)
                    # volumes.append(elem_by_L1)
                    list_L1_ID[vols1] = np.repeat(n1, nn1)
                    # list_L1_ID.append(np.repeat(n1, nn1))
                    n1 += 1
                #2
                list_L2_ID[vols2] = np.repeat(n2, n_vols_2)
                # list_L2_ID.append(np.repeat(n2, n_vols_l3))
                levels[vols2] = np.repeat(level, n_vols_2)
                # list_L3_ID.append(np.repeat(level, n_vols_l3))
                n2 += 1
            #1
            else:
                #2
                vols_1_fora = set(vvv1) - conj_vols_1
                if vols_1_fora:
                    #3
                    for v1 in vols_1_fora:
                        #4
                        vols1 = vols2[gids_1_1==v1]
                        nn1 = len(vols1)
                        level = 1
                        list_L1_ID[vols1] = np.repeat(n1, nn1)
                        # list_L1_ID.append(np.repeat(n1, nn1))
                        list_L2_ID[vols1] = np.repeat(n2, nn1)
                        # list_L2_ID.append(np.repeat(n2, nn1))
                        levels[vols1] = np.repeat(level, nn1)
                        # list_L3_ID.append(np.repeat(level, nn1))
                        n1 += 1
                        n2 += 1

        self.data_impress['LEVEL_ID_1'] = list_L1_ID
        self.data_impress['LEVEL_ID_2'] = list_L2_ID
        self.data_impress['LEVEL'] = levels

        for i in range(self.n_levels):
            self.number_vols_in_levels[i] = len(levels[levels==i])

        self.n1_adm = n1
        self.n2_adm = n2

    def restart_levels(self):
        self.data_impress['LEVEL'] = np.repeat(-1, len(self.data_impress['LEVEL']))

    def restart_levels_2(self):
        self.data_impress['LEVEL'] = self.data_impress['INITIAL_LEVEL'].copy()

    def organize_ops_adm(self, OP_AMS, OR_AMS, level, _pcorr=None):

        gid_0 = self.data_impress['GID_0']
        gid_level = self.data_impress['GID_' + str(level)]
        gid_ant = self.data_impress['GID_' + str(level-1)]
        level_id = self.data_impress['LEVEL_ID_' + str(level)]
        level_id_ant = self.data_impress['LEVEL_ID_' + str(level-1)]
        levels = self.data_impress['LEVEL']
        dual_id = self.data_impress['DUAL_'+str(level)]
        OP_AMS = OP_AMS.tolil()

        if level == 1:
            OP_ADM, OR_ADM, pcorr = self.organize_ops_adm_level_1(OP_AMS, OR_AMS, level, _pcorr=_pcorr)

            self._data[self.adm_op_n + str(level)] = OP_ADM
            self._data[self.adm_rest_n + str(level)] = OR_ADM
            self._data[self.pcorr_n+str(level-1)] = pcorr
            return 0

        n_adm = len(np.unique(level_id))
        n_adm_ant = len(np.unique(level_id_ant))
        n1_adm = n_adm_ant
        n2_adm = n_adm

        if n_adm == n_adm_ant:
            OP_ADM = sp.identity(n_adm)
            self._data[self.adm_op_n + str(level)] = OP_ADM
            self._data[self.adm_rest_n + str(level)] = OP_ADM
            return 0

        lines=[]
        cols=[]
        data=[]

        lines_or=[]
        cols_or=[]
        data_or=[]

        vertices_lv = gid_0[self.data_impress['DUAL_1']==3]
        for i in range(2, self.n_levels+1):
            vertices_lv = np.intersect1d(vertices_lv, gid_0[self.data_impress['DUAL_'+str(i)]==3])

        ids_classic_level = self.data_impress['GID_'+str(level)][vertices_lv]
        ids_adm_level = self.data_impress['LEVEL_ID_'+str(level)][vertices_lv]
        AMS_TO_ADM = np.arange(len(ids_classic_level))
        AMS_TO_ADM[ids_classic_level] = ids_adm_level

        My_IDs_2 = set()
        for gid in gid_0:
            ID_ADM_ANT = level_id_ant[gid]
            if set([ID_ADM_ANT]) & set(My_IDs_2):
                continue
            My_IDs_2.add(ID_ADM_ANT)
            ID_ADM = level_id[gid]
            nivel = levels[gid]
            ID_AMS = gid_ant[gid]

            if nivel<level:
                lines.append([ID_ADM_ANT])
                cols.append([ID_ADM])
                data.append([1])
                lines_or.append(ID_ADM)
                cols_or.append(ID_ADM_ANT)
                data_or.append(1)
                #OP_ADM_2[ID_global][ID_ADM]=1
            else:
                lines_or.append(ID_ADM)
                cols_or.append(ID_ADM_ANT)
                data_or.append(1)
                ff = sp.find(OP_AMS[ID_AMS])
                ids_ADM = AMS_TO_ADM[ff[1]]
                lines.append(np.repeat(ID_ADM_ANT, len(ids_ADM)).astype(int))
                cols.append(ids_ADM)
                data.append(ff[2])

        lines = np.concatenate(lines)
        cols = np.concatenate(cols)
        data = np.concatenate(data)

        OP_ADM = sp.csc_matrix((data,(lines,cols)),shape=(n1_adm,n2_adm))
        OR_ADM = sp.csc_matrix((data_or,(lines_or,cols_or)),shape=(n2_adm,n1_adm))
        import pdb; pdb.set_trace()
        self._data[self.adm_op_n + str(level)] = OP_ADM
        self._data[self.adm_rest_n + str(level)] = OR_ADM

        return 0

    def organize_ops_adm_level_1(self, OP_AMS, OR_AMS, level, _pcorr=None):

        gid_0 = self.data_impress['GID_0']
        gid_level = self.data_impress['GID_' + str(level)]
        gid_ant = self.data_impress['GID_' + str(level-1)]
        level_id = self.data_impress['LEVEL_ID_' + str(level)]
        level_id_ant = self.data_impress['LEVEL_ID_' + str(level-1)]
        levels = self.data_impress['LEVEL']
        vertices = gid_0[self.data_impress['DUAL_1']==3]
        OP_AMS = OP_AMS.copy().tolil()

        AMS_TO_ADM = np.arange(len(gid_level[vertices]))
        AMS_TO_ADM[gid_level[vertices]] = level_id[vertices]

        nivel_0 = gid_0[levels==0]
        ID_global1 = nivel_0
        OP_AMS[nivel_0] = 0

        n1_adm = len(np.unique(level_id))

        ids_adm_nivel0 = level_id[nivel_0]
        IDs_ADM1 = ids_adm_nivel0

        m = sp.find(OP_AMS)
        l1=m[0]
        c1=m[1]
        d1=m[2]
        lines=ID_global1
        cols=IDs_ADM1
        data=np.repeat(1,len(lines))
        ID_ADM1 = AMS_TO_ADM[c1]

        lines = np.concatenate([lines,l1])
        cols = np.concatenate([cols,ID_ADM1])
        data = np.concatenate([data,d1])

        OP_ADM = sp.csc_matrix((data,(lines,cols)),shape=(len(gid_0),n1_adm))

        cols = gid_0
        lines = level_id
        data = np.ones(len(lines))
        OR_ADM = sp.csc_matrix((data,(lines,cols)),shape=(n1_adm,len(gid_0)))

        if self.get_correction_term:
            pcorr = np.zeros(len(gid_0))
            pcorr[levels>0] = _pcorr[levels>0]
        else:
            pcorr = np.array([False])

        return OP_ADM, OR_ADM, pcorr

    def tpfalize_matrix(self, Mat):
        adjs=self.elements_lv0['neig_internal_faces']
        idsc=self.data_impress['LEVEL_ID_1']
        f1=sp.find(Mat)

        lins, diags=f1[0][f1[0]==f1[1]], f1[2][f1[0]==f1[1]]
        control=f1[2]/diags[lins][f1[0]]

        lim=np.load('flying/neta_lim_dual.npy')[0]
        # lim=0.99*control.max()
        adjsc=np.vstack([f1[0][control<=lim],f1[1][control<=lim]]).T
        # adjsc=idsc[adjs]
        adjsc=adjsc[adjsc[:,0]!=adjsc[:,1]]

        ad0=np.concatenate([adjsc[:,0],adjsc[:,1]])
        ad1=np.concatenate([adjsc[:,1],adjsc[:,0]])
        ad0, ad1=ad0[ad0>ad1], ad1[ad0>ad1]
        ad0, ad1=np.unique(np.array([ad0, ad1]), axis=1)
        vals0=np.array(Mat[ad0,ad1])[0]
        vals1=np.array(Mat[ad1,ad0])[0]
        lines = np.concatenate([    ad0,    ad1,    ad0,    ad1])
        cols  = np.concatenate([    ad1,    ad0,    ad0,    ad1])
        data  = np.concatenate([  vals0,  vals1, -vals0, -vals1])
        Mat2   = csc_matrix((data, (lines, cols)),shape=Mat.shape)
        lines_diric=np.arange(Mat.shape[0])[(Mat.min(axis=1).toarray()==0).T[0]]
        Mat2[lines_diric,lines_diric]=1
        return(Mat2)


    def solve_multiscale_pressure(self, T: 'fine transmissibility matrix', b: 'fine source term'):

        T_adm = T.copy()
        b_adm = b.copy()

        n_levels = self.n_levels
        for i in range(1, n_levels):
            level = i
            OP_adm = self._data[self.adm_op_n + str(level)]
            OR_adm = self._data[self.adm_rest_n + str(level)]
            # OR_adm=OP_adm.T
            if self.get_correction_term:
                pcorr_adm = self._data[self.pcorr_n+str(level)]
                b_adm = OR_adm*b_adm - OR_adm*T_adm*pcorr_adm
            else:
                b_adm = OR_adm*b_adm

            T_adm = OR_adm*T_adm*OP_adm


        # CorrectionFunction(data_impress['keq_faces'])
        Cq = CorrectionFunction.solve(b)

        pms = OP_adm*self.solver.direct_solver(T_adm, b_adm)

        # pms = self.smoother_jacobi(OR_adm, T, OP_adm, b, print_errors=True)
        '''
        local_preconditioners=['bosma', 'jacobi', 'olav']
        global_multiscale_preconditioners=['galerkin']
        # compare_smoothers(OR_adm, T, OP_adm, b, global_multiscale_preconditioners, self.elements_lv0, self.data_impress)
        # import pdb; pdb.set_trace()
        ets=[]
        names=[]
        for l in local_preconditioners:
            for g in global_multiscale_preconditioners:
                et = self.test_smoothers(OR_adm, T, OP_adm, b, local_preconditioner=l,\
                global_multiscale_preconditioner=g, n_prec_levels=1)
                ets.append(np.array(et))
                names.append(l+'_'+g)
        plt.close('all')
        self.plot_graphics(ets, names)'''

        self.data_impress['pms'] = pms
        self.data_impress['pressure'] = pms
        # ##########################
        # p2 = self.solver.direct_solver(T, b)
        # self.data_impress['pressure'] = p2
        # self.data_impress['pms'] = p2
        # ##############################
        self.T = T

    def plot_graphics(self, ets, names):
        a=1
        norms=['ep_L2', 'ep_Linf', 'ev_L2', 'ev_Linf']
        for i in range(4): # 4 normas de erro serão plotadas
            plt.close('all')
            fig=plt.figure()

            plt.plot()
            tf=np.inf
            tc=np.inf
            ymin=np.inf
            ymax=0
            xmax=0
            for j in range(len(names)):
                tf2=ets[j][:,4][ets[j][:,4]>0].min()
                tc2=ets[j][:,5][ets[j][:,5]>0].min()
                if tf2<tf:
                    tf=tf2
                if tc2<tc:
                    tc=tc2
                xx=np.cumsum(ets[j][:,6]+ets[j][:,7])+ets[j][:,5].max()
                yy=ets[j][:,i]
                if yy.min()<ymin:
                    ymin=yy.min()
                if yy.max()>ymax and yy.max()<100:
                    ymax=yy.max()
                if xx.max()>xmax and yy.max()<100:
                    xmax=xx.max()
                if yy.max()<100:
                    plt.plot(xx,yy, label=names[j])
                plt.xlabel('Time (s)')
                plt.ylabel(norms[i])
                plt.legend()
            plt.plot(np.repeat(tf,2),np.array([ymin, ymax])) # plots time to solve Tf
            plt.plot(np.repeat(tc,2),np.array([ymin, ymax])) # plots time to solve Tc
            plt.xlim((0,xmax))
            plt.savefig('results_smoothers/'+norms[i]+'.png')

        import pdb; pdb.set_trace()

    def smoother_jacobi(self, R, T, P, b, print_errors=False):

        # try:
        pv=np.load('flying/pms.npy')
        # except:
        #     pv=self.data_impress['pms']
        # pv=np.zeros(len(b))
        dt=np.array(T[range(T.shape[0]),range(T.shape[0])])[0]
        dt1=1/dt
        lc=np.arange(len(dt))
        pl=csc_matrix((dt1,(lc,lc)),shape=T.shape)

        T2=T.copy()
        T2.setdiag(0)
        J=pl*T2
        ap=np.zeros((3,1))
        try:
            ji=np.load("results/jac_iterarion.npy")[0]
        except:
            np.save('results/jac_iterarion.npy',np.array([0]))
            ji=0
        if ji==0:
            pv=np.zeros(len(b))
            nite=10
        else:
            pv=np.load('flying/pms.npy')
            nite=1
        mant=np.inf

        for i in range(nite):
            pv12 = pv + (P*spsolve(csc_matrix(P.T*T*P),P.T.tocsc()*(b-T*pv)))
            for j in range(5):
                pv = pl*b-J*pv12

            if print_errors:
                if i==0:
                    pf = self.solver.direct_solver(T, b)
                ep_l2=np.linalg.norm(pf-pv)/np.linalg.norm(pf)
                adjs=self.elements_lv0['neig_internal_faces']
                internal_faces=self.elements_lv0['internal_faces']
                ts=self.data_impress['transmissibility'][internal_faces]
                hs=self.data_impress['u_normal'][internal_faces].max(axis=1)
                a0=adjs[:,0]
                a1=adjs[:,1]
                vf=ts*(pf[a0]-pf[a1])/hs
                vms=ts*(pv[a0]-pv[a1])/hs
                ev_l2=np.linalg.norm(vf-vms)/np.linalg.norm(vf)
                ap=np.hstack([ap,np.array([[ji+i+1],[ep_l2],[ev_l2]])])
                print(ev_l2,ep_l2)
                # if ev_l2>mant:
                #     print("divergiu")
                #     import pdb; pdb.set_trace()
                # else:
                #     mant=ev_l2
            np.save('results/jac_iterarion.npy',np.array([ji+i+1]))
        np.save('flying/pms.npy',pv)
        pv = pv + (P*spsolve(csc_matrix(R*T*P),R.tocsc()*(b-T*pv)))
        vms=ts*(pv[a0]-pv[a1])/hs
        ev_l2=np.linalg.norm(vf-vms)/np.linalg.norm(vf)
        ap=np.hstack([ap,np.array([[ji+i+2],[ep_l2],[ev_l2]])])
        # plt.plot(ap[0,1:],ap[1,1:])
        plt.yscale('log')
        plt.scatter(ap[0,1:],ap[2,1:])
        plt.savefig("results/notms.png")

        # if print_errors:
        #     try:
        #         jac = np.load('flying/jacobi_smooter.npy')
        #         i_ant = jac[0,-1]
        #         jac.hstack(jac,ap)
        #     except:
        #         jac = np.save('flying/jacobi_smooter.npy')
        #         import pdb; pdb.set_trace()

        return pv

    def test_smoothers(self, R, T, P, b, local_preconditioner, global_multiscale_preconditioner, n_prec_levels=1):
        neta_lim=self.elements_lv0['neta_lim']
        tc=time.time()
        if local_preconditioner=='jacobi':
            alpha_jacobi=0.97
            Ml, pl = self.construct_jacobi_preconsitioner(T)
        elif local_preconditioner=='bosma' or local_preconditioner=='olav':
            Ml=linalg.spilu(T,drop_tol=1e-4,fill_factor=1,permc_spec='NATURAL')
        elif local_preconditioner=='sor':
            Sf=self.get_Sor(T,n_prec_levels)
        tc=time.time()-tc
        if global_multiscale_preconditioner=='galerkin':
            Rg=P.T
        elif global_multiscale_preconditioner=='msfv':
            Rg=R

        Tcg=Rg*T*P #global T coarse

        tf=time.time()
        pf=spsolve(T,b)
        tf=time.time()-tf

        pv=P*spsolve(Rg*T*P,Rg*b)
        # pv=np.zeros_like(b)
        maxiter=20
        er_l2=np.inf
        cont=0
        errors=[]
        times=[]
        ep_l2, ep_linf, ev_l2, ev_linf = self.get_error_norms(pf, pv)
        times.append([0,0,0,0])
        errors.append([ep_l2, ep_linf, ev_l2, ev_linf])
        while er_l2>1e-4 and cont<maxiter:
            tg=time.time()
            pv12 = pv + 1.0*(P*spsolve(Tcg,Rg*(b-T*pv)))
            tg=time.time()-tg
            tl=time.time()
            if local_preconditioner=='jacobi':
                pv = self.iterate_jacobi(pv, pv12,b, pl, Ml, alpha_jacobi=alpha_jacobi)
            elif local_preconditioner=='bosma':
                for i in range(5):
                    pv=pv12+Ml.solve(b-T*pv12)
                    pv12=pv
            elif local_preconditioner=='olav':
                dk=b-T*pv12
                yrf=Ml.solve(dk)
                pv=pv12+P*spsolve(Tcg,Rg*(dk-T*yrf))+yrf
            elif local_preconditioner=='sor':
                pf3 = pv12 + linalg.spsolve(Sf,(b - T*pv12))
                pv = pf3 + linalg.spsolve(Sf,(b - T*pf3))
            tl=time.time()-tl
            er_l2=np.linalg.norm(b-T*pv)/np.linalg.norm(b)
            print(er_l2, 'erro relatitivo no resíduo')
            pms = pv + (P*spsolve(csc_matrix(Rg*T*P),Rg.tocsc()*(b-T*pv)))
            ep_l2, ep_linf, ev_l2, ev_linf = self.get_error_norms(pf, pms)
            errors.append([ep_l2, ep_linf, ev_l2, ev_linf])
            times.append([tf, tc, tg, tl])
            cont+=1
        errors=np.vstack(errors)
        times=np.vstack(times)
        tt=times[:,2].sum()+times[:,3].sum()
        et=np.hstack([errors, times])


        return et

    def get_error_norms(self, pf, pv):
        ep_l2=np.linalg.norm(pf-pv)/np.linalg.norm(pf)
        adjs=self.elements_lv0['neig_internal_faces']
        internal_faces=self.elements_lv0['internal_faces']
        ts=self.data_impress['transmissibility'][internal_faces]
        hs=self.data_impress['u_normal'][internal_faces].max(axis=1)
        a0=adjs[:,0]
        a1=adjs[:,1]
        vf=ts*(pf[a0]-pf[a1])/hs
        vms=ts*(pv[a0]-pv[a1])/hs
        ev_l2=np.linalg.norm(vf-vms)/np.linalg.norm(vf)
        ep_linf=abs(pf-pv).max()/abs(pf).max()
        ev_linf=abs(vf-vms).max()/abs(vf).max()
        return ep_l2, ep_linf, ev_l2, ev_linf

    def get_Sor(self, T, Tcg=None, n_prec_levels=1, Wf=2, Wc=2/3):
        lf, cf, df = find(T)
        l=lf>cf
        # u=lf<cf
        d=lf==cf
        L=csc_matrix((df[l],(lf[l],cf[l])),shape=T.shape)
        # U=csc_matrix((df[u],(lf[u],cf[u])),shape=T.shape)
        D=csc_matrix((df[d],(lf[d],cf[d])),shape=T.shape)
        S=D+Wf*L
        if n_prec_levels==2:
            lf, cf, df = find(Tcg)
            l=lf>cf
            # u=lf<cf
            d=lf==cf
            Lc=csc_matrix((df[l],(lf[l],c[l])),shape=T.shape)
            # Uc=csc_matrix((df[u],(lf[u],c[u])),shape=T.shape)
            Dc=csc_matrix((df[d],(lf[d],c[d])),shape=T.shape)
            Sc=Dc+Wc*Lc
            return S.tocsr(), Sc.tocsr()
        else:
            return S.tocsr()

    def construct_jacobi_preconsitioner(self, T):
        dt=np.array(T[range(T.shape[0]),range(T.shape[0])])[0]
        dt1=1/dt
        lc=np.arange(len(dt))
        pl=csc_matrix((dt1,(lc,lc)),shape=T.shape)
        T2=T.copy()
        T2.setdiag(0)
        J=pl*T2
        return J, pl

    def iterate_jacobi(self, pv, pv12,b, pl, J, alpha_jacobi):
        de_0=1
        de_ant=1
        dp=np.linalg.norm(pv-pv12)
        cj=0
        # while ((de_ant-dp)/(de_0-de_ant)<alpha_jacobi and cj<100) or cj<3:
        for i in range(20):
            de_0=de_ant
            de_ant=dp
            pv = pl*b-J*pv12
            # dp=np.linalg.norm(pv-pv12)
            pv12=pv
            cj+=1
            # print((dp-de_ant)/(de_ant-de_0))
        # print(cj)
        return pv

    def set_pms_flux_intersect_faces_dep0(self):

        levels = self.data_impress['LEVEL']
        faces_intersect_lv1 = np.unique(np.concatenate(self.ml_data['coarse_intersect_faces_level_'+str(1)]))
        neig_intersect_faces_lv1 = self.ml_data['neig_intersect_faces_level_'+str(1)]

        v0 = neig_intersect_faces_lv1
        n_volumes = len(self.elements_lv0['volumes'])
        pms = self.data_impress['pms']
        t_intersect_faces = self.data_impress['transmissibility'][faces_intersect_lv1]
        t0 = t_intersect_faces
        flux_grav_intersect_faces = self.data_impress['flux_grav_faces'][faces_intersect_lv1]

        ps0 = pms[v0[:, 0]]
        ps1 = pms[v0[:, 1]]
        flux_intersect_faces = -((ps1 - ps0) * t0 - flux_grav_intersect_faces)

        lines = np.concatenate([v0[:, 0], v0[:, 1]])
        cols = np.repeat(0, len(lines))
        data = np.concatenate([flux_intersect_faces, -flux_intersect_faces])
        flux_pms_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()

        n = len(self.data_impress['pms_flux_faces'])
        flux = np.zeros(n)
        flux[faces_intersect_lv1] = flux_intersect_faces

        self.data_impress['pms_flux_faces'] = flux
        self.data_impress['pms_flux_interfaces_volumes'] = flux_pms_volumes

    def set_pms_flux_intersect_faces(self):

        levels = self.data_impress['LEVEL']
        n_volumes = len(levels)
        flux_volumes = np.zeros(n_volumes)
        gid0 = self.data_impress['GID_0']
        transmissibility = self.data_impress['transmissibility']
        pms = self.data_impress['pms']
        neig_internal_faces = self.elements_lv0['neig_internal_faces']
        remaped_internal_faces = self.elements_lv0['remaped_internal_faces']
        flux_grav_faces = self.data_impress['flux_grav_faces']

        presc_flux_volumes = np.zeros(len(self.data_impress['pms_flux_interfaces_volumes']))
        pms_flux_faces = np.zeros(len(transmissibility))

        for i in range(self.n_levels):
            level=i+1
            all_gids_coarse = self.data_impress['GID_'+str(level)]
            all_intern_boundary_volumes = self.ml_data['internal_boundary_fine_volumes_level_'+str(level)]
            all_intersect_faces = self.ml_data['coarse_intersect_faces_level_'+str(level)]
            all_intern_faces = self.ml_data['coarse_internal_faces_level_'+str(level)]
            coarse_ids = self.ml_data['coarse_primal_id_level_'+str(level)]
            gids_level = np.unique(all_gids_coarse[levels==level])
            for gidc in gids_level:
                intersect_faces = all_intersect_faces[coarse_ids==gidc][0] # faces na interseccao
                neig_intersect_faces = neig_internal_faces[remaped_internal_faces[intersect_faces]]
                v0 = neig_intersect_faces
                intern_boundary_volumes = all_intern_boundary_volumes[coarse_ids==gidc][0] # volumes internos no contorno
                flux_grav_intersect_faces = flux_grav_faces[intersect_faces]

                pms0 = pms[neig_intersect_faces[:,0]]
                pms1 = pms[neig_intersect_faces[:,1]]
                t0 = transmissibility[intersect_faces]
                flux_intersect_faces = -((pms1 - pms0) * t0 - flux_grav_intersect_faces)
                pms_flux_faces[intersect_faces] = flux_intersect_faces

                lines = np.concatenate([v0[:, 0], v0[:, 1]])
                cols = np.repeat(0, len(lines))
                data = np.concatenate([flux_intersect_faces, -flux_intersect_faces])
                flux_pms_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
                presc_flux_volumes[intern_boundary_volumes] = flux_pms_volumes[intern_boundary_volumes]

        self.data_impress['pms_flux_interfaces_volumes'] = presc_flux_volumes
        self.data_impress['pms_flux_faces'] = pms_flux_faces

    def set_pcorr(self):

        _debug = data_loaded['_debug']
        presc_flux_volumes = self.data_impress['pms_flux_interfaces_volumes'].copy()
        flux_faces = self.data_impress['pms_flux_faces']
        levels = self.data_impress['LEVEL']
        n_volumes = len(levels)
        gid0 = self.data_impress['GID_0']
        transmissibility = self.data_impress['transmissibility']
        pms = self.data_impress['pms']
        neig_internal_faces = self.elements_lv0['neig_internal_faces']
        remaped_internal_faces = self.elements_lv0['remaped_internal_faces']
        flux_grav_faces = self.data_impress['flux_grav_faces']
        flux_grav_volumes = self.data_impress['flux_grav_volumes']

        pcorr = np.zeros(len(pms))
        flux_faces = np.zeros(len(transmissibility))
        flux_volumes = np.zeros(n_volumes)

        for i in range(self.n_levels):
            level=i+1
            all_gids_coarse = self.data_impress['GID_'+str(level)]
            all_local_ids_coarse = self.data_impress['COARSE_LOCAL_ID_'+str(level)]
            all_intern_boundary_volumes = self.ml_data['internal_boundary_fine_volumes_level_'+str(level)]
            all_intersect_faces = self.ml_data['coarse_intersect_faces_level_'+str(level)]
            all_intern_faces = self.ml_data['coarse_internal_faces_level_'+str(level)]
            all_faces = self.ml_data['coarse_faces_level_'+str(level)]
            all_fine_vertex = self.ml_data['fine_vertex_coarse_volumes_level_'+str(level)]
            coarse_ids = self.ml_data['coarse_primal_id_level_'+str(level)]
            gids_level = np.unique(all_gids_coarse[levels==level])
            for gidc in gids_level:
                all_local_faces = all_faces[coarse_ids==gidc][0]
                intersect_faces = all_intersect_faces[coarse_ids==gidc][0] # faces na interseccao
                intern_local_faces = all_intern_faces[coarse_ids==gidc][0] # faces internas
                neig_internal_local_faces = neig_internal_faces[remaped_internal_faces[intern_local_faces]]
                intern_boundary_volumes = all_intern_boundary_volumes[coarse_ids==gidc][0] # volumes internos no contorno
                vertex = all_fine_vertex[coarse_ids==gidc]
                pressure_vertex = pms[vertex]
                volumes = gid0[all_gids_coarse==gidc]

                pdb.set_trace()

                internal_volumes = np.setdiff1d(volumes, np.concatenate([intern_boundary_volumes, vertex]))
                flux_grav_internal_volumes = flux_grav_volumes[internal_volumes]

                neig_intersect_faces = neig_internal_faces[remaped_internal_faces[intersect_faces]]
                v0 = neig_intersect_faces
                transmissibility_intersect_faces = transmissibility[intersect_faces]
                t0 = transmissibility_intersect_faces
                try:

                    pms0 = pms[neig_intersect_faces[:,0]]
                    pms1 = pms[neig_intersect_faces[:,1]]
                except Exception as e:
                    pdb.set_trace()
                flux_grav_intersect_faces = flux_grav_faces[intersect_faces]
                flux_intersect_faces = -((pms1 - pms0) * t0 - flux_grav_intersect_faces)
                flux_faces[intersect_faces] = flux_intersect_faces

                lines = np.concatenate([v0[:, 0], v0[:, 1]])
                cols = np.repeat(0, len(lines))
                data = np.concatenate([flux_intersect_faces, -flux_intersect_faces])
                flux_pms_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
                # presc_flux_volumes[intern_boundary_volumes] = flux_pms_volumes[intern_boundary_volumes]
                values_q = flux_pms_volumes[intern_boundary_volumes]

                local_id_volumes = all_local_ids_coarse[volumes]
                local_neig_internal_local_faces = neig_internal_local_faces.copy()
                local_neig_internal_local_faces[:,0] = all_local_ids_coarse[neig_internal_local_faces[:,0]]
                local_neig_internal_local_faces[:,1] = all_local_ids_coarse[neig_internal_local_faces[:,1]]
                local_intern_boundary_volumes = all_local_ids_coarse[intern_boundary_volumes]
                t0 = transmissibility[intern_local_faces]
                local_vertex = all_local_ids_coarse[vertex]
                t0 = transmissibility[intern_local_faces]

                local_id_internal_volumes = all_local_ids_coarse[internal_volumes]
                indices_q2 = np.concatenate([local_intern_boundary_volumes, local_id_internal_volumes])
                values_q2 = np.concatenate([values_q, flux_grav_internal_volumes])
                # x = solve_local_local_problem(self.solver.direct_solver, local_neig_internal_local_faces, t0, local_id_volumes,
                #     local_vertex, pressure_vertex, local_intern_boundary_volumes, values_q)
                x = solve_local_local_problem(self.solver.direct_solver, local_neig_internal_local_faces, t0, local_id_volumes,
                    local_vertex, pressure_vertex, indices_q2, values_q2)

                pcorr[volumes] = x
                pressure_vols = self.data_impress['pms'][volumes]


                # neig_intersect_faces = neig_internal_faces[remaped_internal_faces[intersect_faces]]
                # transmissibility_intersect_faces = transmissibility[intersect_faces]
                # t0 = transmissibility_intersect_faces
                # pms0 = pms[neig_intersect_faces[:,0]]
                # pms1 = pms[neig_intersect_faces[:,1]]
                # flux_grav_intersect_faces = flux_grav_faces[intersect_faces]
                # flux_intersect_faces = -((pms1 - pms0) * t0 - flux_grav_intersect_faces)
                # flux_faces[intersect_faces] = flux_intersect_faces

                pcorr0 = pcorr[neig_internal_local_faces[:,0]]
                pcorr1 = pcorr[neig_internal_local_faces[:,1]]
                flux_grav_intern_faces = flux_grav_faces[intern_local_faces]
                t0 = transmissibility[intern_local_faces]
                flux_intern_faces = -((pcorr1 - pcorr0) * t0 - flux_grav_intern_faces)
                flux_faces[intern_local_faces] = flux_intern_faces

                v0 = neig_internal_local_faces

                lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
                cols = np.repeat(0, len(lines))
                data = np.array([flux_intern_faces, -flux_intern_faces]).flatten()
                flux_volumes_2 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
                flux_volumes_2[intern_boundary_volumes] += values_q
                flux_volumes[volumes] = flux_volumes_2[volumes]

                # ###
                # ## test
                # all_local_faces = np.setdiff1d(all_local_faces, self.elements_lv0['boundary_faces'])
                # neig_local_faces = neig_internal_faces[remaped_internal_faces[all_local_faces]]
                # v0 = neig_local_faces
                # test_flux_faces = flux_faces[all_local_faces]
                # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
                # cols = np.repeat(0, len(lines))
                # data = np.array([test_flux_faces, -test_flux_faces]).flatten()
                # flux_volumes_2 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
                # flux_volumes_2 = flux_volumes_2[volumes]
                # if not np.allclose(np.absolute(flux_volumes_2), np.absolute(flux_volumes[volumes])):
                #     import pdb; pdb.set_trace()
                # ##################

        volumes_fine = gid0[levels==0]
        intern_faces_volumes_fine = self.mesh.volumes.bridge_adjacencies(volumes_fine, 3, 2)
        intern_faces_volumes_fine = np.setdiff1d(intern_faces_volumes_fine, self.elements_lv0['boundary_faces'])
        neig_intern_faces_volumes_fine = neig_internal_faces[remaped_internal_faces[intern_faces_volumes_fine]]
        v0 = neig_intern_faces_volumes_fine

        pms0 = pms[neig_intern_faces_volumes_fine[:,0]]
        pms1 = pms[neig_intern_faces_volumes_fine[:,1]]
        t0 = transmissibility[intern_faces_volumes_fine]
        flux_grav_faces_volumes_fine = flux_grav_faces[intern_faces_volumes_fine]
        flux_intern_faces_volumes_fine = -((pms1 - pms0) * t0 - flux_grav_faces_volumes_fine)
        flux_faces[intern_faces_volumes_fine] = flux_intern_faces_volumes_fine

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.concatenate([flux_intern_faces_volumes_fine, -flux_intern_faces_volumes_fine])
        flux_volumes_2 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()

        flux_volumes[volumes_fine] = flux_volumes_2[volumes_fine]

        self.data_impress['pcorr'] = pcorr
        self.data_impress['flux_faces'] = flux_faces
        self.data_impress['flux_volumes'] = flux_volumes

        if _debug:
            #######################
            ## test
            v0 = neig_internal_faces
            internal_faces = self.elements_lv0['internal_faces']
            lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
            cols = np.repeat(0, len(lines))
            data = np.array([flux_faces[internal_faces], -flux_faces[internal_faces]]).flatten()
            flux_volumes_2 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
            # if not np.allclose(flux_volumes_2, flux_volumes):
            #     raise ValueError('Diferenca entre fluxo pms local e global')
            self.data_impress['flux_volumes_test'] = flux_volumes_2

    def get_nt_process(self):

        levels = self.data_impress['LEVEL']
        gids = self.data_impress['GID_0']
        nt_process = 0

        for i in range(1, self.n_levels):
            level = i
            gids_lv = self.data_impress['GID_' + str(level)]
            gids_eng = np.unique(gids_lv[levels==level])
            nt_process += len(gids_eng)

        return nt_process

    def get_lists_objects(self, infos):

        levels = self.data_impress['LEVEL']
        pms = self.data_impress['pms']
        gid0 = self.data_impress['GID_0']
        presc_flux_volumes = self.data_impress['pms_flux_interfaces_volumes']

        all_list_volumes = []
        all_list_indices_d = []
        all_list_values_d = []
        all_list_indices_n = []
        all_list_values_n = []
        all_list_faces = []
        all_list_internal_faces = []
        all_list_intersect_faces = []
        cont0 = 0
        cont = 0

        for i in range(self.n_levels):
            level=i+1
            st_level = str(level)
            all_gids_coarse = self.data_impress['GID_'+ st_level]
            all_local_ids_coarse = self.data_impress['COARSE_LOCAL_ID_'+ st_level]
            all_intern_boundary_volumes = self.ml_data['internal_boundary_fine_volumes_level_'+ st_level]
            all_intersect_faces = self.ml_data['coarse_intersect_faces_level_'+ st_level]
            all_intern_faces = self.ml_data['coarse_internal_faces_level_'+ st_level]
            all_faces = self.ml_data['coarse_faces_level_'+ st_level]
            all_fine_vertex = self.ml_data['fine_vertex_coarse_volumes_level_'+ st_level]
            coarse_ids = self.ml_data['coarse_primal_id_level_'+ st_level]
            gids_level = np.unique(all_gids_coarse[levels==level])
            for gidc in gids_level:
                intersect_faces = all_intersect_faces[coarse_ids==gidc][0] # faces na interseccao
                intern_local_faces = all_intern_faces[coarse_ids==gidc][0] # faces internas
                faces = all_faces[coarse_ids==gidc][0] # faces do volume
                intern_boundary_volumes = all_intern_boundary_volumes[coarse_ids==gidc][0] # volumes internos no contorno
                vertex = all_fine_vertex[coarse_ids==gidc]
                pressure_vertex = pms[vertex]
                volumes = gid0[all_gids_coarse==gidc]
                values_q = presc_flux_volumes[intern_boundary_volumes]

                if cont0 < self.n_workers and cont0 >= 0:
                    all_list_volumes.append([volumes])
                    all_list_indices_d.append([vertex])
                    all_list_values_d.append([pressure_vertex])
                    all_list_indices_n.append([intern_boundary_volumes])
                    all_list_values_n.append([values_q])
                    all_list_faces.append([faces])
                    all_list_internal_faces.append([intern_local_faces])
                    all_list_intersect_faces.append([intersect_faces])
                    cont0 += 1
                    if cont0 == self.n_workers:
                        cont0 = -1
                else:
                    all_list_volumes[cont].append(volumes)
                    all_list_indices_d[cont].append(vertex)
                    all_list_values_d[cont].append(pressure_vertex)
                    all_list_indices_n[cont].append(intern_boundary_volumes)
                    all_list_values_n[cont].append(values_q)
                    all_list_faces[cont].append(faces)
                    all_list_internal_faces[cont].append(intern_local_faces)
                    all_list_intersect_faces[cont].append(intersect_faces)
                    cont += 1
                    if cont == self.n_workers:
                        cont = 0

        list_objects = []
        for i in range(self.n_workers):
            list_objects.append(
                LocalSolution(
                    all_list_volumes[i],
                    all_list_indices_d[i],
                    all_list_values_d[i],
                    all_list_indices_n[i],
                    all_list_values_n[i],
                    all_list_faces[i],
                    all_list_internal_faces[i],
                    all_list_intersect_faces[i],
                    infos.copy()
                )
            )

        return list_objects

    def set_paralel_pcorr(self):
        self.set_pms_flux_intersect_faces()
        transmissibility = self.data_impress['transmissibility']
        presc_flux_volumes = self.data_impress['pms_flux_interfaces_volumes']
        levels = self.data_impress['LEVEL']
        gids = self.data_impress['GID_0']
        pms = self.data_impress['pms']
        g_neig_internal_faces = self.elements_lv0['neig_internal_faces']
        remaped_internal_faces = self.elements_lv0['remaped_internal_faces']
        g_flux_grav_faces = self.data_impress['flux_grav_faces']
        g_faces = self.elements_lv0['faces']
        T = self.T
        solver = self.solver.direct_solver
        n_volumes = len(gids)

        _pcorr = np.zeros(n_volumes)
        _flux_faces = np.zeros(len(transmissibility))
        _flux_volumes = np.zeros(n_volumes)

        # nt_process = self.get_nt_process()
        # self.n_workers = nt_process
        # self.n_workers = 1

        # m = mp.Manager()
        # qvolumes = m.Queue()
        # qfaces = m.Queue()
        # qinfos = m.Queue()
        # lock = mp.Lock()

        infos = InfosForProcess(T, pms, g_flux_grav_faces, gids, g_faces, g_neig_internal_faces,
            remaped_internal_faces, solver)

        # for i in range(self.n_workers):
        #     qinfos.put(infos.copy())
        # tt = qinfos.qsize()

        list_objects = self.get_lists_objects(infos)

        def f(local_solution_obj, w2m):
            local_solution_obj.run(w2m)

        master2worker = [mp.Pipe() for _ in range(self.n_workers)]
        m2w, w2m = list(zip(*master2worker))
        procs = [mp.Process(target=f, args=[obj, comm]) for obj, comm in zip(list_objects, w2m)]

        for proc in procs:
            proc.start()

        for comm in m2w:
            msg = comm.recv()
            for resp in msg:
                resp_vols = resp[0]
                resp_faces = resp[1]
                _pcorr[resp_vols['volumes']] = resp_vols['pcorr']
                _flux_volumes[resp_vols['volumes']] = resp_vols['flux_volumes']
                _flux_faces[resp_faces['faces']] = resp_faces['flux_faces']

        for proc in procs:
            proc.join()

        # while not qvolumes.empty():
        #     resp = qvolumes.get()
        #     _pcorr[resp['volumes']] = resp['pcorr']
        #     _flux_volumes[resp['volumes']] = resp['flux_volumes']
        #
        # while not qfaces.empty():
        #     resp = qfaces.get()
        #     _flux_faces[resp['faces']] = resp['flux_faces']

        # for result in results:
        #     for resp in result:
        #         resp_vols = resp[0]
        #         resp_faces = resp[1]
        #         _pcorr[resp_vols['volumes']] = resp_vols['pcorr']
        #         _flux_volumes[resp_vols['volumes']] = resp_vols['flux_volumes']
        #         _flux_faces[resp_faces['faces']] = resp_faces['flux_faces']

        gid0 = gids
        volumes_fine = gid0[levels==0]
        intern_faces_volumes_fine = self.mesh.volumes.bridge_adjacencies(volumes_fine, 3, 2)
        intern_faces_volumes_fine = np.setdiff1d(intern_faces_volumes_fine, self.elements_lv0['boundary_faces'])
        neig_intern_faces_volumes_fine = g_neig_internal_faces[remaped_internal_faces[intern_faces_volumes_fine]]
        v0 = neig_intern_faces_volumes_fine

        pms0 = pms[neig_intern_faces_volumes_fine[:,0]]
        pms1 = pms[neig_intern_faces_volumes_fine[:,1]]
        t0 = transmissibility[intern_faces_volumes_fine]
        flux_grav_faces_volumes_fine = g_flux_grav_faces[intern_faces_volumes_fine]
        flux_intern_faces_volumes_fine = -((pms1 - pms0) * t0 - flux_grav_faces_volumes_fine)
        _flux_faces[intern_faces_volumes_fine] = flux_intern_faces_volumes_fine

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.concatenate([flux_intern_faces_volumes_fine, -flux_intern_faces_volumes_fine])
        flux_volumes_2 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
        _flux_volumes[volumes_fine] = flux_volumes_2[volumes_fine]

        self.data_impress['pcorr'] = _pcorr
        self.data_impress['flux_volumes'] = _flux_volumes
        self.data_impress['flux_faces'] = _flux_faces

    def set_saturation_level(self):

        levels = self.data_impress['LEVEL'].copy()
        gid1 = self.data_impress['GID_1']
        gid0 = self.data_impress['GID_0']
        level_0_ini = set(gid0[levels==0])
        saturation = self.data_impress['saturation']
        all_wells = set(self.all_wells_ids)
        gids_lv1_sat = set()
        gidsc = np.unique(gid1)
        for gidc in gidsc:
            gids0 = gid0[gid1==gidc]
            if set(gids0) & all_wells:
                continue
            sats_local = saturation[gids0]
            dif = sats_local.max() - sats_local.min()
            if dif >= self.delta_sat_max:
                levels[gids0] = 0
                gids_lv1_sat.add(gidc)

        cids_neigh = self.ml_data['coarse_id_neig_face_level_'+str(1)]
        cids_level = self.ml_data['coarse_primal_id_level_'+str(1)]

        for gidc in gids_lv1_sat:
            vizs = cids_neigh[cids_level==gidc]
            for viz in vizs:
                if set([viz]) & gids_lv1_sat:
                    continue
                gids0 = gid0[gid1==gidc]
                if set(gids0) & level_0_ini:
                    continue
                levels[gids0] = np.repeat(1, len(gids0))

        self.data_impress['LEVEL'] = levels.copy()

    def set_initial_mesh(self, mlo, T, b):


        M = self.mesh

        iterar_mono = file_adm_mesh_def['iterar_mono']
        refinar_nv2 = file_adm_mesh_def['refinar_nv2']
        imprimir_a_cada_iteracao = file_adm_mesh_def['imprimir_a_cada_iteracao']
        rel_v2 = file_adm_mesh_def['rel_v2']
        TOL = file_adm_mesh_def['TOL']
        tol_n2 = file_adm_mesh_def['tol_n2']
        Ni = file_adm_mesh_def['Ni']
        calc_tpfa = file_adm_mesh_def['calc_tpfa']
        load_tpfa = file_adm_mesh_def['load_tpfa']
        _dev = file_adm_mesh_def['_dev']
        name = 'flying/SOL_TPFA.npy'
        nfine_vols = len(self.data_impress['LEVEL'])
        GID_0 = self.data_impress['GID_0']
        GID_1 = self.data_impress['GID_1']
        DUAL_1 = self.data_impress['DUAL_1']
        solver = file_adm_mesh_def['solver']
        load_adm_levels = file_adm_mesh_def['load_adm_levels']
        set_initial_mesh = file_adm_mesh_def['set_initial_mesh']

        if load_adm_levels:
            return 0

        if not set_initial_mesh:
            self.restart_levels()
            self.set_level_wells()
            self.set_adm_mesh()
            return 0

        if calc_tpfa:
            SOL_TPFA = self.solver.direct_solver(T, b)
            print("\nresolveu TPFA\n")
            np.save(name, SOL_TPFA)
        elif load_tpfa:
            try:
                SOL_TPFA = np.load(name)
            except:
                raise FileNotFoundError('O aqruivo {} nao existe'.format(name))

        if solver == 'direct':
            solver = linalg.spsolve

        self.restart_levels()
        self.set_level_wells()
        self.set_adm_mesh()

        multilevel_meshes = []

        active_nodes = []
        perro = []
        erro = []

        Nmax = tol_n2*nfine_vols
        finos = self.all_wells_ids.copy()
        primal_finos = np.unique(GID_1[finos])
        pfins = primal_finos
        vertices = GID_0[DUAL_1==3]
        primal_id_vertices = GID_1[vertices]
        dt = [('vertices', np.dtype(int)), ('primal_vertices', np.dtype(int))]
        structured_array = np.zeros(len(vertices), dtype=dt)
        structured_array['vertices'] = vertices
        structured_array['primal_vertices'] = primal_id_vertices
        structured_array = np.sort(structured_array, order='primal_vertices')
        vertices = structured_array['vertices']
        primal_id_vertices = structured_array['primal_vertices']

        nr = int(tol_n2*(len(vertices)-len(primal_finos))/(Ni))
        n1 = self.data_impress['LEVEL_ID_1'].max() + 1
        n2 = self.data_impress['LEVEL_ID_2'].max() + 1

        accum_levels = []

        pseudo_erro=np.repeat(TOL+1,2) #iniciou pseudo_erro
        t0=time.time()
        cont=0
        pos_new_inter=[]
        interm=np.array([])
        continuar = True


        while (pseudo_erro.max()>TOL and n2<Nmax and iterar_mono and continuar) or cont==0:

            if cont>0:

                levels = self.data_impress['LEVEL'].copy()
                # import pdb; pdb.set_trace()
                n1_ant = self.data_impress['LEVEL_ID_1'].max() + 1
                n2_ant = self.data_impress['LEVEL_ID_2'].max() + 1

                lim=np.sort(psr)[len(psr)-nr-1]
                positions=np.where(psr>lim)[0]
                nv_verts=levels[vertices]
                nv_positions=nv_verts[positions]
                pos_new_fines=positions[nv_positions==1]
                pos_new_inter=positions[nv_positions==2]

                interm=np.concatenate([interm,np.array(vertices)[pos_new_inter]]).astype(np.int)
                finos=np.concatenate([finos,np.array(vertices)[pos_new_fines]]).astype(np.int)

                primal_id_interm = np.unique(GID_1[interm])
                interm = np.concatenate([GID_0[GID_1==k] for k in primal_id_interm])
                primal_id_finos = np.unique(GID_1[finos])
                finos = np.concatenate([GID_0[GID_1==k] for k in primal_id_finos])
                pfins=np.unique(GID_1[finos])
                self.restart_levels()
                levels = self.data_impress['LEVEL'].copy()
                levels[finos] = np.zeros(len(finos), dtype=int)
                levels[interm] = np.ones(len(interm), dtype=int)
                self.data_impress['LEVEL'] = levels.copy()
                self.set_adm_mesh()
                n1 = self.data_impress['LEVEL_ID_1'].max() + 1
                n2 = self.data_impress['LEVEL_ID_2'].max() + 1

                if n1 == n1_ant and n2 == n2_ant:
                    continuar = False

                if _dev:
                    print('\n',n1,n2,'n1 e n2\n')

            self.organize_ops_adm(mlo['prolongation_level_1'],
                                  mlo['restriction_level_1'],
                                  1)

            OP_ADM = self._data[self.adm_op_n + str(1)]
            OR_ADM = self._data[self.adm_rest_n + str(1)]

            if (len(pos_new_inter)>0 or cont==0) and refinar_nv2:
                self.organize_ops_adm(mlo['prolongation_level_2'],
                                      mlo['restriction_level_2'],
                                      2)

                OP_ADM_2 = self._data[self.adm_op_n + str(2)]
                OR_ADM_2 = self._data[self.adm_rest_n + str(2)]

                SOL_ADM=solver(OR_ADM_2*OR_ADM*T*OP_ADM*OP_ADM_2,OR_ADM_2*OR_ADM*b)
                SOL_ADM_fina=OP_ADM*OP_ADM_2*SOL_ADM
            else:
                SOL_ADM=solver(OR_ADM*T*OP_ADM,OR_ADM*b)
                SOL_ADM_fina=OP_ADM*SOL_ADM
            self.data_impress['pressure'] = SOL_ADM_fina
            x0=Jacobi(SOL_ADM_fina, T, b)
            pseudo_erro=abs((SOL_ADM_fina-x0))

            if calc_tpfa or load_tpfa:
                erro.append(abs((SOL_TPFA-SOL_ADM_fina)/SOL_TPFA).max())
            else:
                erro.append(abs(pseudo_erro/x0).max())
                SOL_TPFA=x0
            OR_AMS = mlo['restriction_level_1']
            psr=(OR_AMS*abs(pseudo_erro))
            psr[pfins]=0

            perro.append(abs((SOL_ADM_fina-x0)/x0).max())
            active_nodes.append(n2/nfine_vols)

            if imprimir_a_cada_iteracao:
                # M1.mb.tag_set_data(Pseudo_ERRO_tag,M1.all_volumes,abs(pseudo_erro/x0)[GIDs])
                #
                # M1.mb.tag_set_data(ERRO_tag,M1.all_volumes,abs((SOL_ADM_fina-SOL_TPFA)/SOL_TPFA)[GIDs])
                # M1.mb.tag_set_data(P_ADM_tag,M1.all_volumes,SOL_ADM_fina[GIDs])
                # M1.mb.tag_set_data(P_TPFA_tag,M1.all_volumes,SOL_TPFA[GIDs])
                # ext_vtk = 'testes_MAD'  + str(cont) + '.vtk'
                # M1.mb.write_file(ext_vtk,[av])
                self.data_impress.update_variables_to_mesh(['LEVEL', 'pressure'])
                M.core.print(folder='results', file='test'+ str(cont), extension='.vtk', config_input='input_cards/print_settings0.yml')
            cont+=1

            accum_levels.append(self.data_impress['LEVEL'].copy())


        plt.plot(active_nodes,perro, marker='o')
        plt.yscale('log')
        plt.savefig('results/initial_adm_mesh/hist.png')

        n = int(input('\nQual a malha adm que deseja utilizar?\nDigite o numero da iteracao.\n'))

        self.data_impress['INITIAL_LEVEL'] = accum_levels[n]

        self.data_impress.update_variables_to_mesh()
        self.data_impress.export_to_npz()

    def so_nv1():
        doc = "The so_nv1 property."
        def fget(self):
            return self._so_nv1
        def fset(self, value):
            self._so_nv1 = value
            if self._so_nv1:
                self.n_levels = 2
        def fdel(self):
            del self._so_nv1
        return locals()
    so_nv1 = property(**so_nv1())

    def print_test(self):
        self.data_impress.update_variables_to_mesh()
