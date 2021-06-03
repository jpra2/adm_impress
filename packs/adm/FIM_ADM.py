from implicit_impress.jacobian.impress_assembly import assembly
import numpy as np
import scipy.sparse as sp
@profile
def newton_iteration_ADM(vec,ffv, vff, mlo, adm_method, F_Jacobian, Ts, adjs, data_impress, time_step, wells, rel_tol=1e-3):
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
        delta_sat_max=0.03
        adm_method.set_saturation_level_homogeneo(delta_sat_max)
        adm_method.set_adm_mesh_non_nested(np.arange(len(data_impress['swn1s']))[adm_method.data_impress['LEVEL']==0])
        OR_ADM, OP_ADM, _ = adm_method.get_OR_and_OP(mlo['prolongation_lcd_level_1'])
        # OP_ADM = adm_method._data['adm_prolongation_level_1']
        # OR_ADM = adm_method._data['adm_restriction_level_1']
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
    sats=get_sat_averager(sats,data_impress, vff,ffv, vec)
    # sats=OR_ADM.T*(OR_ADM*sats/np.array(OR_ADM.sum(axis=1)).T[0])
    data_impress['swns']=sats.copy()
    data_impress['swns'][sats>1]=1.0
    data_impress['swns'][sats<0]=0.0
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

def ADM_solver(J, q, R, P):
    sol=sp.linalg.spsolve(R*J*P,R*q)
    return sol

def get_sat_averager(sats, data_impress, vff,ffv, vec):
    gids0=data_impress['GID_0']
    gids1=data_impress['GID_1']
    level=data_impress['LEVEL']
    gids_adm_c=data_impress['LEVEL_ID_1']
    for i in np.unique(data_impress['LEVEL_ID_1'][level==1]):
        vols=gids0[gids_adm_c==i]
        faces=np.unique(np.concatenate(vff[vols]))
        faces=faces[vec[faces]==1]
        adjs=np.vstack(ffv[faces])
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

def newton_iteration_finescale(F_Jacobian, Ts, adjs, data_impress, time_step, wells, rel_tol=1e-3):
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
