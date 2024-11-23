import matplotlib.pyplot as plt
import numpy as np
from packs.adm.smoothers.ilu import ilu_smoother
from packs.adm.smoothers.jacobi import jacobi_smoother

def get_error_norms(pf, pv, elements_lv0, data_impress):
    ep_l2=np.linalg.norm(pf-pv)/np.linalg.norm(pf)
    adjs=elements_lv0['neig_internal_faces']
    internal_faces=elements_lv0['internal_faces']
    ts=data_impress['transmissibility'][internal_faces]
    hs=data_impress['u_normal'][internal_faces].max(axis=1)
    a0=adjs[:,0]
    a1=adjs[:,1]
    vf=ts*(pf[a0]-pf[a1])/hs
    vms=ts*(pv[a0]-pv[a1])/hs
    ev_l2=np.linalg.norm(vf-vms)/np.linalg.norm(vf)
    ep_linf=abs(pf-pv).max()/abs(pf).max()
    ev_linf=abs(vf-vms).max()/abs(vf).max()
    return ep_l2, ep_linf, ev_l2, ev_linf

def plot_graphics(ets, names):
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

def compare_smoothers(R, T, P, b, global_multiscale_preconditioners, elements_lv0, data_impress):
    ets=[]
    names=[]
    local_preconditioners = ilu_smoother, jacobi_smoother
    for smoother in local_preconditioners:
        for g in global_multiscale_preconditioners:
            et = smoother(R, T, P, b, g, get_error_norms, elements_lv0, data_impress)
            ets.append(np.array(et))
            names.append(l+'_'+g)
    plt.close('all')
    self.plot_graphics(ets, names)
