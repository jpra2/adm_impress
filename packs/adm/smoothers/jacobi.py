import time
def jacobi_smoother(R, T, P, b, global_multiscale_preconditioner,get_error_norms):
    neta_lim=elements_lv0['neta_lim']
    tc=time.time()
    alpha_jacobi=0.97
    Ml, pl = construct_jacobi_preconsitioner(T)
    tc=time.time()-tc
    if global_multiscale_preconditioner=='galerkin':
        Rg=P.T
    elif global_multiscale_preconditioner=='msfv':
        Rg=R

    Tcg=Rg*T*P #global T coarse

    tf=time.time()
    pf=spsolve(T,b)
    tf=time.time()-tf
    pv=P*spsolve(R*T*P,R*b)
    maxiter=50
    ep_l2=np.inf
    cont=0
    errors=[]
    times=[]
    ep_l2, ep_linf, ev_l2, ev_linf = get_error_norms(pf, pv, elements_lv0, data_impress)
    times.append([0,0,0,0])
    errors.append([ep_l2, ep_linf, ev_l2, ev_linf])
    while ep_l2>0.01 and cont<maxiter:
        tg=time.time()
        pv12 = pv + 1.0*(P*spsolve(Tcg,Rg*(b-T*pv)))
        tg=time.time()-tg
        tl=time.time()
        pv = iterate_jacobi(pv, pv12,b, pl, Ml, alpha_jacobi=alpha_jacobi)
        tl=time.time()-tl
        pms = pv + (P*spsolve(csc_matrix(R*T*P),R.tocsc()*(b-T*pv)))
        ep_l2, ep_linf, ev_l2, ev_linf = get_error_norms(pf, pms, elements_lv0, data_impress)
        errors.append([ep_l2, ep_linf, ev_l2, ev_linf])
        times.append([tf, tc, tg, tl])
        cont+=1
    errors=np.vstack(errors)
    times=np.vstack(times)
    tt=times[:,2].sum()+times[:,3].sum()
    et=np.hstack([errors, times])
    return et

def construct_jacobi_preconsitioner(T):
    dt=np.array(T[range(T.shape[0]),range(T.shape[0])])[0]
    dt1=1/dt
    lc=np.arange(len(dt))
    pl=csc_matrix((dt1,(lc,lc)),shape=T.shape)
    T2=T.copy()
    T2.setdiag(0)
    J=pl*T2
    return J, pl
