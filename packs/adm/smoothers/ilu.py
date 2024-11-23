import time
def ilu_smoother(R, T, P, b, global_multiscale_preconditioner, get_error_norms, elements_lv0, data_impress):
      neta_lim=elements_lv0['neta_lim']
      tc=time.time()
      Ml=linalg.spilu(T,drop_tol=1e-4,fill_factor=1,permc_spec='NATURAL')
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
          tl=time.time()
          pv=pv12+Ml.solve(b-T*pv12)
          pv12=pv
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
