import scipy.sparse as sp
import numpy as np
def monotonize_adm(mlo, T, neta_lim):
    T=T.copy()
    levels=np.ones(T.shape[0],dtype=int)
    op_ams = mlo['prolongation_level_1']
    or_ams = mlo['restriction_level_1']
    Tc=or_ams*T*op_ams
    lcd=sp.find(Tc)
    lc=lcd[0]
    cc=lcd[1]
    dc=lcd[2]
    diags=lc==cc
    off_diag=lc!=cc
    val_diag=dc[diags]
    lcc_diag=lc[diags]
    val_off_diag=dc[off_diag]

    lc_off_diag=lc[off_diag]
    val_diag_off_diag=val_diag[lcc_diag][lc_off_diag]
    netas=val_off_diag/val_diag_off_diag
    netasp=netas>neta_lim
    netas=netas[netasp]
    I=lc[off_diag][netasp]
    J=cc[off_diag][netasp]
    gids_to_monotonize=[]
    if len(netas)>0:
        print(netas.max(), "ams")
        # import pdb; pdb.set_trace()
        if netas.max()>neta_lim:
            for i in range(len(I)):
                tpJ=T*op_ams[:,J[i]]
                rItpJ=sp.csr_matrix.multiply(or_ams[I[i],:],tpJ.T)
                diagI=Tc[I[i],I[i]]
                try:
                    vec=-np.array(rItpJ[rItpJ<0])
                    if vec.shape[1]>0:
                        srItpJ=-np.sort(vec[0])
                    else:
                        srItpJ=np.array([])
                except:
                    import pdb; pdb.set_trace()
                cs=np.cumsum(srItpJ)
                val_lim=srItpJ[cs<diagI*neta_lim]
                if len(val_lim)>0:
                    val_lim=val_lim.max()
                else:
                    val_lim=-np.inf
                gids_to_monotonize.append(sp.find(rItpJ<val_lim)[1])
                # gids_to_monotonize.append(sp.find(rItpJ<0)[1])
                # import pdb; pdb.set_trace()
    if len(gids_to_monotonize)>1:
        gids_to_monotonize=np.concatenate(gids_to_monotonize)
    else:
        gids_to_monotonize=np.array(gids_to_monotonize).astype(int)
    return gids_to_monotonize

def verify_monotonize_adm(or_ams, T, op_ams, neta_lim,ids_ams_1):
    Tc=or_ams*T*op_ams
    lcd=sp.find(Tc)
    lc=lcd[0]
    cc=lcd[1]
    dc=lcd[2]
    diags=lc==cc
    off_diag=lc!=cc
    val_diag=dc[diags]
    lcc_diag=lc[diags]
    val_off_diag=dc[off_diag]

    lc_off_diag=lc[off_diag]
    val_diag_off_diag=val_diag[lcc_diag][lc_off_diag]
    netas=val_off_diag/val_diag_off_diag
    netasp=netas>10
    netas=netas[netasp]
    I=lc[off_diag][netasp]
    J=cc[off_diag][netasp]
    # tpJ=T*op_ams[:,J[1]]
    # rItpJ=sp.csr_matrix.multiply(or_ams[I[1],:],tpJ.T)


    if len(netas)>0:
        if netas.max()>10:
            print(netas.max(),'adm')
            print(ids_ams_1[I],ids_ams_1[J],netas)
        # import pdb; pdb.set_trace()
