import scipy.sparse as sp
import numpy as np
def tpfalize(M,mlo,data_impress):
    int_faces=M.faces.internal
    adjs=M.faces.bridge_adjacencies(int_faces,2,3)
    ad0=adjs[:,0]
    ad1=adjs[:,1]

    op=mlo['prolongation_level_1']
    id_ams1=data_impress['GID_1']
    gid0=data_impress['GID_0']
    a0=id_ams1[ad0]
    a1=id_ams1[ad1]
    nvols=len(gid0)
    lines=np.concatenate([a0,a1,a0,a1])
    cols=np.concatenate([a1,a0,a0,a1])
    data=np.ones(len(lines))
    mc=sp.csc_matrix((data,(lines,cols)),shape=(nvols,nvols))
    fmc=sp.find(mc)
    ll=fmc[0]
    cc=fmc[1]
    lo=[]
    co=[]
    for i in np.unique(ll):
        cols=cc[ll==i]
        for c in cols:
            ad_lines=gid0[id_ams1==c]
            lo.append(ad_lines)
            co.append(np.repeat(i,len(ad_lines)))
    lo=np.concatenate(lo)
    co=np.concatenate(co)
    do=np.ones(len(lo))
    oo=sp.csc_matrix((do,(lo,co)),shape=(gid0.max()+1,id_ams1.max()+1))
    ot=sp.csr_matrix.multiply(oo,op)

    diag=sp.csc_matrix((1/np.array(ot.sum(axis=1)).T[0],(np.arange(nvols),np.arange(nvols))),shape=(nvols,nvols))
    op=diag*ot
    mlo['prolongation_level_1']=op
    return mlo
