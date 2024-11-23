import numpy as np
from scipy.sparse import csc_matrix, find, csgraph
import time

def get_coupled_dual_volumes(mlo,T, M, data_impress,neta_lim=0.0, ind=0):
    OP_AMS=mlo['prolongation_level_1']
    OR_AMS=mlo['restriction_level_1']
    Tc=OR_AMS*T*OP_AMS
    Tc2=Tc.copy()
    Tc2.setdiag(0)
    DTc=1/np.array(Tc[range(Tc.shape[0]),range(Tc.shape[0])])[0]
    # DTc*=np.log10(M.maxs)
    if (DTc>0).sum()>0:# and abs(Tc[DTc>0].sum())<0.01:
        print((DTc>0).sum(),"diagonais positivas !!!!!!!!!!!")
        DTc[DTc>0]=-abs(DTc).max()

    lines=np.arange(Tc.shape[0])
    dia=csc_matrix((DTc,(lines,lines)),shape=Tc.shape)
    netas=dia*Tc2
    fn=find(netas)
    superates_tol=fn[2]>neta_lim

    nsp=fn[2][superates_tol]
    i=fn[1][superates_tol]
    j=fn[0][superates_tol]

    ####################### new_feature

    maxs, ads = M.edges_and_vals
    beta_lim=np.load('flying/beta_lim_dual.npy')[0]
    bigs=maxs>beta_lim
    nsp=np.concatenate([nsp, np.zeros_like(ads[:,0][bigs])])
    i=np.concatenate([i,ads[:,0][bigs]])
    j=np.concatenate([j,ads[:,1][bigs]])

    ######################

    internal_faces=M.faces.internal
    adjs=M.faces.bridge_adjacencies(internal_faces,2,3)
    adjs0=adjs[:,0]
    adjs1=adjs[:,1]
    ii=data_impress['GID_1'][adjs0]
    jj=data_impress['GID_1'][adjs1]
    positives=fn[2]>0.0
    nsp_all=fn[2][positives]
    i_all=fn[1][positives]
    j_all=fn[0][positives]
    for k in range(len(nsp_all)):
        _non = (ii==i_all[k]) & (jj==j_all[k]) | (ii==j_all[k]) & (jj==i_all[k])
        ad0=adjs0[_non]
        ad1=adjs1[_non]
        neta_p=nsp_all[k]
        setar=np.concatenate([ad0,ad1])
        vals_set=data_impress["non_physical_value_"+str(ind)][setar[data_impress["DUAL_1"][setar]==2]] #Evita que I, J sobreponha J, I
        try:
            neta_p=max(neta_p,vals_set.max())
        except:
            neta_p=max(neta_p,vals_set)
        data_impress["non_physical_value_"+str(ind)][setar]=np.repeat(neta_p,len(setar))
    global dual_volumes
    dual_volumes = np.array(M.multilevel_data['dual_structure_level_1'])

    # dual_volumes = [dd['volumes'] for dd in dual_structure]

    dual_lines = [np.repeat(i,len(dual_volumes[i])) for i in range(len(dual_volumes))]
    # for dvv in range(len(dual_volumes)): #id_dual
    #     data_impress["perm_z"][dual_volumes[dvv]]=dual_lines[dvv]
    dvs=np.concatenate(dual_volumes)

    pvs=data_impress['GID_1'][dvs]
    dls=np.concatenate(dual_lines)
    data=np.repeat(1,len(pvs))
    dp=csc_matrix((data,(dls,pvs)),shape=(dls.max()+1, pvs.max()+1))
    dp[dp>1]=1
    cds=[]
    for k in range(len(i)):
        ddp_i=dp[:,i[k]]
        ddp_j=dp[:,j[k]]
        ddp=(ddp_i.sum(axis=1)>0) & (ddp_j.sum(axis=1)>0)
        duais_coup=np.arange(len(ddp))[np.array(ddp).T[0]]
        if len(duais_coup)>2:
            import pdb; pdb.set_trace()
        if len(duais_coup)==1:
            duais_coup=np.repeat(duais_coup[0],2)
        if len(duais_coup)==2:
            cds.append(duais_coup)

    cds=np.array(cds)
    if len(cds)>0:
        values=np.unique(np.concatenate(cds))
        mapd=np.arange(len(dual_volumes))
        mapd[values]=np.arange(len(values))

        lines=np.concatenate([mapd[cds[:,0]],mapd[cds[:,1]]])
        cols=np.concatenate([mapd[cds[:,1]],mapd[cds[:,0]]])

        data=np.ones(len(lines))
        graph=csc_matrix((data,(lines,cols)),shape=(len(values),len(values)))

        n_l,labels=csgraph.connected_components(graph,connection='strong')
        groups=[]
        for k in range(n_l):
            groups.append(values[labels==k])
        return groups
    else:
        return []

def get_dual_subdomains(groups):
    juntares=groups
    dv=[]
    for juntar in juntares:
        todos=np.arange(len(dual_volumes))
        keep_dual=np.setdiff1d(todos,juntar[1:])

        # dual_volumes = np.array(dual_volumes)
        dual_volumes2=dual_volumes[keep_dual]

        new_volume=np.unique(np.hstack(dual_volumes[juntar]))
        dv.append(new_volume)
    return(dv)

def group_dual_volumes_and_get_OP(mlo, T, M, data_impress, T_without_boundary, neta_lim):
    t0=time.time()

    OP_AMS=mlo['prolongation_level_1'].copy().tolil()
    groups = get_coupled_dual_volumes(mlo,T,M,data_impress,neta_lim, ind=0)

    dv=get_dual_subdomains(groups)
    if len(dv)>0:
        # mlo.run_paralel(tpfa_solver['Tini'],dv,1,False)
        mlo.run_paralel(T_without_boundary,dv,1,False)
        OP_AMS_groups=mlo['prolongation_level_1']
        lins_par=np.unique(np.concatenate(dv))
        OP_AMS[lins_par]=OP_AMS_groups[lins_par]
        mlo['prolongation_level_1']=OP_AMS
        multilevel_operators=mlo
        ###############################new_feature
        OP_AMS=mlo['prolongation_level_1']
        phiks=OP_AMS[data_impress['GID_0'],data_impress['GID_1']].toarray()[0]
        r_phiks=(1-phiks)/phiks
        maxs=csc_matrix((r_phiks,(data_impress['GID_0'],data_impress['GID_1'])),shape=OP_AMS.shape).max(axis=0).toarray()[0]
        M.maxs=maxs
        ########################end
    ######################### # #############################
    old_groups=groups.copy()
    if len(old_groups)==0:
        val=1
    else:
        val=6
    for ind in range(1,val):
        groups2 = get_coupled_dual_volumes(mlo,T,M,data_impress,neta_lim, ind=ind)
        # neta_lim/=2
        lgs=[np.repeat(i,len(old_groups[i])) for i in range(len(old_groups))]
        gs=np.concatenate(old_groups)
        lgs=np.concatenate(lgs)
        all_joined=np.zeros(len(gs))
        new_groups=[]
        for g2 in groups2:
            joins=np.zeros(len(gs))
            for g in g2:
                joins+= gs==g
                all_joined+= gs==g
            neighs=lgs[joins>0]
            if len(neighs)>0:
                aglomerated_dual=np.unique(np.concatenate([np.concatenate(np.array(old_groups)[neighs]),g2]))
            else:
                aglomerated_dual=g2
            new_groups.append(aglomerated_dual)
        inds_manteined=np.setdiff1d(lgs,np.unique(lgs[all_joined>0]))
        groups_manteined=np.array(old_groups)[inds_manteined]
        atualized_groups=groups_manteined.tolist()+new_groups

        print(ind,len(old_groups), len(atualized_groups),len(np.concatenate(atualized_groups)),"dsjjjjjja")
        dv=get_dual_subdomains(new_groups)
        if len(dv)>0:
            # mlo.run_paralel(tpfa_solver['Tini'],dv,1,False)
            mlo.run_paralel(T_without_boundary,dv,1,False)
            OP_AMS_groups=mlo['prolongation_level_1']
            lins_par=np.unique(np.concatenate(dv))
            OP_AMS[lins_par]=OP_AMS_groups[lins_par]
            mlo['prolongation_level_1']=OP_AMS
            multilevel_operators=mlo
        old_groups=atualized_groups.copy()
    mlo['prolongation_level_1']=OP_AMS
    print("Time to adapt RBC: {} seconds".format(time.time()-t0))
