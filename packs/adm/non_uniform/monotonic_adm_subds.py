import numpy as np
import scipy.sparse as sp
import time
def get_preprossed_monotonic_primal_objects(data_impress,elements_lv0, OP_AMS,neumann_subds,phiK_raz_lim=5):
    OP_AMS_c=OP_AMS.tolil()
    gids1=data_impress['GID_1']
    ffv=np.vstack(elements_lv0['faces_face_volumes'][elements_lv0['internal_faces']])

    gids0_0=ffv[:,0]
    gids0_1=ffv[:,1]
    gids1_0=gids1[ffv[:,0]]
    gids1_1=gids1[ffv[:,1]]

    preprocessed_primal_objects=[]
    critical_groups=[]
    for primal in neumann_subds:
        gid1=primal.gid1
        intersect_faces=primal.intersect_faces
        adjs_intersect_faces=primal.adjs_intersect_faces
        adjs_intern_faces=primal.adj_intern_local_faces


        intern_volumes=np.concatenate(adjs_intern_faces)
        ids1_intersect=gids1[adjs_intersect_faces]
        vols_inter=adjs_intersect_faces[ids1_intersect==gid1]
        vols_exter=adjs_intersect_faces[ids1_intersect!=gid1]

        volumes=np.sort(np.unique(np.concatenate([intern_volumes, vols_inter, vols_exter])))
        map_l=np.zeros(volumes.max()+1,dtype=int)
        map_l[volumes]=range(len(volumes))
        local_vertex=map_l[primal.local_vertex]
        l=map_l[adjs_intern_faces[:,0]]
        c=map_l[adjs_intern_faces[:,1]]
        l_intern=np.concatenate([l, c, l, c])
        c_intern=np.concatenate([c, l, l, c])

        l_exter=np.concatenate([map_l[vols_exter],map_l[vols_exter]])
        c_exter=np.concatenate([map_l[vols_exter],map_l[vols_inter]])
        lc_exter=np.vstack([l_exter,c_exter])

        internal_faces=primal.intern_local_faces
        lines=np.concatenate([l_intern, map_l[vols_inter], map_l[vols_inter]]) #[l, c, l, c, inter, exter]
        cols=np.concatenate([c_intern, map_l[vols_inter], map_l[vols_exter]]) #[c, l, l, c, inter, exter]
        intern_global_faces=primal.intern_local_faces

        l0=gids0_0[(gids1_0==gid1) | (gids1_1==gid1)]
        l1=gids0_1[(gids1_0==gid1) | (gids1_1==gid1)]
        ls=np.unique(np.concatenate([l0,l1]))
        OP_AMS_local=OP_AMS_c[ls]
        gids1_exter=gids1[np.concatenate([vols_exter,vols_exter])]

        phi_k_exter=OP_AMS_local[l_exter,gids1_exter].toarray()[0]

        lcd_OP_local=sp.find(OP_AMS_local)
        prep_primal=PrepMonotonicPrimal(lines, cols, intern_global_faces, intersect_faces, lcd_OP_local, gid1, volumes, lc_exter, phi_k_exter, phiK_raz_lim, data_impress, local_vertex)
        preprocessed_primal_objects.append(prep_primal)
        for cg in prep_primal.critical_groups:
            critical_groups.append(np.unique(cg))
    return preprocessed_primal_objects, critical_groups

class PrepMonotonicPrimal:
    def __init__(self, lines, cols, intern_global_faces, intersect_faces, lcd_OP_local, gid1, volumes, lc_exter, phi_k_exter, phiK_raz_lim, data_impress, local_vertex):
        self.create_primal_subds(lines, cols, intern_global_faces, intersect_faces, lcd_OP_local, gid1, volumes, lc_exter, phi_k_exter,phiK_raz_lim, data_impress, local_vertex)

    def create_primal_subds(self, lines, cols, intern_global_faces, intersect_faces, lcd_OP_local, gid1, volumes, lc_exter, phi_k_exter,phiK_raz_lim, data_impress, local_vertex):
        # import pdb; pdb.set_trace()
        self.lc_exter=lc_exter
        self.equals=lines==cols
        self.different=lines!=cols
        self.le=lines[self.equals]
        self.ce=cols[self.equals]
        self.ld=lines[self.different]
        self.cd=cols[self.different]
        self.local_vertex=local_vertex

        self.ls=np.unique(lines[self.equals])
        self.lines=lines
        self.cols=cols
        self.intern_global_faces=intern_global_faces
        self.intersect_faces=intersect_faces
        lcd=lcd_OP_local
        cols=lcd[1]
        map_c=np.zeros(cols.max()+1,dtype=int)
        map_c[np.unique(cols)]=range(len(np.unique(cols)))

        self.OP_local=sp.csc_matrix((lcd[2],(lcd[0], map_c[lcd[1]])),shape=(lcd[0].max()+1,lcd[1].max()+1))
        self.OP_numpy=self.OP_local.toarray()
        self.T_numpy=np.zeros((len(volumes), len(volumes)))
        self.gid1_local=map_c[gid1]
        self.volumes=volumes
        # import pdb; pdb.set_trace()
        # phi_k=self.OP_local[:,self.gid1_local]
        phi_k=self.OP_local[:,self.gid1_local].T.toarray()[0]
        # phi_k[np.setdiff1d(np.arange(len(volumes)),np.unique(self.lines))]=1
        # import pdb; pdb.set_trace()
        phi_k[lc_exter[0]]=phi_k_exter
        raz_phi=(1-phi_k)/phi_k

        actual_phi=data_impress['raz_phi'][volumes]
        efective_phi=np.maximum(actual_phi,raz_phi)
        data_impress['raz_phi'][volumes]=efective_phi

        phi_subs=1
        critical_groups=[]
        lc_exter=self.lc_exter

        ls=np.concatenate([self.lines,lc_exter[0]])
        cs=np.concatenate([self.cols,lc_exter[1]])
        # import pdb; pdb.set_trace()
        if raz_phi.max()>phiK_raz_lim:
            critical_volumes=np.arange(len(volumes))[raz_phi>phiK_raz_lim]
            map_g=np.repeat(-1,len(volumes))
            map_g[critical_volumes]=range(len(critical_volumes))
            pos_crit=(map_g[ls]>-1) & (map_g[cs]>-1)
            lg=map_g[np.arange(len(volumes))[ls[pos_crit]]]
            cg=map_g[np.arange(len(volumes))[cs[pos_crit]]]
            dg=np.ones(len(lg))
            graph=sp.csc_matrix((dg, (lg, cg)), shape=(len(critical_volumes), len(critical_volumes)))
            nc, labels = sp.csgraph.connected_components(csgraph=graph, directed=False, return_labels=True)
            critical_groups=[]
            for i in range(nc):
                pos=labels==i
                if pos.sum()>1:
                    critical_groups.append(volumes[critical_volumes[pos]])
        self.critical_groups=critical_groups

def get_monotonizing_volumes(preprocessed_primal_objects, transmissibility):
    volumes=[]
    netasp_list=[]
    for primal in preprocessed_primal_objects:
        t_internal=transmissibility[primal.intern_global_faces]
        t_intersect=transmissibility[primal.intersect_faces]
        lines=primal.lines
        cols=primal.cols
        data=np.concatenate([t_internal, t_internal, -t_internal, -t_internal, -t_intersect, t_intersect])
        T_numpy=primal.T_numpy
        equals=primal.equals
        different=primal.different
        le=primal.le
        ce=primal.ce
        ld=primal.ld
        cd=primal.cd

        ls=primal.ls
        T_numpy[ld,cd]=data[different]

        T_numpy[ls, ls]=np.bincount(le,weights=data[equals])[ls]

        ###############Just for the example
        if len(T_numpy)==8:
            T_numpy[primal.volumes==15]=0
            T_numpy[primal.volumes==15,primal.volumes==15]=1
        ############################

        OP_numpy=primal.OP_numpy
        gid1=primal.gid1_local
        tpn=np.dot(T_numpy,OP_numpy)


        dcn=tpn[:,gid1].sum()
        nfn=(1/dcn)*tpn

        nfn[:,gid1]=0 #

        # netaspn=nfn.max(axis=1)
        # lines_out_prim=np.arange(len(T_numpy))[T_numpy.max(axis=1)==0]
        # nfn[lines_out_prim]=0

        netaspn=nfn.max(axis=1)

        netaspn[primal.local_vertex]=0

        volumes.append(primal.volumes)
        netasp_list.append(netaspn)
    netasp_array=np.hstack(netasp_list)
    volumes=np.concatenate(volumes)
    return volumes, netasp_array

def incorporate_fs_to_groups(l_groups, groups_c, critical_groups, data_impress, elements_lv0):
    print('start to incorporate')
    t0=time.time()
    v0s=data_impress['GID_0'][data_impress['LEVEL']==0]
    data_impress['nfp'][:]=-1
    data_impress['nfp'][groups_c]=l_groups.astype('int')
    viz=elements_lv0['volumes_face_volumes'][v0s]
    for i in range(len(v0s)):
        for j in data_impress['nfp'][viz[i]]:
            if j>=0:
                try:
                    critical_groups[int(j)]=np.append(critical_groups[int(j)],v0s[i])
                    l_groups=np.append(l_groups,int(j))
                    groups_c=np.append(groups_c,v0s[i])
                except:
                    pass

    print('finish incorporate', time.time()-t0)
    return l_groups, groups_c, critical_groups


def get_monotonizing_level(l_groups, groups_c, critical_groups,data_impress,elements_lv0, volumes,netasp_array, tol):
    maxs=np.zeros(len(np.unique(volumes)))
    np.maximum.at(maxs,volumes,netasp_array)
    data_impress['ref_param'][np.unique(volumes)]=maxs
    vols_orig=np.unique(volumes)[maxs>tol]

    data_impress['LEVEL'][vols_orig]=0

    # l_groups, groups_c, critical_groups=incorporate_fs_to_groups(l_groups, groups_c, critical_groups, data_impress, elements_lv0)
    t1=time.time()
    teste=True
    vols_0=np.array([])
    while teste:
        lv=len(vols_0)
        groups_lv0=np.unique(l_groups[data_impress['LEVEL'][groups_c]==0])
        if len(groups_lv0)>0 and len(critical_groups) > 0:
            vols_lv0=np.concatenate(np.array(critical_groups)[groups_lv0])
            vols_0=np.unique(np.append(vols_0,vols_lv0))
        if len(vols_0)==lv:
            teste=False
        else:
            vols_lv0=np.concatenate(np.array(critical_groups)[groups_lv0])
            data_impress['LEVEL'][vols_lv0]=0
    vols_orig=np.unique(np.concatenate([vols_orig,vols_0])).astype(int)
    data_impress['LEVEL'][vols_orig]=0
    return vols_orig
