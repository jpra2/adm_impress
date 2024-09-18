import numpy as np
import scipy.sparse as sp
from typing import Sequence
from packs.manager.boundary_conditions import BoundaryConditions


def define_NU_ADM_mesh(DUAL_1: np.ndarray, GID_0: np.ndarray, GID_1: np.ndarray, fs_vols: np.ndarray):
    levels = np.ones_like(GID_1)
    NU_ADM_ID = -levels
    levels[fs_vols] = 0
    coarse_volumes = levels == 1
    NU_ADM_ID[coarse_volumes]=GID_1[coarse_volumes]
    all_cvs=np.unique(NU_ADM_ID)
    if all_cvs.min()==-1:
        all_cvs=all_cvs[1:]
    remaining_ids=np.setdiff1d(np.unique(GID_1), all_cvs)
    nids=len(fs_vols)-len(remaining_ids)
    ids=np.concatenate([remaining_ids, GID_1.max()+np.arange(nids)+1])
    # import pdb; pdb.set_trace()
    NU_ADM_ID[fs_vols]=ids

    vertices=GID_0[DUAL_1==3]
    gid1=GID_1[vertices]
    for rgid in remaining_ids:
        if rgid in gid1:
            gid1[gid1==rgid]=NU_ADM_ID[vertices[GID_1[vertices]==rgid]]
    coarse_id_NU_ADM=gid1

    return coarse_id_NU_ADM, NU_ADM_ID


def update_averager_v0(GID_0, GID_1, fs_vols, adjs):
    # self.fs_vols=fs_vols
    levels=np.ones_like(GID_1)
    NU_ADM_ID = -levels
    levels[fs_vols]=0
    gid1_adjs=GID_1[adjs]
    same_gid=gid1_adjs[:,0]==gid1_adjs[:,1]
    cc_adjs=levels[adjs].sum(axis=1)==2
    adjs2=adjs[same_gid & cc_adjs]
    fines=np.tile(fs_vols,(2,1)).T
    adjs2=np.vstack([fines,adjs2])
    adjs2=np.tile(adjs2,(2,1))
    data = np.ones(len(adjs2))
    n=len(levels)
    graph = sp.csc_matrix((data, (adjs2[:,0], adjs2[:,1])),shape=(n,n))
    n,labels=sp.csgraph.connected_components(graph)
    '''
    self.NU_ADM_ID=labels
    gid1=self.GID_1[self.GID_0[self.DUAL_1==3]]
    self.coarse_id_NU_ADM=gid1
    '''
    cols=GID_0
    lines=labels
    data=np.ones_like(cols)
    # averager=sp.csc_matrix(())
    # import pdb; pdb.set_trace()
    averager=sp.csc_matrix((data, (lines, cols)), shape=(lines.max()+1, cols.max()+1))
    return averager



def update_NU_ADM_operators_v0(OP, levels, coarse_id_NU_ADM, GID_1, GID_0, NU_ADM_ID, fs_vols):
    
    
    l, c, d=OP
    coarse=levels[l]==1
    # import pdb; pdb.set_trace()
    # mapc = self.NU_ADM_ID[self.DUAL_1==3] #aqui trocar por linha abaixo
    mapc = coarse_id_NU_ADM
    lines = l[coarse]
    cols = mapc[c[coarse]]
    # import pdb; pdb.set_trace()
    same=GID_1[lines]==c[coarse]
    cols[same]=NU_ADM_ID[lines[same]]
    # import pdb; pdb.set_trace()
    # cols = self.NU_ADM_ID[lines]
    data = d[coarse]
    ls = fs_vols
    cs = NU_ADM_ID[fs_vols]
    ds = np.ones_like(cs)

    lines = np.concatenate([lines, ls])
    cols = np.concatenate([cols, cs])
    data = np.concatenate([data, ds])

    NU_ADM_OP = [lines, cols, data]
    # import pdb; pdb.set_trace()
    # visualize.plot_labels(self.OP[:,4].T.toarray()[0])
    # visualize.plot_labels(self.NU_ADM_OP[:,12].T.toarray()[0])
    # visualize.plot_labels(self.levels)
    # import pdb; pdb.set_trace()
    cols = GID_0
    lines = NU_ADM_ID
    data = np.ones(len(lines))
    NU_ADM_OR = [lines,cols, data]

    lp, cp, dp = NU_ADM_OP
    lr, cr, dr = NU_ADM_OR
    n_f, n_ADM=lp.max()+1, cp.max()+1

    OP_NU_ADM = sp.csc_matrix((dp, (lp, cp)), shape=(n_f, n_ADM))
    OR_NU_ADM = sp.csc_matrix((dr, (lr, cr)), shape=(n_ADM, n_f))

    return OP_NU_ADM, OR_NU_ADM, None

def get_beta_groups(GID_0, GID_1, OP, internal_adjacencies, beta_lim=3.0):
        adjs = internal_adjacencies
        pos=GID_1[OP[0]]==OP[1]
        index = OP[0][pos]
        index2 = np.setdiff1d(GID_0, index)
        # v1 = GID_1[OP[0]][pos]
        # v2 = OP[1][pos]
        phis=OP[2][pos][np.argsort(OP[0][pos])]
        # betas=(1-phis)/phis
        betas = np.zeros(pos.shape[0])
        betas[index] = (1-phis)/phis
        betas[index2] = np.inf
        beta_facs=betas[adjs].max(axis=1)
        ads=adjs[beta_facs>beta_lim]
        map=np.arange(adjs.max()+1)
        uads=np.unique(ads)
        n=len(uads)
        map[uads]=np.arange(n)
        adjs=map[ads]
        adjs=np.vstack([adjs,np.array([adjs[:,1],adjs[:,0]]).T])
        graph = sp.csc_matrix((np.ones_like(adjs[:,0]), (adjs[:,0], adjs[:,1])),shape=(n,n))
        n_l,labels=sp.csgraph.connected_components(graph)
        beta_ind=-np.ones_like(GID_0)
        beta_ind[uads]=labels
        beta_groups=np.array([uads[labels==l] for l in range(n_l)], dtype='O')
        return beta_groups, beta_ind, betas

def get_finescale_vols(inital_level0_vols:np.ndarray, alpha_vols_finescale:np.ndarray, beta_ind:np.ndarray, beta_groups:Sequence[np.ndarray]) -> np.ndarray:
    fs_vs = np.unique(np.concatenate([inital_level0_vols, alpha_vols_finescale]))
    bs=beta_ind[fs_vs]
    binds=np.unique(bs[bs>-1])

    if len(binds)>0:
        bvols=np.concatenate(beta_groups[binds])
        fs_vs=np.unique(np.concatenate([fs_vs,bvols]))
    
    return fs_vs

def set_adm_mesh_non_nested(
        v0:np.ndarray,
        levels,
        GID_0,
        GID_1,
        DUAL_1
    ):
    # levels = self.data_impress['LEVEL'].copy()
    gids_0 = GID_0
    gids_1 = GID_1
    n1 = 0
    n0 = len(levels)
    list_L1_ID = np.repeat(-1, n0)
    # list_L2_ID = np.repeat(-1, n0)
    list_L1_ID[v0] = np.arange(len(v0))
    # list_L2_ID[v0] = np.arange(len(v0))

    LEVEL_ID_1 = np.repeat(-1, n0)
    ADM_COARSE_ID_LEVEL_1 = np.repeat(-1, n0)

    # self.data_impress['LEVEL_ID_1'][v0] = list_L1_ID[v0]
    LEVEL_ID_1[v0] = list_L1_ID[v0]
    # self.data_impress['LEVEL_ID_2'][v0] = list_L2_ID[v0]
    # self.data_impress['ADM_COARSE_ID_LEVEL_2'][:] = -1
    # self.data_impress['ADM_COARSE_ID_LEVEL_1'][:] = -1

    n1+=len(v0)
    # n2+=len(v0)
    # ids_ms_2 = range(len(np.unique(gids_2)))

    ids_ms_1 = np.unique(GID_1)
    for vol1 in ids_ms_1:
        vols1 = gids_0[gids_1==vol1]
        levels_vols_1 = levels[vols1]
        vols_ms1_lv1 = vols1[levels_vols_1>=1]
        if len(vols_ms1_lv1)>0:
            list_L1_ID[vols_ms1_lv1] = np.repeat(n1,len(vols_ms1_lv1))
            ADM_COARSE_ID_LEVEL_1[vols1] = np.repeat(n1, len(vols1))
            n1+=1
        else:
            vertex = vols1[DUAL_1[vols1]==3]
            gid1_vertex = LEVEL_ID_1[vertex]
            ADM_COARSE_ID_LEVEL_1[vols1] = np.repeat(gid1_vertex, len(vols1))
    
    LEVEL_ID_1[:] = list_L1_ID
    n1_adm = n1
    return LEVEL_ID_1, ADM_COARSE_ID_LEVEL_1

def organize( 
        levels, 
        mm,
        GID_0,
        GID_1,
        LEVEL_ID_1,
        ADM_COARSE_ID_LEVEL_1,
        DUAL_1
    ):
    gid_0 = GID_0
    gid_level = GID_1
    adm_id = LEVEL_ID_1
    level_adm_coarse_id = ADM_COARSE_ID_LEVEL_1
    vertices = gid_0[DUAL_1==3]

    AMS_TO_ADM = np.arange(len(gid_level[vertices]))
    AMS_TO_ADM[gid_level[vertices]] = level_adm_coarse_id[vertices]

    gid_vols_nv0 = gid_0[levels==0]
    adm_vols_nv0 = adm_id[gid_vols_nv0]

    l1=mm[0]
    c1=mm[1]
    d1=mm[2].copy()
    lvs=levels[l1]
    d1[lvs==0]=0

    lines=gid_vols_nv0
    cols=adm_vols_nv0
    data=np.repeat(1,len(lines))
    c1_adm = AMS_TO_ADM[c1]

    lines = np.concatenate([lines,l1])
    cols = np.concatenate([cols,c1_adm])
    data = np.concatenate([data,d1])

    n1_adm = c1_adm.max()+1

    OP_ADM = sp.csc_matrix((data,(lines,cols)),shape=(len(gid_0),n1_adm))

    cols = gid_0
    lines = adm_id
    data = np.ones(len(lines))
    OR_ADM = sp.csc_matrix((data,(lines,cols)),shape=(n1_adm,len(gid_0)))

    return OP_ADM, OR_ADM

def calculate_fine_flux(
        GID_0: np.ndarray,
        GID_1: np.ndarray,
        pms_pressure: np.ndarray,
        fine_edges: np.ndarray,
        fine_edges_flux: np.ndarray,
        fine_boundary_edges: np.ndarray,
        fine_weights: np.ndarray,
        fine_transmissibility: sp.csc_matrix,
        DUAL_1: np.ndarray
):
    pass
    
    
