from ..adm_method import AdmMethod, Jacobi
import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg
import time
from ...directories import file_adm_mesh_def
import matplotlib.pyplot as plt
from pymoab import types
from scipy.sparse import csc_matrix, csgraph
# from .paralel_neuman import masterNeumanNonNested
# from .paralel_neuman_new1 import masterNeumanNonNested
from .paralell_neumann_numpy import masterNeumanNonNested

class AdmNonNested(AdmMethod):

    def set_adm_mesh_non_nested(self, v0=[], v1=[], pare=False):
        levels = self.data_impress['LEVEL'].copy()
        gids_0 = self.data_impress['GID_0']
        gids_1 = self.data_impress['GID_1']
        gids_2 = self.data_impress['GID_2']
        n1 = 0
        n2 = 0
        n0 = len(levels)
        list_L1_ID = np.repeat(-1, n0)
        list_L2_ID = np.repeat(-1, n0)
        list_L1_ID[v0] = np.arange(len(v0))
        list_L2_ID[v0] = np.arange(len(v0))

        self.data_impress['LEVEL_ID_1'][v0] = list_L1_ID[v0]
        self.data_impress['LEVEL_ID_2'][v0] = list_L2_ID[v0]
        self.data_impress['ADM_COARSE_ID_LEVEL_2'][:] = -1
        self.data_impress['ADM_COARSE_ID_LEVEL_1'][:] = -1

        n1+=len(v0)
        n2+=len(v0)
        ids_ms_2 = range(len(np.unique(gids_2)))

        for vol2 in ids_ms_2:
            vols2 = gids_0[gids_2==vol2]
            levels_vols_2 = levels[vols2]
            vols_ms2_lv2 = vols2[levels_vols_2==2]
            list_L2_ID[vols_ms2_lv2] = np.repeat(n2,len(vols_ms2_lv2))
            self.data_impress['ADM_COARSE_ID_LEVEL_2'][vols2] = np.repeat(n2, len(vols2))
            if len(vols_ms2_lv2)>0:
                n2+=1
            gids_1_1 = gids_1[vols2]
            ids_ms_1 = np.unique(gids_1_1)
            for vol1 in ids_ms_1:
                vols1 = vols2[gids_1_1==vol1]
                levels_vols_1 = levels_vols_2[gids_1_1==vol1]
                vols_ms1_lv1 = vols1[levels_vols_1>=1]

                if len(vols_ms1_lv1)>0:
                    list_L1_ID[vols_ms1_lv1] = np.repeat(n1,len(vols_ms1_lv1))
                    self.data_impress['ADM_COARSE_ID_LEVEL_1'][vols1] = np.repeat(n1, len(vols1))
                    n1+=1
                else:
                    vertex = vols1[self.data_impress['DUAL_1'][vols1]==3]
                    gid1_vertex = self.data_impress['LEVEL_ID_1'][vertex]
                    self.data_impress['ADM_COARSE_ID_LEVEL_1'][vols1] = np.repeat(gid1_vertex, len(vols1))

                vols_ms2_lv1 = vols1[levels_vols_1==1]
                if len(vols_ms2_lv1)>0:
                    list_L2_ID[vols_ms2_lv1] = n2
                    n2+=1

        self.data_impress['LEVEL_ID_1'] = list_L1_ID
        self.data_impress['LEVEL_ID_2'] = list_L2_ID

        for i in range(self.n_levels):
            self.number_vols_in_levels[i] = len(levels[levels==i])
        self.n1_adm = n1
        self.n2_adm = n2

    def equalize_levels(self):

        levels = self.data_impress['LEVEL'].copy()
        gid0 = self.data_impress['GID_0']
        # ############
        # coarse_id_wells = np.unique(self.data_impress['GID_1'][self.all_wells_ids])
        # for cid in coarse_id_wells:
        #     levels[self.data_impress['GID_1']==cid] = 0
        # ############
        gids_with_level = gid0[levels >= 0]

        for level in range(1, self.n_levels):
            cids_level = self.data_impress['GID_'+str(level)]
            all_cids = np.unique(cids_level[levels == level])

            for cid in all_cids:
                gids_coarse = gid0[cids_level==cid]
                n0 = gids_coarse.shape[0]
                gids_coarse_without_level = np.setdiff1d(gids_coarse, gids_with_level)
                n1 = gids_coarse_without_level.shape[0]
                if n1 == n0:
                    continue
                elif n1 == 0:
                    continue
                else:
                    # n1 < n0:
                    levels_gids_coarse = levels[gids_coarse]
                    gids_coarse_with_level = np.setdiff1d(gids_coarse, gids_coarse_without_level)
                    levels_gids_coarse_with_level = levels[gids_coarse_with_level]
                    for llevel in np.unique(levels_gids_coarse_with_level):
                        if llevel == self.n_levels-1:
                            continue
                        vols = gids_coarse_with_level[levels_gids_coarse_with_level==llevel]
                        neig_vols = np.unique(np.concatenate(self.elements_lv0['volumes_face_volumes'][vols]))
                        levels_neig_vols = levels[neig_vols]
                        vols_without_level = neig_vols[levels_neig_vols < 0]
                        n3 = vols_without_level.shape[0]
                        if n3 > 0:
                            levels[vols_without_level] = llevel + 1

        levels[levels < 0] = self.n_levels-1
        self.data_impress['LEVEL'] = levels

    def set_saturation_level(self):

        levels = self.data_impress['LEVEL'].copy()
        import pdb; pdb.set_trace()
        gid1 = self.data_impress['GID_1']
        gid0 = self.data_impress['GID_0']
        level_0_ini = set(gid0[levels==0])
        saturation = self.data_impress['saturation']
        all_wells = set(self.all_wells_ids)
        gids_lv1_sat = set()
        gidsc = np.unique(gid1)
        internal_faces = self.elements_lv0['internal_faces']
        v0 = self.elements_lv0['neig_internal_faces']

        ds = saturation[v0]
        ds = np.absolute(ds[:,1] - ds[:,0])

        inds = ds >= self.delta_sat_max

        levels[v0[inds][:,0]] = 0
        levels[v0[inds][:,1]] = 0

        all_lv0 = set(gid0[levels==0])

        for gidc in gidsc:
            gids0 = gid0[gid1==gidc]
            if set(gids0) & all_lv0:
                gids_fora = np.array(list(set(gids0) - all_lv0))
                if len(gids_fora) > 0:
                    levels[gids_fora] = 1
                gids_lv1_sat.add(gidc)

        cids_neigh = self.ml_data['coarse_id_neig_face_level_'+str(1)]
        cids_level = self.ml_data['coarse_primal_id_level_'+str(1)]

        for gidc in gids_lv1_sat:
            vizs = cids_neigh[cids_level==gidc][0]
            for viz in vizs:
                if set([viz]) & gids_lv1_sat:
                    continue
                gids0 = gid0[gid1==viz]
                if set(gids0) & all_lv0:
                    gids_fora = np.array(list(set(gids0) - all_lv0))
                    levels[gids_fora] = 1
                else:
                    levels[gids0] = 1

        self.data_impress['LEVEL'] = levels.copy()

    def set_saturation_level_simple(self,delta_sat_max=0.1):
        levels = self.data_impress['LEVEL'].copy()
        saturation = self.data_impress['saturation']
        internal_faces = self.elements_lv0['internal_faces']
        v0 = self.elements_lv0['neig_internal_faces']
        gids1=self.data_impress['GID_1'][v0]
        ds = saturation[v0]
        # ds = ds.sum(axis=1)
        ds = np.absolute(ds[:,1] - ds[:,0])

        inds = (ds >= delta_sat_max) & (saturation[v0].min(axis=1)<min(saturation[v0].min(),0)+0.5)
        # inds = ds >= delta_sat_max
        levels[v0[inds][:,0]] = 0
        levels[v0[inds][:,1]] = 0
        self.data_impress['LEVEL'] = levels.copy()

    def set_saturation_level_homogeneo(self,delta_sat_max=0.1):
        levels = self.data_impress['LEVEL'].copy()
        saturation = self.data_impress['saturation']
        internal_faces = self.elements_lv0['internal_faces']
        v0 = self.elements_lv0['neig_internal_faces']
        gids1=self.data_impress['GID_1'][v0]
        gids_coarse=self.data_impress['GID_1']
        ds = saturation[v0]
        # ds = ds.sum(axis=1)
        ds = np.absolute(ds[:,1] - ds[:,0])


        inds = (ds >= delta_sat_max) & (saturation[v0].min(axis=1)<min(saturation[v0].min(),0)+0.5)
        # inds = ds >= delta_sat_max

        levels[v0[inds][:,0]] = 0
        levels[v0[inds][:,1]] = 0

        gids_fin=np.unique(gids_coarse[levels==0])
        for gid in gids_fin:
            levels[gids_coarse==gid]=0

        self.data_impress['LEVEL'] = levels.copy()

    def set_saturation_level_uniform(self, delta_max=0.1):
        levels = self.data_impress['LEVEL'].copy()
        gids1 = self.data_impress['GID_1']
        saturation = self.data_impress['saturation']
        mins=np.ones(len(np.unique(gids1)))
        maxs=np.zeros(len(np.unique(gids1)))
        np.maximum.at(maxs,gids1,saturation)
        np.minimum.at(mins,gids1,saturation)
        deltas=maxs-mins
        levels[deltas[gids1]>delta_max]=0
        self.data_impress['LEVEL'] = levels.copy()

    def set_saturation_level_imposed_joined_coarse(self):

        levels = self.data_impress['LEVEL'].copy()
        dual_flag = self.data_impress['DUAL_1'].copy()
        gid1 = self.data_impress['GID_1']
        gid0 = self.data_impress['GID_0']
        level_0_ini = set(gid0[levels==0])
        saturation = self.data_impress['saturation']
        all_wells = set(self.all_wells_ids)
        gids_lv1_sat = set()
        gidsc = np.unique(gid1)
        internal_faces = self.elements_lv0['internal_faces']
        v0 = self.elements_lv0['neig_internal_faces']

        ds = saturation[v0]
        ds = np.absolute(ds[:,1] - ds[:,0])

        inds = ds >= self.delta_sat_max

        levels[v0[inds][:,0]] = 0
        levels[v0[inds][:,1]] = 0

        all_lv0 = set(gid0[levels==0])
        for gidc in gidsc:
            gids0 = gid0[gid1==gidc]
            vertex = gids0[dual_flag[gids0]==3]
            if (levels[gids0].max()-levels[gids0].min())>0:
                facs=np.unique(np.concatenate(self.elements_lv0['volumes_face_faces'][gids0]))
                facs=np.intersect1d(facs,internal_faces)
                ad=np.vstack(self.elements_lv0['faces_face_volumes'][facs])
                ad0=ad[:,0]
                ad1=ad[:,1]
                l0=levels[ad0]
                l1=levels[ad1]
                map_lid=-np.ones(max(gid0)+1)
                map_lid[gids0]=np.arange(len(gids0))
                l0[map_lid[ad0]<0]=-1
                l1[map_lid[ad1]<0]=-1
                fadj1=l0+l1>=0
                # import pdb; pdb.set_trace()
                lines=map_lid[ad0[fadj1]].astype(int)
                cols=map_lid[ad1[fadj1]].astype(int)
                data=np.ones(len(lines))
                graph=csc_matrix((data,(lines,cols)),shape=(len(gids0),len(gids0)))
                n_l,labels=csgraph.connected_components(graph,connection='weak')
                groups=[gids0[labels==k] for k in range(n_l)]
                ls=np.array([len(g) for g in groups])
                print(ls.max())
                if ls.max()>1:
                    vols_nv1=np.array(groups)[ls==ls.max()][0]
                    levels[np.setdiff1d(gids0,vols_nv1)]=0

                # import pdb; pdb.set_trace()


        for gidc in gidsc:
            gids0 = gid0[gid1==gidc]
            if set(gids0) & all_lv0:
                gids_fora = np.array(list(set(gids0) - all_lv0))
                if len(gids_fora) > 0:
                    levels[gids_fora] = 1
                gids_lv1_sat.add(gidc)

        cids_neigh = self.ml_data['coarse_id_neig_face_level_'+str(1)]
        cids_level = self.ml_data['coarse_primal_id_level_'+str(1)]

        for gidc in gids_lv1_sat:
            vizs = cids_neigh[cids_level==gidc][0]
            for viz in vizs:
                if set([viz]) & gids_lv1_sat:
                    continue
                gids0 = gid0[gid1==viz]
                if set(gids0) & all_lv0:
                    gids_fora = np.array(list(set(gids0) - all_lv0))
                    levels[gids_fora] = 1
                else:
                    levels[gids0] = 1

        self.data_impress['LEVEL'] = levels.copy()

    def set_saturation_level_imposed_bound_level_continuity(self):

        levels = self.data_impress['LEVEL'].copy()
        dual_flag = self.data_impress['DUAL_1'].copy()
        gid1 = self.data_impress['GID_1']
        gid0 = self.data_impress['GID_0']
        level_0_ini = set(gid0[levels==0])
        saturation = self.data_impress['saturation']
        all_wells = set(self.all_wells_ids)
        gids_lv1_sat = set()
        gidsc = np.unique(gid1)
        internal_faces = self.elements_lv0['internal_faces']
        v0 = self.elements_lv0['neig_internal_faces']
        ds = saturation[v0]
        ds = np.absolute(ds[:,1] - ds[:,0])
        inds = ds >= self.delta_sat_max
        levels[v0[inds][:,0]] = 0
        levels[v0[inds][:,1]] = 0
        c0=gid1[v0[:,0]]
        c1=gid1[v0[:,1]]
        fb=c0!=c1
        c0b=v0[:,0][fb]
        c1b=v0[:,1][fb]
        levels[c1b[levels[c0b]==0]]=0
        levels[c0b[levels[c1b]==0]]=0

        c0b=v0[:,0][fb]
        c1b=v0[:,1][fb]
        levels[c1b[levels[c0b]==0]]=0
        levels[c0b[levels[c1b]==0]]=0

        all_lv0 = set(gid0[levels==0])
        # for gidc in gidsc:
        #     gids0 = gid0[gid1==gidc]




        # for gidc in gidsc:
        #     gids0 = gid0[gid1==gidc]
        #     if set(gids0) & all_lv0:
        #         gids_fora = np.array(list(set(gids0) - all_lv0))
        #         if len(gids_fora) > 0:
        #             levels[gids_fora] = 1
        #         gids_lv1_sat.add(gidc)

        cids_neigh = self.ml_data['coarse_id_neig_face_level_'+str(1)]
        cids_level = self.ml_data['coarse_primal_id_level_'+str(1)]

        for gidc in gids_lv1_sat:
            vizs = cids_neigh[cids_level==gidc][0]
            for viz in vizs:
                if set([viz]) & gids_lv1_sat:
                    continue
                gids0 = gid0[gid1==viz]
                if set(gids0) & all_lv0:
                    gids_fora = np.array(list(set(gids0) - all_lv0))
                    levels[gids_fora] = 1
                else:
                    levels[gids0] = 1

        self.data_impress['LEVEL'] = levels.copy()

    def set_saturation_level_original(self):
        levels = self.data_impress['LEVEL'].copy()
        dual_flag = self.data_impress['DUAL_1'].copy()
        gid1 = self.data_impress['GID_1']
        gid0 = self.data_impress['GID_0']
        level_0_ini = set(gid0[levels==0])
        saturation = self.data_impress['saturation']
        all_wells = set(self.all_wells_ids)
        gids_lv1_sat = set()
        gidsc = np.unique(gid1)
        internal_faces = self.elements_lv0['internal_faces']
        v0 = self.elements_lv0['neig_internal_faces']
        # ds = saturation[v0]
        # ds = np.absolute(ds[:,1] - ds[:,0])
        # inds = ds >= self.delta_sat_max
        # levels[v0[inds][:,0]] = 0
        # levels[v0[inds][:,1]] = 0
        # c0=gid1[v0[:,0]]
        # c1=gid1[v0[:,1]]
        # fb=c0!=c1
        # c0b=v0[:,0][fb]
        # c1b=v0[:,1][fb]
        # levels[c1b[levels[c0b]==0]]=0
        # levels[c0b[levels[c1b]==0]]=0
        #
        # c0b=v0[:,0][fb]
        # c1b=v0[:,1][fb]
        # levels[c1b[levels[c0b]==0]]=0
        # levels[c0b[levels[c1b]==0]]=0

        all_lv0 = set(gid0[levels==0])
        for gidc in gidsc:
            gids0 = gid0[gid1==gidc]
            ds = saturation[gids0].max()-saturation[gids0].min()
            if ds >= self.delta_sat_max:
                print("aqui",gidc, gids0)
                levels[gids0[dual_flag[gids0]>1]]=0

        # for gidc in gidsc:
        #     gids0 = gid0[gid1==gidc]
        #     if set(gids0) & all_lv0:
        #         gids_fora = np.array(list(set(gids0) - all_lv0))
        #         if len(gids_fora) > 0:
        #             levels[gids_fora] = 1
        #         gids_lv1_sat.add(gidc)

        cids_neigh = self.ml_data['coarse_id_neig_face_level_'+str(1)]
        cids_level = self.ml_data['coarse_primal_id_level_'+str(1)]

        for gidc in gids_lv1_sat:
            vizs = cids_neigh[cids_level==gidc][0]
            for viz in vizs:
                if set([viz]) & gids_lv1_sat:
                    continue
                gids0 = gid0[gid1==viz]
                if set(gids0) & all_lv0:
                    gids_fora = np.array(list(set(gids0) - all_lv0))
                    levels[gids_fora] = 1
                else:
                    levels[gids0] = 1

        self.data_impress['LEVEL'] = levels.copy()

    def set_saturation_level_new0(self):

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
            sats_local = saturation[gids0]
            dif = sats_local.max() - sats_local.min()
            if dif >= self.delta_sat_max:
                levels[gids0] = 0
                gids_lv1_sat.add(gidc)

        cids_neigh = self.ml_data['coarse_id_neig_face_level_'+str(1)]
        cids_level = self.ml_data['coarse_primal_id_level_'+str(1)]

        all_lv0 = set(gid0[levels==0])

        for gidc in gids_lv1_sat:
            vizs = cids_neigh[cids_level==gidc][0]
            for viz in vizs:
                if set([viz]) & gids_lv1_sat:
                    continue
                gids0 = gid0[gid1==viz]
                if set(gids0) & all_lv0:
                    continue
                else:
                    levels[gids0] = 1

        self.data_impress['LEVEL'] = levels.copy()

    def set_level_wells_3(self):
        self.data_impress['LEVEL'][:] = 1
        gid1_wells=self.data_impress['GID_1'][self.all_wells_ids]
        pos=np.zeros_like(self.data_impress['LEVEL'])
        for g1 in gid1_wells:
            pos+=self.data_impress['GID_1']==g1
        vols_to_lv0=self.data_impress['GID_0'][pos>0]

        # vols_to_lv0=np.unique(np.hstack(self.elements_lv0['volumes_face_volumes'][vols_to_lv0]))
        self.data_impress['LEVEL'][self.all_wells_ids] = np.zeros(len(self.all_wells_ids))
        self.data_impress['LEVEL'][vols_to_lv0] = 0

    def set_level_wells_only(self):

        self.data_impress['LEVEL'][:] = 1
        self.data_impress['LEVEL'][self.all_wells_ids] = 0

    def set_monotonizing_level(self, gids_to_monotonize):
        self.data_impress['LEVEL'][gids_to_monotonize]=0

    def organize_ops_adm(self, mlo, level, _pcorr=None):

        OP_AMS_lcd=mlo['prolongation_lcd_level_1']
        # OP_AMS=mlo['prolongation_level_1'].tolil()
        OR_AMS=mlo['restriction_level_1']
        t0=time.time()
        gid_0 = self.data_impress['GID_0']
        gid_level = self.data_impress['GID_' + str(level)]
        gid_ant = self.data_impress['GID_' + str(level-1)]
        level_id = self.data_impress['LEVEL_ID_' + str(level)]
        level_id_ant = self.data_impress['LEVEL_ID_' + str(level-1)]
        levels = self.data_impress['LEVEL']
        dual_id = self.data_impress['DUAL_'+str(level)]
        # OP_AMS = OP_AMS.tolil()

        if level == 1:
            lcd=OP_AMS_lcd
            # t0=time.time()
            # OP_ADM1, OR_ADM1, pcorr1 = self.organize_ops_adm_level_1(OP_AMS, OR_AMS, level, _pcorr=_pcorr)
            OP_ADM, OR_ADM, pcorr =self.organize(level,lcd, _pcorr=_pcorr)

            self._data[self.adm_op_n + str(level)] = OP_ADM
            self._data[self.adm_rest_n + str(level)] = OR_ADM
            self._data[self.pcorr_n+str(level-1)] = pcorr
            return 0

        import pdb; pdb.set_trace()

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


        self._data[self.adm_op_n + str(level)] = OP_ADM
        self._data[self.adm_rest_n + str(level)] = OR_ADM

        return 0

    def organize_ops_adm_level_1(self, OP_AMS, OR_AMS, level, _pcorr=None):

        gid_0 = self.data_impress['GID_0']
        gid_level = self.data_impress['GID_' + str(level)]
        gid_ant = self.data_impress['GID_' + str(level-1)]
        level_id = self.data_impress['LEVEL_ID_' + str(level)]
        level_adm_coarse_id = self.data_impress['ADM_COARSE_ID_LEVEL_1']
        level_id_ant = self.data_impress['LEVEL_ID_' + str(level-1)]
        levels = self.data_impress['LEVEL']
        vertices = gid_0[self.data_impress['DUAL_1']==3]

        OP_AMS_local = OP_AMS.copy().tolil()

        # import pdb; pdb.set_trace()
        AMS_TO_ADM = np.arange(len(gid_level[vertices]))
        AMS_TO_ADM[gid_level[vertices]] = level_adm_coarse_id[vertices]

        nivel_0 = gid_0[levels==0]
        ID_global1 = nivel_0

        OP_AMS_local[nivel_0] = 0

        n1_adm = len(np.unique(level_id))

        ids_adm_nivel0 = level_id[nivel_0]
        IDs_ADM1 = ids_adm_nivel0

        m = sp.find(OP_AMS_local)
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
        try:
            OP_ADM = sp.csc_matrix((data,(lines,cols)),shape=(len(gid_0),n1_adm))

            # plot_operator(self, OP_ADM, np.arange())
            # import pdb; pdb.set_trace()
        except:
            import pdb; pdb.set_trace()

        cols = gid_0
        lines = level_id
        data = np.ones(len(lines))
        OR_ADM = sp.csc_matrix((data,(lines,cols)),shape=(n1_adm,len(gid_0)))

        if self.get_correction_term:
            pcorr = np.zeros(len(gid_0))
            pcorr[levels>0] = _pcorr[levels>0]
        else:#np.arange(729)[np.array(OP_ADM.sum(axis=1)==0).T[0]]
            pcorr = np.array([False])
        if OP_ADM.sum()<OP_ADM.shape[0]-0.1:
            print("verify ADM prolongation operator organization")
            import pdb; pdb.set_trace()
        return OP_ADM, OR_ADM, pcorr

    def organize(self, level, mm,_pcorr=None):
        gid_0 = self.data_impress['GID_0']
        gid_level = self.data_impress['GID_' + str(level)]
        adm_id = self.data_impress['LEVEL_ID_' + str(level)]
        level_adm_coarse_id = self.data_impress['ADM_COARSE_ID_LEVEL_1']
        levels = self.data_impress['LEVEL']
        vertices = gid_0[self.data_impress['DUAL_1']==3]

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

        # AMS_TO_ADM = np.arange(len(gid_level[vertices]))
        # AMS_TO_ADM[gid_level[vertices]] = level_adm_coarse_id[vertices]

        n1_adm = c1_adm.max()+1

        OP_ADM = sp.csc_matrix((data,(lines,cols)),shape=(len(gid_0),n1_adm))

        cols = gid_0
        lines = adm_id
        data = np.ones(len(lines))
        OR_ADM = sp.csc_matrix((data,(lines,cols)),shape=(n1_adm,len(gid_0)))

        if self.get_correction_term:
            pcorr = np.zeros(len(gid_0))
            pcorr[levels>0] = _pcorr[levels>0]
        else:
            pcorr = np.array([False])
        if OP_ADM.sum()<OP_ADM.shape[0]-0.1:
            print("verify ADM prolongation operator organization")
            import pdb; pdb.set_trace()
        return OP_ADM, OR_ADM, pcorr

    def organize_ops_adm_level_1_dep(self, OP_AMS, OR_AMS, level, _pcorr=None):
        gid_0 = self.data_impress['GID_0']
        gid_level = self.data_impress['GID_' + str(level)]
        gid_ant = self.data_impress['GID_' + str(level-1)]
        level_id = self.data_impress['LEVEL_ID_' + str(level)]
        level_adm_coarse_id = self.data_impress['ADM_COARSE_ID_LEVEL_1']
        level_id_ant = self.data_impress['LEVEL_ID_' + str(level-1)]
        levels = self.data_impress['LEVEL']
        vertices = gid_0[self.data_impress['DUAL_1']==3]
        OP_AMS_local = OP_AMS.copy().tolil()

        AMS_TO_ADM = np.arange(len(gid_level[vertices]))
        AMS_TO_ADM[gid_level[vertices]] = level_adm_coarse_id[vertices]

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
        try:
            OP_ADM = sp.csc_matrix((data,(lines,cols)),shape=(len(gid_0),n1_adm))
            # plot_operator(self, OP_ADM, np.arange())
            # import pdb; pdb.set_trace()
        except:
            import pdb; pdb.set_trace()

        cols = gid_0
        lines = level_id
        data = np.ones(len(lines))
        OR_ADM = sp.csc_matrix((data,(lines,cols)),shape=(n1_adm,len(gid_0)))

        if self.get_correction_term:
            pcorr = np.zeros(len(gid_0))
            pcorr[levels>0] = _pcorr[levels>0]
        else:#np.arange(729)[np.array(OP_ADM.sum(axis=1)==0).T[0]]
            pcorr = np.array([False])

        return OP_ADM, OR_ADM, pcorr

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

        OP_AMS_1 = mlo['prolongation_level_1']
        OR_AMS_1 = mlo['restriction_level_1']

        if load_adm_levels:
            return 0

        if not set_initial_mesh:
            self.restart_levels()
            self.set_level_wells()
            gid_0 = self.data_impress['GID_0'][self.data_impress['LEVEL']==0]
            self.set_adm_mesh_non_nested(gid_0)
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
        gid_0 = self.data_impress['GID_0'][self.data_impress['LEVEL']==0]
        gid_1 = self.data_impress['GID_0'][self.data_impress['LEVEL']==1]
        self.set_adm_mesh_non_nested(gid_0, gid_1)

        multilevel_meshes = []

        active_nodes = []
        perro = []
        erro = []

        Nmax = tol_n2*nfine_vols
        finos = self.all_wells_ids.copy()

        vertices = GID_0[DUAL_1==3]
        primal_id_vertices = GID_1[vertices]
        dt = [('vertices', np.dtype(int)), ('primal_vertices', np.dtype(int))]
        structured_array = np.zeros(len(vertices), dtype=dt)
        structured_array['vertices'] = vertices
        structured_array['primal_vertices'] = primal_id_vertices
        structured_array = np.sort(structured_array, order='primal_vertices')
        vertices = structured_array['vertices']
        primal_id_vertices = structured_array['primal_vertices']

        nr = int(tol_n2*(len(vertices)-len(finos))/(Ni))
        n1 = self.data_impress['LEVEL_ID_1'].max() + 1
        n2 = self.data_impress['LEVEL_ID_2'].max() + 1
        nr = 10

        accum_levels = []

        pseudo_erro=np.repeat(TOL+1,2) #iniciou pseudo_erro
        t0=time.time()
        cont=0
        pos_new_inter=[]
        interm=np.array([])
        continuar = True


        while ((pseudo_erro.max()>TOL and n2<Nmax and iterar_mono and continuar) or cont==0) and cont<1:

            if cont>0:
                levels = self.data_impress['LEVEL'].copy()
                n1_ant = self.data_impress['LEVEL_ID_1'].max() + 1
                n2_ant = self.data_impress['LEVEL_ID_2'].max() + 1

                lim=np.sort(psr)[len(psr)-nr-1]
                positions=np.where(psr>lim)[0]

                # nv_verts=levels[vertices]
                # nv_positions=nv_verts[positions]
                # pos_new_fines=positions[nv_positions==1]
                # pos_new_inter=positions[nv_positions==2]

                # pos_new_fines=positions
                # interm=np.concatenate([interm,pos_new_inter]).astype(np.int)

                finos=np.concatenate([finos,positions]).astype(np.int)
                interm = np.setdiff1d(np.concatenate(M.volumes.bridge_adjacencies(finos, 2, 3)), finos)
                # primal_id_interm = np.unique(GID_1[interm])
                # interm = np.concatenate([GID_0[GID_1==k] for k in primal_id_interm])
                # primal_id_finos = np.unique(GID_1[finos])
                # finos = np.concatenate([GID_0[GID_1==k] for k in primal_id_finos])

                self.restart_levels()
                levels = self.data_impress['LEVEL'].copy()
                levels[finos] = 0
                levels[interm] = 1
                self.data_impress['LEVEL'] = levels.copy()
                gid_0 = self.data_impress['GID_0'][self.data_impress['LEVEL']==0]
                gid_1 = self.data_impress['GID_0'][self.data_impress['LEVEL']==1]
                self.set_adm_mesh_non_nested(gid_0, gid_1)
                # self.set_adm_mesh()
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
            psr=abs(pseudo_erro)
            psr[finos]=0

            perro.append(abs((SOL_ADM_fina-x0)/x0).max())
            active_nodes.append(n2/nfine_vols)

            if imprimir_a_cada_iteracao:
            #     # M1.mb.tag_set_data(Pseudo_ERRO_tag,M1.all_volumes,abs(pseudo_erro/x0)[GIDs])
            #     #
            #     # M1.mb.tag_set_data(ERRO_tag,M1.all_volumes,abs((SOL_ADM_fina-SOL_TPFA)/SOL_TPFA)[GIDs])
            #     # M1.mb.tag_set_data(P_ADM_tag,M1.all_volumes,SOL_ADM_fina[GIDs])
            #     # M1.mb.tag_set_data(P_TPFA_tag,M1.all_volumes,SOL_TPFA[GIDs])
            #     # ext_vtk = 'testes_MAD'  + str(cont) + '.vtk'
            #     # M1.mb.write_file(ext_vtk,[av])

                self.plot_operator(OP_ADM, OP_AMS_1, 0)
                self.data_impress.update_variables_to_mesh()
            #     # import pdb; pdb.set_trace()
            #     M.core.print(file='testt'+ str(cont), extension='.vtk', config_input='input_cards/print_settings0.yml')
            cont+=1

            accum_levels.append(self.data_impress['LEVEL'].copy())


        plt.plot(active_nodes,perro, marker='o')
        plt.yscale('log')
        plt.savefig('results/initial_adm_mesh/hist.png')



        # n = int(input('\nQual a malha adm que deseja utilizar?\nDigite o numero da iteracao.\n'))
        n=0

        if len(accum_levels)>0:
            self.data_impress['INITIAL_LEVEL'] = accum_levels[n]
        else:
            self.data_impress['INITIAL_LEVEL'] = self.data_impress['LEVEL'].copy()
        import pdb; pdb.set_trace()
        self.data_impress.update_variables_to_mesh()
        self.data_impress.export_to_npz()

    def set_pms_flux(self, wells, neumann_subds):

        master = masterNeumanNonNested(self.mesh, self.data_impress, self.elements_lv0, self.ml_data, self.n_levels, wells, neumann_subds)
        ms_flux_faces, pcorr = master.run()
        # del master
        self.data_impress['flux_faces'] = ms_flux_faces
        self.data_impress['pcorr'] = pcorr

    def plot_operator(self, OP_ADM, OP_AMS, v):

        vertices = self.elements_lv0['volumes'][self.data_impress['DUAL_1']==3]
        # primal_ids_vertices = self.data_impress['GID_1'][vertices]

        tags_adm = []
        tags_ams = []
        for i, v in enumerate(vertices):
            primal_id = self.data_impress['GID_1'][v]
            corresp = self.data_impress['ADM_COARSE_ID_LEVEL_1'][v]

            tags_ams.append(self.mesh.core.mb.tag_get_handle("OP_AMS"+str(i), 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True))
            tags_adm.append(self.mesh.core.mb.tag_get_handle("OP_ADM"+str(i), 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True))

            fb_ams = OP_AMS[:, primal_id].toarray()
            fb_adm = OP_ADM[:, corresp].toarray()

            self.mesh.core.mb.tag_set_data(tags_ams[i], self.mesh.core.all_volumes, fb_ams)
            self.mesh.core.mb.tag_set_data(tags_adm[i], self.mesh.core.all_volumes, fb_adm)



        # col = self.data_impress['GID_1'][self.all_wells_ids[v]]
        # fb_ms = OP_AMS[:, col].toarray()
        # corresp = self.data_impress['ADM_COARSE_ID_LEVEL_1'][self.all_wells_ids[v]]
        # fb_adm = OP_ADM[:, corresp].toarray()
        # self.data_impress['verif_po'] = fb_ms
        # self.data_impress['verif_rest'] = fb_adm

        self.print_test()
