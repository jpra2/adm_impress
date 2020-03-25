from ..adm_method import AdmMethod
import numpy as np
import scipy.sparse as sp

class AdmNonNested(AdmMethod):

    def set_adm_mesh_non_nested(self, v0=[], v1=[]):

        levels = np.repeat(-1, len(self.data_impress['LEVEL']))
        gids_0 = self.data_impress['GID_0']
        gids_1 = self.data_impress['GID_1']
        gids_2 = self.data_impress['GID_2']
        v2=np.setdiff1d(gids_0,np.concatenate([v0, v1]))

        levels[v0]=np.repeat(0, len(v0))
        levels[v1]=np.repeat(1, len(v1))
        levels[v2]=np.repeat(2, len(v2))

        all_wells = self.all_wells_ids
        if self.so_nv1==True:
            v1=np.setdiff1d(np.arange(len(levels)),all_wells)
        else:
            v1 = np.setdiff1d(np.concatenate(self.mesh.volumes.bridge_adjacencies(all_wells, 2, 3)), all_wells)
        levels[v1] = np.repeat(1, len(v1))

        n1 = 0
        n2 = 0
        self.data_impress['LEVEL'] = levels.copy()

        # vols_level_0 = gid_0[levels==0]
        # for n_level in range(1, self.n_levels):
        #     n = 0
        #     self.data_impress['LEVEL_ID_'+str(n_level)][vols_level_0] = np.arange(len(vols_level_0))
        #     n += len(vols_level_0)
        #     if n_level > 1:
        #         adm_ant_id = self.data_impress['LEVEL_ID_'+str(n_level-1)][levels<n_level & levels>0]
        #         unique_adm_ant_id = np.unique(adm_ant_id)
        #         for idd in unique_adm_ant_id:
        #             gid = gid_0[self.data_impress['LEVEL_ID_'+str(n_level-1)]==idd]
        #             self.data_impress['LEVEL_ID_'+str(n_level)][gid] = np.repeat(n, len(gid))
        #             n += 1
        #
        #     all_ids_coarse_level = self.data_impress['GID_'+str(n_level)][levels>=n_level]
        #     meshsets_ids = np.unique(all_ids_coarse_level)
        #     # v_level = gids_0[levels>=n_level]
        #     for id_coarse in meshsets_ids:
        #         volumes_in_meshset = gid_0[all_ids_coarse_level == id_coarse]
        #         # v_meshset_level=np.intersect1d(volumes_in_meshset,v_level)
        #         self.data_impress['LEVEL_ID_'+str(n_level)][volumes_in_meshset] = np.repeat(n, len(volumes_in_meshset))
        #         n += 1

        n0 = len(levels)
        list_L1_ID = np.repeat(-1, n0)
        list_L2_ID = np.repeat(-1, n0)

        list_L1_ID[v0] = np.arange(len(v0))
        list_L2_ID[v0] = np.arange(len(v0))
        n1+=len(v0)
        n2+=len(v0)

        ids_ms_2 = range(len(np.unique(gids_2)))


        print('\n')
        print("INICIOU GERACAO DA MALHA ADM")
        print('\n')

        for vol2 in ids_ms_2:
            #1
            # n_vols_l3 = 0
            vols2 = gids_0[gids_2==vol2]
            levels_vols_2 = levels[vols2]
            vols_ms2_lv2 = vols2[levels_vols_2==2]
            list_L2_ID[vols_ms2_lv2] = np.repeat(n2,len(vols_ms2_lv2))
            self.data_impress['ADM_COARSE_ID_LEVEL_2'][vols2] = np.repeat(n2, len(vols2))
            if len(vols_ms2_lv2)>0:
                n2+=1

            # gids_1_1 = gids_1[gids_2==v2]
            gids_1_1 = gids_1[vols2]
            ids_ms_1 = np.unique(gids_1_1)

            for vol1 in ids_ms_1:
                #2
                # elem_by_L1 = mb.get_entities_by_handle(m1)
                vols1 = vols2[gids_1_1==vol1]
                levels_vols_1 = levels_vols_2[gids_1_1==vol1]
                vols_ms1_lv1 = vols1[levels_vols_1>=1]
                list_L1_ID[vols_ms1_lv1] = np.repeat(n1,len(vols_ms1_lv1))
                self.data_impress['ADM_COARSE_ID_LEVEL_1'][vols1] = np.repeat(n1, len(vols1))
                n1+=1

                vols_ms2_lv1 = vols1[levels_vols_1==1]
                if len(vols_ms2_lv1)>0:
                    list_L2_ID[vols_ms2_lv1] = np.repeat(n2,len(vols_ms2_lv1))
                    n2+=1


        self.data_impress['LEVEL_ID_1'] = list_L1_ID
        self.data_impress['LEVEL_ID_2'] = list_L2_ID

        for i in range(self.n_levels+1):
            self.number_vols_in_levels[i] = len(levels[levels==i])

        self.n1_adm = n1
        self.n2_adm = n2

    def organize_ops_adm(self, OP_AMS, OR_AMS, level, _pcorr=None):

        gid_0 = self.data_impress['GID_0']
        gid_level = self.data_impress['GID_' + str(level)]
        gid_ant = self.data_impress['GID_' + str(level-1)]
        level_id = self.data_impress['LEVEL_ID_' + str(level)]
        level_id_ant = self.data_impress['LEVEL_ID_' + str(level-1)]
        levels = self.data_impress['LEVEL']
        dual_id = self.data_impress['DUAL_'+str(level)]
        OP_AMS = OP_AMS.tolil()

        if level == 1:
            OP_ADM, OR_ADM, pcorr = self.organize_ops_adm_level_1(OP_AMS, OR_AMS, level, _pcorr=_pcorr)
            self._data[self.adm_op_n + str(level)] = OP_ADM
            self._data[self.adm_rest_n + str(level)] = OR_ADM
            self._data[self.pcorr_n+str(level-1)] = pcorr
            return 0

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

        OP_ADM = sp.csc_matrix((data,(lines,cols)),shape=(len(gid_0),n1_adm))

        cols = gid_0
        lines = level_id
        data = np.ones(len(lines))
        OR_ADM = sp.csc_matrix((data,(lines,cols)),shape=(n1_adm,len(gid_0)))

        if self.get_correction_term:
            pcorr = np.zeros(len(gid_0))
            pcorr[levels>0] = _pcorr[levels>0]
        else:
            pcorr = np.array([False])

        return OP_ADM, OR_ADM, pcorr

    def plot_operator(self, OP_ADM, OP_AMS, v):

        col = self.data_impress['GID_1'][self.all_wells_ids[0]]
        fb_ms = OP_AMS[:, col].toarray()
        corresp = self.data_impress['ADM_COARSE_ID_LEVEL_1'][self.all_wells_ids[0]]
        fb_adm = OP_ADM[:, corresp].toarray()
        self.data_impress['verif_po'] = fb_ms
        self.data_impress['verif_rest'] = fb_adm
        self.print_test()
