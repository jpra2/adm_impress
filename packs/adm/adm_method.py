# from ..data_class.data_manager import DataManager
import numpy as np

class AdmMethod:

    def __init__(self, all_wells_ids, n_levels, M, data_impress, load=False):
        self.mesh = M
        self.all_wells_ids = all_wells_ids
        self.n_levels = n_levels
        self.data_impress = data_impress
        self.number_vols_in_levels = np.zeros(self.n_levels+1, dtype=int)

        if load == False:
            self.set_initial_mesh()

    def set_level_wells(self):
        levels = np.repeat(-1, len(self.data_impress['GID_0']))
        levels[self.all_wells_ids] = np.zeros(len(self.all_wells_ids))
        self.data_impress['LEVEL'] = levels

    def set_adm_mesh(self):

        levels = self.data_impress['LEVEL']
        gids_0 = self.data_impress['GID_0']
        gids_1 = self.data_impress['GID_1']
        gids_2 = self.data_impress['GID_2']
        vvv2 = range(len(np.unique(gids_2)))
        n0 = len(levels)

        list_L1_ID = np.repeat(-1, n0)
        list_L2_ID = np.repeat(-1, n0)

        n1=0
        n2=0
        n_vols = 0
        # meshset_by_L2 = mb.get_child_meshsets(self.L2_meshset)
        print('\n')
        print("INICIOU GERACAO DA MALHA ADM")
        print('\n')
        import pdb; pdb.set_trace()

        for v2 in vvv2:
            #1
            # n_vols_l3 = 0
            nivel3 = True
            nivel2 = False
            nivel1 = False
            vols2 = gids_0[gids_2==v2]
            gids_1_1 = gids_1[gids_2==v2]
            vvv1 = range(len(np.unique(gids_1_1)))
            n_vols_2 = len(vols2)
            conj_vols_1 = set()

            for v1 in vvv1:
                #2
                # elem_by_L1 = mb.get_entities_by_handle(m1)
                vols1 = vols2[gids_1_1==v1]
                nn1 = len(vols1)
                # n_vols += nn1
                # n_vols_l3 += nn1
                levels_vols_1 = levels[vols1]
                set_verif = set(levels_vols_1)

                if set([0]) & set_verif: # se houver volumes no nivel 0
                    #3
                    # volumes.append(elem_by_L1)
                    # meshsets_nv1.add(m1)
                    conj_vols_1.add(v1)
                    nivel3 = False
                    nivel1 = True
                    level = 0
                    list_L1_ID[vols1] = np.arange(n1, n1+nn1)
                    # list_L1_ID.append(np.arange(n1, n1+nn1))
                    list_L2_ID[vols1] = np.arange(n2, n2+nn1)
                    # list_L2_ID.append(np.arange(n2, n2+nn1))
                    levels[vols1] = np.repeat(level, nn1)
                    # list_L3_ID.append(np.repeat(level, nn1))
                    n1 += nn1
                    n2 += nn1
                #2
                elif set([1]) & set_verif: # se houver volumes no nivel 1
                    #3
                    # volumes.append(elem_by_L1)
                    # meshsets_nv2.add(m1)
                    conj_vols_1.add(v1)
                    nivel3 = False
                    nivel2 = True
                    level = 1
                    list_L1_ID[vols1] = np.repeat(n1, nn1)
                    # list_L1_ID.append(np.repeat(n1, nn1))
                    list_L2_ID[vols1] = np.repeat(n2, nn1)
                    # list_L2_ID.append(np.repeat(n2, nn1))
                    levels[vols1] = np.repeat(level, nn1)
                    # list_L3_ID.append(np.repeat(level, nn1))
                    n1 += 1
                    n2 += 1
            #1
            if nivel3:
                #2
                level = 2
                for v1 in vvv1:
                    #3
                    vols1 = vols2[gids_1_1==v1]
                    nn1 = len(vols1)
                    # volumes.append(elem_by_L1)
                    list_L1_ID[vols1] = np.repeat(n1, nn1)
                    # list_L1_ID.append(np.repeat(n1, nn1))
                    n1 += 1
                #2
                list_L2_ID[vols2] = np.repeat(n2, n_vols_2)
                # list_L2_ID.append(np.repeat(n2, n_vols_l3))
                levels[vols2] = np.repeat(level, n_vols_2)
                # list_L3_ID.append(np.repeat(level, n_vols_l3))
                n2 += 1
            #1
            else:
                #2
                vols_1_fora = set(vvv1) - conj_vols_1
                if vols_1_fora:
                    #3
                    for v1 in vols_1_fora:
                        #4
                        vols1 = vols2[gids_1_1==v1]
                        nn1 = len(vols1)
                        level = 1
                        list_L1_ID[vols1] = np.repeat(n1, nn1)
                        # list_L1_ID.append(np.repeat(n1, nn1))
                        list_L2_ID[vols1] = np.repeat(n2, nn1)
                        # list_L2_ID.append(np.repeat(n2, nn1))
                        levels[vols1] = np.repeat(level, nn1)
                        # list_L3_ID.append(np.repeat(level, nn1))
                        n1 += 1
                        n2 += 1

        self.data_impress['LEVEL_ID_1'] = list_L1_ID
        self.data_impress['LEVEL_ID_2'] = list_L2_ID
        self.data_impress['LEVEL'] = levels

        for i in range(self.n_levels+1):
            self.number_vols_in_levels[i] = len(levels[levels==i])

        self.n1_adm = n1
        self.n2_adm = n2

    def set_initial_mesh(self):
        pass
