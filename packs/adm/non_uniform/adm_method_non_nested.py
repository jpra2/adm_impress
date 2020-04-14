from ..adm_method import AdmMethod, Jacobi
import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg
import time
from ...directories import file_adm_mesh_def
import matplotlib.pyplot as plt
from pymoab import types

class AdmNonNested(AdmMethod):

    def set_adm_mesh_non_nested(self, v0=[], v1=[]):

        print("\nINICIOU GERACAO DA MALHA ADM\n")

        levels = np.repeat(-1, len(self.data_impress['LEVEL']))
        gids_0 = self.data_impress['GID_0']
        gids_1 = self.data_impress['GID_1']
        gids_2 = self.data_impress['GID_2']
        v2=np.setdiff1d(gids_0,np.concatenate([v0, v1]))

        levels[v0]=0
        levels[v1]=1
        levels[v2]=2


        if self.so_nv1==True:
            v1=np.setdiff1d(np.arange(len(levels)),v0)
        else:
            v1 = np.setdiff1d(np.concatenate(self.mesh.volumes.bridge_adjacencies(v0, 2, 3)), v0)
        levels[v1] = 1
        self.data_impress['LEVEL'] = levels.copy()

        n1 = 0
        n2 = 0

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
                # list_L1_ID[vols_ms1_lv1] = np.repeat(n1,len(vols_ms1_lv1))

                self.data_impress['ADM_COARSE_ID_LEVEL_1'][vols1] = np.repeat(n1, len(vols1))
                # n1+=1
                vols_ms2_lv1 = vols1[levels_vols_1==1]
                if len(vols_ms2_lv1)>0:
                    list_L2_ID[vols_ms2_lv1] = np.repeat(n2,len(vols_ms2_lv1))
                    list_L1_ID[vols_ms1_lv1] = np.repeat(n1,len(vols_ms1_lv1))
                    n1+=1
                    n2+=1

        self.data_impress['LEVEL_ID_1'] = list_L1_ID
        self.data_impress['LEVEL_ID_2'] = list_L2_ID

        for i in range(self.n_levels):
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


        while (pseudo_erro.max()>TOL and n2<Nmax and iterar_mono and continuar) and cont==0:

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

                np.savetxt("results/OP_ADM.csv",OP_ADM.toarray(),delimiter=",")
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

            # if imprimir_a_cada_iteracao:
            #     # M1.mb.tag_set_data(Pseudo_ERRO_tag,M1.all_volumes,abs(pseudo_erro/x0)[GIDs])
            #     #
            #     # M1.mb.tag_set_data(ERRO_tag,M1.all_volumes,abs((SOL_ADM_fina-SOL_TPFA)/SOL_TPFA)[GIDs])
            #     # M1.mb.tag_set_data(P_ADM_tag,M1.all_volumes,SOL_ADM_fina[GIDs])
            #     # M1.mb.tag_set_data(P_TPFA_tag,M1.all_volumes,SOL_TPFA[GIDs])
            #     # ext_vtk = 'testes_MAD'  + str(cont) + '.vtk'
            #     # M1.mb.write_file(ext_vtk,[av])
            #     self.data_impress.update_variables_to_mesh()
            #     self.plot_operator(OP_ADM, OP_AMS_1, 0)
            #     # import pdb; pdb.set_trace()
            #     M.core.print(file='testt'+ str(cont), extension='.vtk', config_input='input_cards/print_settings0.yml')
            cont+=1

            accum_levels.append(self.data_impress['LEVEL'].copy())


        plt.plot(active_nodes,perro, marker='o')
        plt.yscale('log')
        plt.savefig('results/initial_adm_mesh/hist.png')



        # n = int(input('\nQual a malha adm que deseja utilizar?\nDigite o numero da iteracao.\n'))
        n=0

        self.data_impress['INITIAL_LEVEL'] = accum_levels[n]

        self.data_impress.update_variables_to_mesh()
        self.data_impress.export_to_npz()

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
