from ..data_class.data_manager import DataManager
import numpy as np
import scipy.sparse as sp
from ..solvers.solvers_scipy.solver_sp import SolverSp
from ..flux_calculation.flux_tpfa import TpfaFlux2
import multiprocessing as mp
from .local_solution import LocalSolution
from. obj_infos import InfosForProcess
from ..directories import data_loaded
import time


def get_levelantids_levelids(level_ids_ant, level_ids):

    gids2 = np.unique(level_ids_ant)
    level_ids2 = []
    for i in gids2:
        test = level_ids_ant==i
        level_id = np.unique(level_ids[test])
        if len(level_id) > 1:
            raise ValueError('erro get_level_id')
        level_ids2.append(level_id[0])

    level_ids2 = np.array(level_ids2)

    return gids2, level_ids2

def solve_local_local_problem(solver, neigh_intern_faces, transmissibility, volumes,
                                        indices_p, values_p, indices_q=[], values_q=[]):

    # t0 = transmissibility
    v02 = neigh_intern_faces
    indices_p2 = indices_p
    t0 = transmissibility
    n = len(volumes)
    # local_ids = np.arange(n)
    # map_volumes = dict(zip(volumes, local_ids))
    # v02 = np.zeros(v0.shape)
    # v02[:,0] = np.array([map_volumes[k] for k in v0[:, 0]])
    # v02[:,1] = np.array([map_volumes[k] for k in v0[:, 1]])
    # indices_p2 = np.array(map_volumes[k] for k in indices_p)
    # indices_q2 = np.array(map_volumes[k] for k in indices_p)

    lines = np.concatenate([v02[:, 0], v02[:, 1], v02[:, 0], v02[:, 1]])
    cols = np.concatenate([v02[:, 1], v02[:, 0], v02[:, 0], v02[:, 1]])
    data = np.concatenate([t0, t0, -t0, -t0])
    T = sp.csc_matrix((data, (lines, cols)), shape=(n, n))
    T = T.tolil()

    T[indices_p2] = sp.lil_matrix((len(indices_p2), n))
    T[indices_p2, indices_p2] = np.ones(len(indices_p2))

    b = np.zeros(n)
    b[indices_p] = values_p
    b[indices_q] += values_q
    x = solver(T.tocsc(), b)
    return x

def set_boundary_conditions(T: 'transmissibility matrix',
                            b: 'source term',
                            indices_q: 'indices flux prescription',
                            values_q: 'flux prescription',
                            indices_p: 'indices pressure prescription',
                            values_p: 'values pressure prescription'):
    n = T.shape[0]

    T = T.tolil()
    np = len(indices_p)
    nq = len(indices_q)

    T[indices_p] = sp.lil_matrix((np, n))
    T[indices_p, indices_p] = np.ones(np)

    b[indices_q] += values_q

    return T.tocsc(), b


class AdmMethod(DataManager, TpfaFlux2):

    def __init__(self, all_wells_ids, n_levels, M, data_impress, elements_lv0, load=False):
        data_name = 'AdmMethod.npz'
        super().__init__(data_name=data_name, load=load)
        self.mesh = M
        self.elements_lv0 = elements_lv0
        self.ml_data = M.multilevel_data
        self.all_wells_ids = all_wells_ids
        self.n_levels = n_levels
        self.data_impress = data_impress
        self.number_vols_in_levels = np.zeros(self.n_levels+1, dtype=int)
        gids_0 = self.data_impress['GID_0']
        self.data_impress['LEVEL_ID_0'] = gids_0.copy()
        self.solver = SolverSp()

        self.adm_op_n = 'adm_prolongation_level_'
        self.adm_rest_n = 'adm_restriction_level_'

        self.n_cpu = mp.cpu_count()
        self.n_workers = self.n_cpu

        if load == False:
            self.set_initial_mesh()

    def set_level_wells(self):
        self.data_impress['LEVEL'][self.all_wells_ids] = np.zeros(len(self.all_wells_ids))

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

        for v2 in vvv2:
            #1
            # n_vols_l3 = 0
            nivel3 = True
            nivel2 = False
            nivel1 = False
            vols2 = gids_0[gids_2==v2]
            # gids_1_1 = gids_1[gids_2==v2]
            gids_1_1 = gids_1[vols2]
            vvv1 = np.unique(gids_1_1)
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

    def restart_levels(self):
        self.data_impress['LEVEL'] = np.repeat(-1, len(self.data_impress['LEVEL']))

    def organize_ops_adm(self, OP_AMS, OR_AMS, level):

        gid_0 = self.data_impress['GID_0']
        gid_level = self.data_impress['GID_' + str(level)]
        gid_ant = self.data_impress['GID_' + str(level-1)]
        level_id = self.data_impress['LEVEL_ID_' + str(level)]
        level_id_ant = self.data_impress['LEVEL_ID_' + str(level-1)]
        levels = self.data_impress['LEVEL']
        OP_AMS = OP_AMS.tolil()

        n_adm = len(np.unique(level_id))
        n_adm_ant = len(np.unique(level_id_ant))

        gids_nivel_n_engrossados = gid_0[levels<level]
        classic_ids_n_engrossados = set(gid_ant[gids_nivel_n_engrossados])
        adm_ids_ant_n_engrossados = level_id_ant[gids_nivel_n_engrossados]
        adm_ids_level_n_engrossados = level_id[gids_nivel_n_engrossados]
        if level > 1:
            adm_ids_ant_n_engrossados, adm_ids_level_n_engrossados = get_levelantids_levelids(adm_ids_ant_n_engrossados, adm_ids_level_n_engrossados)

        lines_op = adm_ids_ant_n_engrossados
        cols_op = adm_ids_level_n_engrossados
        data_op = np.ones(len(adm_ids_ant_n_engrossados))

        adm_ids_ant_gids = level_id_ant
        adm_ids_level = level_id
        classic_ids_ant = gid_ant
        classic_ids_level = gid_level

        ams_to_adm_coarse = dict(zip(classic_ids_level, adm_ids_level))
        ams_to_adm_fine = dict(zip(classic_ids_ant, adm_ids_ant_gids))

        if level > 1:
            adm_ids_ant_gids, adm_ids_level = get_levelantids_levelids(adm_ids_ant_gids, adm_ids_level)

        lines_2_op = []
        cols_2_op = []
        data_2_op = []

        lines_or = adm_ids_level
        cols_or = adm_ids_ant_gids
        data_or = np.repeat(1.0, len(adm_ids_level))

        data_op_ams = sp.find(OP_AMS)

        for l, c, d, in zip(data_op_ams[0], data_op_ams[1], data_op_ams[2]):
            if set([l]) & classic_ids_n_engrossados:
                continue
            lines_2_op.append(ams_to_adm_fine[l])
            cols_2_op.append(ams_to_adm_coarse[c])
            data_2_op.append(d)

        lines_2_op = np.array(lines_2_op)
        cols_2_op = np.array(cols_2_op)
        data_2_op = np.array(data_2_op)

        lines_op = np.concatenate([lines_op, lines_2_op])
        cols_op = np.concatenate([cols_op, cols_2_op])
        data_op = np.concatenate([data_op, data_2_op])

        OP_ADM = sp.csc_matrix((data_op, (lines_op, cols_op)), shape=(n_adm_ant, n_adm))
        OR_ADM = sp.csc_matrix((data_or, (lines_or, cols_or)), shape=(n_adm, n_adm_ant))

        self._data[self.adm_op_n + str(level)] = OP_ADM
        self._data[self.adm_rest_n + str(level)] = OR_ADM

    def solve_multiscale_pressure(self, T: 'fine transmissibility matrix', b: 'fine source term'):

        T_adm = T.copy()
        b_adm = b.copy()

        n_levels = self.n_levels
        for i in range(n_levels):
            level = i+1
            # op_adm = self._data[self.adm_op_n + str(level)]
            # rest_adm = self._data[self.adm_rest_n + str(level)]
            T_adm = self._data[self.adm_rest_n + str(level)]*T_adm*self._data[self.adm_op_n + str(level)]
            b_adm = self._data[self.adm_rest_n + str(level)]*b_adm

        pms = self.solver.direct_solver(T_adm, b_adm)
        # p_adm = pms.copy()

        for i in range(n_levels):
            level = self.n_levels - i
            pms = self._data[self.adm_op_n + str(level)]*pms

        self.data_impress['pms'] = pms
        self.data_impress['pressure'] = pms
        self.T = T

    def set_pms_flux_intersect_faces(self):

        levels = self.data_impress['LEVEL']
        faces_intersect_lv1 = np.unique(np.concatenate(self.ml_data['coarse_intersect_faces_level_'+str(1)]))
        neig_intersect_faces_lv1 = self.ml_data['neig_intersect_faces_level_'+str(1)]

        v0 = neig_intersect_faces_lv1
        n_volumes = len(self.elements_lv0['volumes'])
        pms = self.data_impress['pms']
        t_intersect_faces = self.data_impress['transmissibility'][faces_intersect_lv1]
        t0 = t_intersect_faces
        flux_grav_intersect_faces = self.data_impress['flux_grav_faces'][faces_intersect_lv1]

        ps0 = pms[v0[:, 0]]
        ps1 = pms[v0[:, 1]]
        flux_intersect_faces = -((ps1 - ps0) * t0 - flux_grav_intersect_faces)

        lines = np.concatenate([v0[:, 0], v0[:, 1]])
        cols = np.repeat(0, len(lines))
        data = np.concatenate([flux_intersect_faces, -flux_intersect_faces])
        flux_pms_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()

        n = len(self.data_impress['pms_flux_faces'])
        flux = np.zeros(n)
        flux[faces_intersect_lv1] = flux_intersect_faces

        self.data_impress['pms_flux_faces'] = flux
        self.data_impress['pms_flux_interfaces_volumes'] = flux_pms_volumes

    def set_pcorr(self):
        presc_flux_volumes = self.data_impress['pms_flux_interfaces_volumes'].copy()
        levels = self.data_impress['LEVEL']
        n_volumes = len(levels)
        gid0 = self.data_impress['GID_0']
        transmissibility = self.data_impress['transmissibility']
        pms = self.data_impress['pms']
        neig_internal_faces = self.elements_lv0['neig_internal_faces']
        remaped_internal_faces = self.elements_lv0['remaped_internal_faces']
        flux_grav_faces = self.data_impress['flux_grav_faces']

        pcorr = np.zeros(len(pms))
        flux_faces = np.zeros(len(transmissibility))
        flux_volumes = np.zeros(n_volumes)

        for i in range(self.n_levels):
            level=i+1
            all_gids_coarse = self.data_impress['GID_'+str(level)]
            all_local_ids_coarse = self.data_impress['COARSE_LOCAL_ID_'+str(level)]
            all_intern_boundary_volumes = self.ml_data['internal_boundary_fine_volumes_level_'+str(level)]
            all_intersect_faces = self.ml_data['coarse_intersect_faces_level_'+str(level)]
            all_intern_faces = self.ml_data['coarse_internal_faces_level_'+str(level)]
            # all_faces = self.ml_data['coarse_faces_level_'+str(level)]
            all_fine_vertex = self.ml_data['fine_vertex_coarse_volumes_level_'+str(level)]
            coarse_ids = self.ml_data['coarse_primal_id_level_'+str(level)]
            gids_level = np.unique(all_gids_coarse[levels==level])
            for gidc in gids_level:
                intersect_faces = all_intersect_faces[coarse_ids==gidc][0] # faces na interseccao
                intern_local_faces = all_intern_faces[coarse_ids==gidc][0] # faces internas
                neig_internal_local_faces = neig_internal_faces[remaped_internal_faces[intern_local_faces]]
                intern_boundary_volumes = all_intern_boundary_volumes[coarse_ids==gidc][0] # volumes internos no contorno
                vertex = all_fine_vertex[coarse_ids==gidc]
                pressure_vertex = pms[vertex]
                volumes = gid0[all_gids_coarse==gidc]

                local_id_volumes = all_local_ids_coarse[volumes]
                local_neig_internal_local_faces = neig_internal_local_faces.copy()
                local_neig_internal_local_faces[:,0] = all_local_ids_coarse[neig_internal_local_faces[:,0]]
                local_neig_internal_local_faces[:,1] = all_local_ids_coarse[neig_internal_local_faces[:,1]]
                local_intern_boundary_volumes = all_local_ids_coarse[intern_boundary_volumes]
                values_q = presc_flux_volumes[intern_boundary_volumes]
                local_vertex = all_local_ids_coarse[vertex]
                t0 = transmissibility[intern_local_faces]
                x = solve_local_local_problem(self.solver.direct_solver, local_neig_internal_local_faces, t0, local_id_volumes,
                    local_vertex, pressure_vertex, local_intern_boundary_volumes, values_q)

                pcorr[volumes] = x

                neig_intersect_faces = neig_internal_faces[remaped_internal_faces[intersect_faces]]
                transmissibility_intersect_faces = transmissibility[intersect_faces]
                t0 = transmissibility_intersect_faces
                pms0 = pms[neig_intersect_faces[:,0]]
                pms1 = pms[neig_intersect_faces[:,1]]
                flux_grav_intersect_faces = flux_grav_faces[intersect_faces]
                flux_intersect_faces = -((pms1 - pms0) * t0 - flux_grav_intersect_faces)
                flux_faces[intersect_faces] = flux_intersect_faces

                pcorr0 = pcorr[neig_internal_local_faces[:,0]]
                pcorr1 = pcorr[neig_internal_local_faces[:,1]]
                flux_grav_intern_faces = flux_grav_faces[intern_local_faces]
                t0 = transmissibility[intern_local_faces]
                flux_intern_faces = -((pcorr1 - pcorr0) * t0 - flux_grav_intern_faces)
                flux_faces[intern_local_faces] = flux_intern_faces

                v0 = neig_internal_local_faces

                lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
                cols = np.repeat(0, len(lines))
                data = np.array([flux_intern_faces, -flux_intern_faces]).flatten()
                flux_volumes_2 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
                flux_volumes_2[intern_boundary_volumes] += values_q
                flux_volumes[volumes] = flux_volumes_2[volumes]

        volumes_fine = gid0[levels==0]
        intern_faces_volumes_fine = self.mesh.volumes.bridge_adjacencies(volumes_fine, 3, 2)
        intern_faces_volumes_fine = np.setdiff1d(intern_faces_volumes_fine, self.elements_lv0['boundary_faces'])
        neig_intern_faces_volumes_fine = neig_internal_faces[remaped_internal_faces[intern_faces_volumes_fine]]
        v0 = neig_intern_faces_volumes_fine

        pms0 = pms[neig_intern_faces_volumes_fine[:,0]]
        pms1 = pms[neig_intern_faces_volumes_fine[:,1]]
        t0 = transmissibility[intern_faces_volumes_fine]
        flux_grav_faces_volumes_fine = flux_grav_faces[intern_faces_volumes_fine]
        flux_intern_faces_volumes_fine = -((pms1 - pms0) * t0 - flux_grav_faces_volumes_fine)
        flux_faces[intern_faces_volumes_fine] = flux_intern_faces_volumes_fine

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.concatenate([flux_intern_faces_volumes_fine, -flux_intern_faces_volumes_fine])
        flux_volumes_2 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()

        flux_volumes[volumes_fine] = flux_volumes_2[volumes_fine]

        self.data_impress['pcorr'] = pcorr
        self.data_impress['flux_faces'] = flux_faces
        self.data_impress['flux_volumes'] = flux_volumes

        # #######################
        # ## test
        # v0 = neig_internal_faces
        # internal_faces = self.elements_lv0['internal_faces']
        # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        # cols = np.repeat(0, len(lines))
        # data = np.array([flux_faces[internal_faces], -flux_faces[internal_faces]]).flatten()
        # flux_volumes_2 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
        # self.data_impress['flux_volumes_test'] = flux_volumes_2
        # ######################################

    def get_nt_process(self):

        levels = self.data_impress['LEVEL']
        gids = self.data_impress['GID_0']
        nt_process = 0

        for i in range(self.n_levels):
            level = i+1
            gids_lv = self.data_impress['GID_' + str(level)]
            gids_eng = np.unique(gids_lv[levels==level])
            nt_process += len(gids_eng)

        return nt_process

    def get_lists_objects(self, infos):

        levels = self.data_impress['LEVEL']
        pms = self.data_impress['pms']
        gid0 = self.data_impress['GID_0']
        presc_flux_volumes = self.data_impress['pms_flux_interfaces_volumes']

        all_list_volumes = []
        all_list_indices_d = []
        all_list_values_d = []
        all_list_indices_n = []
        all_list_values_n = []
        all_list_faces = []
        all_list_internal_faces = []
        all_list_intersect_faces = []
        cont0 = 0
        cont = 0

        for i in range(self.n_levels):
            level=i+1
            st_level = str(level)
            all_gids_coarse = self.data_impress['GID_'+ st_level]
            all_local_ids_coarse = self.data_impress['COARSE_LOCAL_ID_'+ st_level]
            all_intern_boundary_volumes = self.ml_data['internal_boundary_fine_volumes_level_'+ st_level]
            all_intersect_faces = self.ml_data['coarse_intersect_faces_level_'+ st_level]
            all_intern_faces = self.ml_data['coarse_internal_faces_level_'+ st_level]
            all_faces = self.ml_data['coarse_faces_level_'+ st_level]
            all_fine_vertex = self.ml_data['fine_vertex_coarse_volumes_level_'+ st_level]
            coarse_ids = self.ml_data['coarse_primal_id_level_'+ st_level]
            gids_level = np.unique(all_gids_coarse[levels==level])
            for gidc in gids_level:
                intersect_faces = all_intersect_faces[coarse_ids==gidc][0] # faces na interseccao
                intern_local_faces = all_intern_faces[coarse_ids==gidc][0] # faces internas
                faces = all_faces[coarse_ids==gidc][0] # faces do volume
                intern_boundary_volumes = all_intern_boundary_volumes[coarse_ids==gidc][0] # volumes internos no contorno
                vertex = all_fine_vertex[coarse_ids==gidc]
                pressure_vertex = pms[vertex]
                volumes = gid0[all_gids_coarse==gidc]
                values_q = presc_flux_volumes[intern_boundary_volumes]

                if cont0 < self.n_workers and cont0 >= 0:
                    all_list_volumes.append([volumes])
                    all_list_indices_d.append([vertex])
                    all_list_values_d.append([pressure_vertex])
                    all_list_indices_n.append([intern_boundary_volumes])
                    all_list_values_n.append([values_q])
                    all_list_faces.append([faces])
                    all_list_internal_faces.append([intern_local_faces])
                    all_list_intersect_faces.append([intersect_faces])
                    cont0 += 1
                    if cont0 == self.n_workers:
                        cont0 = -1
                else:
                    all_list_volumes[cont].append(volumes)
                    all_list_indices_d[cont].append(vertex)
                    all_list_values_d[cont].append(pressure_vertex)
                    all_list_indices_n[cont].append(intern_boundary_volumes)
                    all_list_values_n[cont].append(values_q)
                    all_list_faces[cont].append(faces)
                    all_list_internal_faces[cont].append(intern_local_faces)
                    all_list_intersect_faces[cont].append(intersect_faces)
                    cont += 1
                    if cont == self.n_workers:
                        cont = 0

        list_objects = []
        for i in range(self.n_workers):
            list_objects.append(
                LocalSolution(
                    all_list_volumes[i],
                    all_list_indices_d[i],
                    all_list_values_d[i],
                    all_list_indices_n[i],
                    all_list_values_n[i],
                    all_list_faces[i],
                    all_list_internal_faces[i],
                    all_list_intersect_faces[i],
                    infos.copy()
                )
            )

        return list_objects

    def set_paralel_pcorr(self):
        transmissibility = self.data_impress['transmissibility']
        presc_flux_volumes = self.data_impress['pms_flux_interfaces_volumes']
        levels = self.data_impress['LEVEL']
        gids = self.data_impress['GID_0']
        pms = self.data_impress['pms']
        g_neig_internal_faces = self.elements_lv0['neig_internal_faces']
        remaped_internal_faces = self.elements_lv0['remaped_internal_faces']
        g_flux_grav_faces = self.data_impress['flux_grav_faces']
        g_faces = self.elements_lv0['faces']
        T = self.T
        solver = self.solver.direct_solver
        n_volumes = len(gids)

        _pcorr = np.zeros(n_volumes)
        _flux_faces = np.zeros(len(transmissibility))
        _flux_volumes = np.zeros(n_volumes)

        # nt_process = self.get_nt_process()
        # self.n_workers = nt_process
        # self.n_workers = 1

        # m = mp.Manager()
        # qvolumes = m.Queue()
        # qfaces = m.Queue()
        # qinfos = m.Queue()
        # lock = mp.Lock()

        infos = InfosForProcess(T, pms, g_flux_grav_faces, gids, g_faces, g_neig_internal_faces,
            remaped_internal_faces, solver)

        # for i in range(self.n_workers):
        #     qinfos.put(infos.copy())
        # tt = qinfos.qsize()

        list_objects = self.get_lists_objects(infos)

        def f(local_solution_obj, w2m):
            local_solution_obj.run(w2m)

        master2worker = [mp.Pipe() for _ in range(self.n_workers)]
        m2w, w2m = list(zip(*master2worker))
        procs = [mp.Process(target=f, args=[obj, comm]) for obj, comm in zip(list_objects, w2m)]

        # def f(local_solution_obj, qinfos, qvolumes, qfaces):
        #     local_solution_obj.run(qinfos, qvolumes, qfaces)
        #     # return 0

            # return 0

        # results = Parallel(n_jobs=self.n_workers, require='sharedmem')(delayed(f)(i) for i in list_objects)
        # procs = []

        # for i in range(self.n_workers):
        #     proc = mp.Process(target=f, args=(list_objects[i], qinfos, qvolumes, qfaces))
        #     procs.append(proc)

        for proc in procs:
            proc.start()

        for comm in m2w:
            msg = comm.recv()
            for resp in msg:
                resp_vols = resp[0]
                resp_faces = resp[1]
                _pcorr[resp_vols['volumes']] = resp_vols['pcorr']
                _flux_volumes[resp_vols['volumes']] = resp_vols['flux_volumes']
                _flux_faces[resp_faces['faces']] = resp_faces['flux_faces']

        for proc in procs:
            proc.join()

        # while not qvolumes.empty():
        #     resp = qvolumes.get()
        #     _pcorr[resp['volumes']] = resp['pcorr']
        #     _flux_volumes[resp['volumes']] = resp['flux_volumes']
        #
        # while not qfaces.empty():
        #     resp = qfaces.get()
        #     _flux_faces[resp['faces']] = resp['flux_faces']

        # for result in results:
        #     for resp in result:
        #         resp_vols = resp[0]
        #         resp_faces = resp[1]
        #         _pcorr[resp_vols['volumes']] = resp_vols['pcorr']
        #         _flux_volumes[resp_vols['volumes']] = resp_vols['flux_volumes']
        #         _flux_faces[resp_faces['faces']] = resp_faces['flux_faces']

        gid0 = gids
        volumes_fine = gid0[levels==0]
        intern_faces_volumes_fine = self.mesh.volumes.bridge_adjacencies(volumes_fine, 3, 2)
        intern_faces_volumes_fine = np.setdiff1d(intern_faces_volumes_fine, self.elements_lv0['boundary_faces'])
        neig_intern_faces_volumes_fine = g_neig_internal_faces[remaped_internal_faces[intern_faces_volumes_fine]]
        v0 = neig_intern_faces_volumes_fine

        pms0 = pms[neig_intern_faces_volumes_fine[:,0]]
        pms1 = pms[neig_intern_faces_volumes_fine[:,1]]
        t0 = transmissibility[intern_faces_volumes_fine]
        flux_grav_faces_volumes_fine = g_flux_grav_faces[intern_faces_volumes_fine]
        flux_intern_faces_volumes_fine = -((pms1 - pms0) * t0 - flux_grav_faces_volumes_fine)
        _flux_faces[intern_faces_volumes_fine] = flux_intern_faces_volumes_fine

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.concatenate([flux_intern_faces_volumes_fine, -flux_intern_faces_volumes_fine])
        flux_volumes_2 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
        _flux_volumes[volumes_fine] = flux_volumes_2[volumes_fine]

        self.data_impress['pcorr'] = _pcorr
        self.data_impress['flux_volumes'] = _flux_volumes
        self.data_impress['flux_faces'] = _flux_faces

        _debug = data_loaded['_debug']
        if _debug:

            #######################
            ## test
            v0 = g_neig_internal_faces
            internal_faces = self.elements_lv0['internal_faces']
            lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
            cols = np.repeat(0, len(lines))
            data = np.array([_flux_faces[internal_faces], -_flux_faces[internal_faces]]).flatten()
            flux_volumes_2 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
            self.data_impress['flux_volumes_test'] = flux_volumes_2
            ######################################

    def set_initial_mesh(self):
        # TODO: atualizar
        pass
