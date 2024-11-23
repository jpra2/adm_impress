import scipy.sparse as sp
import numpy as np
from ..directories import data_loaded
from ..errors.err import EmptyQueueError

def run_local_solution(local_solution_obj,
    T: 'global transmissibility matrix without boundary conditions',
    pms: 'global multiscale presure',
    g_flux_grav_faces,
    gids: 'global gids',
    g_faces: 'global_faces',
    g_neig_internal_faces: 'all neig internal faces',
    remaped_internal_faces,
    solver):
    return local_solution_obj.run(T, pms, g_flux_grav_faces, gids, g_faces,
           g_neig_internal_faces, remaped_internal_faces, solver)

def get_list_n_process(n_cpu, n_total_process):

    n_process_per_cpu = int(round(n_total_process/n_cpu))
    list_of_n_process = np.repeat(n_process_per_cpu, n_cpu)
    cont = 0

    if n_cpu*n_process_per_cpu > n_total_process:
        while list_of_n_process.sum() > n_total_process:
            list_of_n_process[cont] -= 1
            cont += 1
    elif n_cpu*n_process_per_cpu < n_total_process:
        while list_of_n_process.sum() < n_total_process:
            list_of_n_process[cont] += 1
            cont += 1

    return list_of_n_process


class LocalSolution:
    def __init__(self,
        volumes: 'list of gids of local volumes',
        indices_d: 'list of gids of indices with dirichlet prescription',
        values_d: 'list of pressure values',
        indices_n: 'list of gids of indices with neumman prescription',
        values_n: 'list of neumman values',
        faces: 'list of faces',
        internal_faces: 'list of internal faces',
        intersect_faces: 'list of intersect faces',
        infos):

        self.n_problems = len(volumes)
        self.volumes = volumes
        self.indices_d = indices_d
        self.values_d = values_d
        self.indices_n = indices_n
        self.values_n = values_n
        self.faces = faces
        self.internal_faces = internal_faces
        self.intersect_faces = intersect_faces
        self.infos = infos

    def get_remaped_gids(self, gids, volumes, g_faces, faces):
        n_vols = len(volumes)
        local_ids = np.arange(n_vols)
        remaped_gids = gids.copy()
        remaped_gids[volumes] = local_ids
        local_faces = np.arange(len(faces))
        remaped_faces = g_faces.copy()
        remaped_faces[faces] = local_faces

        return remaped_gids, remaped_faces

    def get_local_t(self, volumes, T):
        '''
            get the local transmissibility
        '''
        T2_w = T[volumes][:,volumes]
        data = np.array(T2_w.sum(axis=1).transpose())[0]
        diag = T2_w.diagonal()
        diag -= data
        T2_w.setdiag(diag)
        return  T2_w

    def get_local_problem(self, T2_w, indices_d, values_d, indices_n, values_n, remaped_gids):
        '''
            get the local transmissibility and source
            term with boundary conditions
        '''

        if data_loaded['_debug']:
            soma = abs(values_n.sum())
            if soma > 1e-5:
                raise ValueError(f'Volume nao conservativo {soma}, {values_n}')

        n1 = len(indices_d)
        T2 = T2_w.copy().tolil()
        T2[remaped_gids[indices_d]] = sp.lil_matrix((n1, n1))
        T2[remaped_gids[indices_d], remaped_gids[indices_d]] = np.ones(n1)
        n_vols = T2.shape[0]

        b = np.zeros(n_vols)
        if len(values_n) > 0:
            b[remaped_gids[indices_n]] += values_n
        b[remaped_gids[indices_d]] = values_d

        return T2.tocsc(), b

    def solve_local_problem(self, T2, b, solver):
        return solver(T2, b)

    def local_flux(self,
        T,
        T2_w,
        pcorr,
        pms,
        values_n,
        indices_n,
        remaped_gids,
        g_flux_grav_faces,
        internal_faces,
        g_neig_internal_faces,
        faces,
        remaped_faces,
        intersect_faces,
        remaped_internal_faces):

        # internal faces
        ids_l0 = remaped_gids[g_neig_internal_faces[remaped_internal_faces[internal_faces]][:, 0]]
        ids_l1 = remaped_gids[g_neig_internal_faces[remaped_internal_faces[internal_faces]][:, 1]]
        ps0 = pcorr[ids_l0]
        ps1 = pcorr[ids_l1]
        t0 = np.array(T2_w[ids_l0, ids_l1])[0]
        flux_faces = np.zeros(len(faces))
        n_volumes = T2_w.shape[0]
        grav_flux_int_faces = g_flux_grav_faces[internal_faces]

        flux_internal_faces = -((ps1 - ps0) * t0 - grav_flux_int_faces)

        lines = np.array([ids_l0, ids_l1]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.concatenate([flux_internal_faces, -flux_internal_faces])
        flux_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
        flux_volumes[remaped_gids[indices_n]] += values_n

        # intersect faces
        ids_g0 = g_neig_internal_faces[remaped_internal_faces[intersect_faces]][:,0]
        ids_g1 = g_neig_internal_faces[remaped_internal_faces[intersect_faces]][:,1]
        pms0 = pms[ids_g0]
        pms1 = pms[ids_g1]
        t0 = np.array(T[ids_g0, ids_g1])[0]
        grav_flux_intersect_faces = g_flux_grav_faces[intersect_faces]

        flux_intersect_faces = -((pms1 - pms0) * t0 - grav_flux_intersect_faces)

        flux_faces[remaped_faces[intersect_faces]] = flux_intersect_faces
        flux_faces[remaped_faces[internal_faces]] = flux_internal_faces

        return flux_faces, flux_volumes

    def run(self, w2m):

        # # lock.acquire()
        # try:
        #     # infos = qinfos.get(block=True, timeout=1)
        #     infos = qinfos.get()
        # except:
        #     raise EmptyQueueError('Empty Queue')
        # lock.release()

        infos = self.infos

        # solution = 0
        solution = []
        dtvolumes = [('volumes', np.dtype(int)), ('pcorr', np.dtype(float)), ('flux_volumes', np.dtype(float))]
        dtfaces = [('faces', np.dtype(int)), ('flux_faces', np.dtype(float))]

        for i in range(self.n_problems):
            volumes = self.volumes[i]
            indices_d = self.indices_d[i]
            indices_n = self.indices_n[i]
            values_d = self.values_d[i]
            values_n = self.values_n[i]
            faces = self.faces[i]
            internal_faces = self.internal_faces[i]
            intersect_faces = self.intersect_faces[i]
            sarray_vols = np.zeros(len(volumes), dtype=dtvolumes)
            sarray_faces = np.zeros(len(faces), dtype=dtfaces)

            remaped_gids, remaped_faces = self.get_remaped_gids(infos.gids, volumes, infos.g_faces, faces)
            T2_w = self.get_local_t(volumes, infos.T)
            T2, b = self.get_local_problem(T2_w, indices_d, values_d, indices_n, values_n, remaped_gids)
            pcorr = self.solve_local_problem(T2, b, infos.solver)
            flux_faces, flux_volumes = self.local_flux(infos.T, T2_w, pcorr, infos.pms, values_n, indices_n,
                                                       remaped_gids, infos.g_flux_grav_faces, internal_faces,
                                                       infos.g_neig_internal_faces, faces, remaped_faces,
                                                       intersect_faces, infos.remaped_internal_faces)
            sarray_vols['volumes'] = volumes
            sarray_vols['pcorr'] = pcorr
            sarray_vols['flux_volumes'] = flux_volumes
            sarray_faces['faces'] = faces
            sarray_faces['flux_faces'] = flux_faces
            solution.append(np.array([sarray_vols, sarray_faces]))
            # # lock.acquire()
            # qvolumes.put(sarray_vols)
            # qfaces.put(sarray_faces)
            # # lock.release()

        w2m.send(np.array(solution))
