import numpy as np
import scipy.sparse as sp
from ...common_files.common_infos import CommonInfos
import multiprocessing as mp
from ...solvers.solvers_scipy.solver_sp import SolverSp

from ...solvers.solvers_trilinos.solvers_tril import solverTril
import time

class masterNeumanNonNested:

    def __init__(self,M, data_impress, elements_lv0, ml_data, n_levels, T_without, wells, pare=False):
        self.data_impress = data_impress
        self.elements_lv0 = elements_lv0
        self.ml_data = ml_data
        self.n_levels = n_levels
        self.T_without = T_without
        self.wells = wells
        self.pare = pare
        self.one_worker = True
        self.mesh = M

    def get_n_workers(self, list_of_subdomains):

        if self.one_worker:
            n_cpu = 1
        else:
            n_cpu = mp.cpu_count()//2 - 1

        self.n_workers = n_cpu
        list_of_process_per_cpu = []

        n_subdomains = len(list_of_subdomains)
        resto = n_subdomains % self.n_workers
        n_process_per_cpu = n_subdomains//self.n_workers

        if n_process_per_cpu > 0:

        	for i in range(self.n_workers):
        		list_of_process_per_cpu.append(list_of_subdomains[i*n_process_per_cpu:n_process_per_cpu*(i+1)])

        	if resto != 0:
        		for i in range(resto):
        			list_of_process_per_cpu[i].append(list_of_subdomains[-i])

        else:
        	self.n_workers = resto

        	for i in range(resto):
        		list_of_process_per_cpu[i].append(list_of_subdomains[-i])

        return list_of_process_per_cpu

    def get_subdomains(self):

        '''
            ordem de envio:

                volumes: global id dos volumes locais
                ind_diric: indice dos volumes com pressao prescrita
                ind_neum: indice dos volumes com vazao prescrita
                val_diric: valores de pressao prescrita
                val_neum: valores de vazao prescrita
                local_transm: transmissibilidade local

                all_faces: todas faces do coarse volume
                intern_faces: faces internas do coarse volume
                intersect_faces: faces na interseccao
        '''

        list_of_subdomains = []
        pms_flux_faces = np.zeros(len(self.elements_lv0['faces']))
        levels = self.data_impress['LEVEL']
        pms = self.data_impress['pms']
        remaped_internal_faces = self.elements_lv0['remaped_internal_faces']
        neig_internal_faces = self.elements_lv0['neig_internal_faces']
        gid0 = self.data_impress['GID_0']
        n_volumes = len(levels)
        vols_lv0 = set(gid0[levels==0])

        if self.pare:
            import pdb; pdb.set_trace()
        self.data_impress['val_diric'][:]=0
        self.data_impress['val_neum'][:]=0
        for level in range(1, self.n_levels):
            str_level = str(level)
            set_level = set([level])

            all_gids_coarse = self.data_impress['GID_'+ str_level]
            # all_local_ids_coarse = self.data_impress['COARSE_LOCAL_ID_'+ str_level]
            all_intern_boundary_volumes = self.ml_data['internal_boundary_fine_volumes_level_'+ str_level]
            all_intersect_faces = self.ml_data['coarse_intersect_faces_level_'+ str_level]
            all_intern_faces = self.ml_data['coarse_internal_faces_level_'+ str_level]
            all_faces = self.ml_data['coarse_faces_level_'+ str_level]
            all_fine_vertex = self.ml_data['fine_vertex_coarse_volumes_level_'+ str_level]
            coarse_ids = self.ml_data['coarse_primal_id_level_'+ str_level]
            gids_level = np.unique(all_gids_coarse)

            for gidc in gids_level:

                intersect_faces = all_intersect_faces[coarse_ids==gidc][0] # faces na interseccao
                intern_local_faces = all_intern_faces[coarse_ids==gidc][0] # faces internas
                faces = all_faces[coarse_ids==gidc][0] # faces do volume
                intern_boundary_volumes = all_intern_boundary_volumes[coarse_ids==gidc][0] # volumes internos no contorno
                vertex = all_fine_vertex[coarse_ids==gidc]
                pressure_vertex = pms[vertex]
                volumes = self.elements_lv0['volumes'][all_gids_coarse==gidc]
                level_volumes = levels[volumes]
                volumes_dirichlet = set(volumes) & set(self.wells['ws_p'])
                volumes_neuman = set(volumes) & set(self.wells['ws_q'])

                adjs_intersect_faces = neig_internal_faces[remaped_internal_faces[intersect_faces]]
                adj_intern_local_faces = neig_internal_faces[remaped_internal_faces[intern_local_faces]]
                v0_new = adjs_intersect_faces.copy()
                intersect_faces_new = intersect_faces.copy()
                intern_boundary_volumes_new = intern_boundary_volumes.copy()

                ind_diric = []
                ind_neum = []
                val_diric = []
                val_neum = []

                if volumes_dirichlet:
                    ind_diric += list(volumes_dirichlet)
                    for v in ind_diric:
                        ## adicionar em ind_dirich os valores de pressao prescrita
                        val_diric += [self.wells['values_p'][self.wells['ws_p']==v][0]]
                        ## pegando os indices nas faces de interseccao local e removendo caso haja volumes de dirichlet
                        inds = ~((v0_new[:,0]==v) | (v0_new[:,1]==v))
                        intersect_faces_new = intersect_faces_new[inds]
                        intern_boundary_volumes_new = intern_boundary_volumes_new[~(intern_boundary_volumes_new==v)]

                volumes_dirichlet_2 = (set(volumes) & vols_lv0) - set(self.wells['ws_p'])
                if volumes_dirichlet_2:
                    for v in volumes_dirichlet_2:
                        ind_diric.append(v)
                        val_diric.append(pms[v])
                        inds = ~((v0_new[:,0]==v) | (v0_new[:,1]==v))
                        v0_new = v0_new[inds]
                        intersect_faces_new = intersect_faces_new[inds]
                        intern_boundary_volumes_new = intern_boundary_volumes_new[~(intern_boundary_volumes_new==v)]

                if volumes_neuman:
                    ind_neum += list(volumes_neuman)
                    for v in ind_neum:
                        val_neum += [self.wells['values_q'][self.wells['ws_q']==v][0]]
                        inds = ~((v0_new[:,0]==v) | (v0_new[:,1]==v))
                        v0_new = v0_new[inds]
                        intersect_faces_new = intersect_faces_new[inds]
                        intern_boundary_volumes_new = intern_boundary_volumes_new[~(intern_boundary_volumes_new==v)]

                if len(intern_boundary_volumes_new) > 0:
                    v0 = v0_new
                    pms0 = pms[v0[:,0]]
                    pms1 = pms[v0[:,1]]
                    t0 = self.data_impress['transmissibility'][intersect_faces_new]
                    pms_flux_faces_local = get_flux_faces(pms1, pms0, t0)
                    pms_flux_faces[intersect_faces_new] = pms_flux_faces_local

                    lines = np.concatenate([v0[:, 0], v0[:, 1]])
                    cols = np.repeat(0, len(lines))
                    data = np.concatenate([pms_flux_faces_local, -pms_flux_faces_local])
                    flux_pms_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
                    presc_flux_intern_boundary_volumes = flux_pms_volumes[intern_boundary_volumes_new]

                    ind_neum += list(intern_boundary_volumes_new)
                    val_neum += list(presc_flux_intern_boundary_volumes)

                if len(ind_diric) == 0:
                    if set(vertex) & set(ind_neum):
                        candidatos = set(volumes) - set(ind_neum)
                        vol = candidatos.pop()
                        ind_diric += [vol]
                        val_diric += [pms[vol]]
                    else:
                        ind_diric += list(vertex)
                        val_diric += list(pressure_vertex)

                # ind_diric=volumes[level_volumes==1]
                # if len(ind_diric)>0:
                #
                #     ind_diric=np.setdiff1d(ind_diric,intern_boundary_volumes_new)
                #     if len(ind_diric)>0:
                #         ind_diric=ind_diric[0]
                #     else:
                #         ind_diric=vertex
                # else:
                ind_diric=vertex
                val_diric=pms[ind_diric]
                self.data_impress['val_diric'][ind_diric]=val_diric
                self.data_impress['val_neum'][ind_neum]=val_neum

                list_of_subdomains.append(Subdomain(volumes, ind_diric, ind_neum, val_diric, val_neum, intern_local_faces, adj_intern_local_faces, self.T_without))

        return list_of_subdomains, pms_flux_faces

    def get_subdomains_2(self):
        '''
            ordem de envio:

                volumes: global id dos volumes locais
                ind_diric: indice dos volumes com pressao prescrita
                ind_neum: indice dos volumes com vazao prescrita
                val_diric: valores de pressao prescrita
                val_neum: valores de vazao prescrita
                local_transm: transmissibilidade local

                all_faces: todas faces do coarse volume
                intern_faces: faces internas do coarse volume
                intersect_faces: faces na interseccao
        '''
        list_of_subdomains = []
        pms_flux_faces = np.zeros(len(self.elements_lv0['faces']))
        levels = self.data_impress['LEVEL']
        pms = self.data_impress['pms']
        remaped_internal_faces = self.elements_lv0['remaped_internal_faces']
        neig_internal_faces = self.elements_lv0['neig_internal_faces']
        gid0 = self.data_impress['GID_0']
        n_volumes = len(gid0)
        for level in range(1, self.n_levels):
            str_level = str(level)
            set_level = set([level])

            all_gids_coarse = self.data_impress['GID_'+ str_level]
            all_intern_boundary_volumes = self.ml_data['internal_boundary_fine_volumes_level_'+ str_level]
            all_intersect_faces = self.ml_data['coarse_intersect_faces_level_'+ str_level]
            all_intern_faces = self.ml_data['coarse_internal_faces_level_'+ str_level]
            all_faces = self.ml_data['coarse_faces_level_'+ str_level]
            all_fine_vertex = self.ml_data['fine_vertex_coarse_volumes_level_'+ str_level]
            coarse_ids = self.ml_data['coarse_primal_id_level_'+ str_level]
            gids_level = np.unique(all_gids_coarse)

            import pdb; pdb.set_trace()

    def preprocess(self):
        list_of_subdomains, self.global_ms_flux_faces = self.get_subdomains()
        list_of_process_per_cpu = self.get_n_workers(list_of_subdomains)
        return list_of_process_per_cpu

    def run(self):
        list_of_process_per_cpu = self.preprocess()
        master2worker = [mp.Pipe() for _ in range(self.n_workers)]
        m2w, w2m = list(zip(*master2worker))
        procs = [mp.Process(target=run_thing, args=[LocalSolution(obj, comm)]) for obj, comm in zip(list_of_process_per_cpu, w2m)]
        del list_of_process_per_cpu
        global_pcorr = np.zeros(len(self.data_impress['GID_0']))

        for proc in procs:
        	proc.start()

        for comm in m2w:
            msg = comm.recv()
            for resp in msg:
                faces = resp[0]['faces']
                ms_flux = resp[0]['ms_flux_faces']
                self.global_ms_flux_faces[faces] = ms_flux
                volumes = resp[1]['volumes']
                pcorr = resp[1]['pcorr']
                global_pcorr[volumes] = pcorr

        for proc in procs:
        	proc.join()

        return self.global_ms_flux_faces.copy(), global_pcorr


def get_flux_faces(p1, p0, t0, flux_grav_faces=None):

    if flux_grav_faces != None:
        flux = -((p1 - p0) * t0 - flux_grav_faces)
    else:
        flux = -((p1 - p0) * t0)

    return flux

class Subdomain(CommonInfos):

    def __init__(self, volumes, ind_diric, ind_neum, val_diric, val_neum,
        intern_faces, adjs_intern_faces, T_global):

        self.T_local = self.get_local_t(T_global, volumes).tolil()
        self.volumes = volumes
        self.ids_local = np.arange(len(volumes))
        self.ind_diric = ind_diric
        self.ind_neum = ind_neum
        self.val_diric = val_diric
        self.val_neum = val_neum
        self.intern_faces = intern_faces
        self.adjs_intern_faces = adjs_intern_faces
        self.map_gid_in_lid = np.repeat(-1, volumes.max()+1)
        self.map_gid_in_lid[volumes] = self.ids_local

class LocalSolution:

    def __init__(self, subdomains, comm):
        self.subdomains = subdomains
        self.comm = comm

    def run(self):

        data = []
        dt = [('faces', int), ('ms_flux_faces', float)]
        dt_vol = [('volumes', int), ('pcorr', float)]
        solver = solverTril()

        for subd in self.subdomains:
            volumes = subd.volumes
            T_local = subd.T_local
            ids_local = subd.ids_local
            ind_diric = subd.ind_diric
            ind_neum = subd.ind_neum
            val_diric = subd.val_diric
            val_neum = subd.val_neum
            intern_faces = subd.intern_faces
            adjs_intern_faces = subd.adjs_intern_faces
            map_gid_in_lid = subd.map_gid_in_lid

            ind_diric_local = map_gid_in_lid[ind_diric]
            T_local_2 = T_local.copy()
            T_local_2[ind_diric_local] = 0
            T_local_2[ind_diric_local, ind_diric_local] = 1

            b = np.zeros(len(volumes))

            b[map_gid_in_lid[ind_neum]] = val_neum
            b[map_gid_in_lid[ind_diric]] = val_diric

            T_local_2 = T_local_2.tocsc()
            x = solver.solve_linear_problem(T_local_2,b)

            # x=a.solve(b)
            # print('\n')
            # print('pcorr')
            # print(x)
            # print('val_diric')
            # print(val_diric)
            # print('\n')
            del T_local_2

            t0 = T_local[map_gid_in_lid[adjs_intern_faces[:,0]], map_gid_in_lid[adjs_intern_faces[:,1]]].toarray().flatten()
            p0 = x[map_gid_in_lid[adjs_intern_faces[:,0]]]
            p1 = x[map_gid_in_lid[adjs_intern_faces[:,1]]]
            ms_flux = get_flux_faces(p1, p0, t0)

            sarray = np.zeros(len(intern_faces), dtype=dt)
            sarray['faces'] = intern_faces
            sarray['ms_flux_faces'] = ms_flux
            sarray_vol = np.zeros(len(volumes), dtype=dt_vol)
            sarray_vol['volumes'] = volumes
            sarray_vol['pcorr'] = x

            data.append([sarray, sarray_vol])

        self.comm.send(data)

def run_thing(local_solution_obj):
    local_solution_obj.run()
