from packs.solvers.solvers_scipy.solver_sp import SolverSp
from multiprocessing import Queue
import os
from packs.compositional.IMPEC.global_pressure_solver import GlobalIMPECPressureSolver as Gips
from packs.multiscale.ms_utils.multiscale_functions import update_local_transmissibility, set_matrix_pressure_prescription

def run_thing(local_solver_obj):
    local_solver_obj.run()

def run_thing2(local_solver_obj):
    local_solver_obj.run2()


class InitAndFinishProcess:
    def finish(self):
        print(f'\nProcess {self.id_process} finished \n')
        self.finished.value = True

    def initialize(self):
        print(f'\nProcess {self.id_process} start \n')
    


class GlobalLocalSolver(InitAndFinishProcess):
    def __init__(self, subdomains, queue: Queue, comm, finished, id_process, **kwargs):
        """
        @param subdomains: list of subdomains
        @param queue: queue object
        @param comm: communicator
        @param finished: shared boolean value: True = process finished
        @param id_process: process id
        """
        self.subdomains = subdomains
        self.queue = queue
        self.comm = comm
        self.finished = finished
        self.id_process = id_process
        self.update_FC = kwargs.get('update_FC')
        self.global_vector_update = kwargs.get('global_vector_update')
        self.T_fine_without_bc = kwargs.get('T_fine_without_bc')
        self.global_diagonal_term = kwargs.get('global_diagonal_term')
        self.global_source_term = kwargs.get('global_source_term')


class LocalSolver1(GlobalLocalSolver):

    def run(self):
        self.initialize()
        solver = SolverSp()
        for subd in self.subdomains:
            resp = solver.direct_solver(subd.Tlocal, subd.local_rhs)
            self.queue.put([subd.volumes, resp])

        self.finish()


class LocalSolver2(InitAndFinishProcess):
    
    def __init__(self, subdomains, queue: Queue, comm, finished, id_process, local_params, **kwargs):
        """
        @param subdomains: list of subdomains
        @param queue: queue object
        @param comm: communicator
        @param finished: shared boolean value: True = process finished
        @param id_process: process id
        """
        self.subdomains = subdomains
        self.queue = queue
        self.comm = comm
        self.finished = finished
        self.id_process = id_process
        
        self.Tglobal = kwargs['fine_scale_transmissibility_no_bc']
        self.map_internal_faces = kwargs['map_internal_faces'] ## equals to elements_lv0['remaped_internal_faces]
        self.adm_pressure = kwargs['adm_pressure']
        self.diagonal_term = kwargs['diagonal_term']
        self.global_molar_flux_prescription = kwargs['global_molar_flux_prescription']
        
        self.params = local_params
    
    def run(self):
        
        self.initialize()
        solver = SolverSp()
        
        for subd in self.subdomains:
            subd.Tlocal_no_bc = update_local_transmissibility(self.Tglobal, subd.volumes, self.diagonal_term[subd.volumes])
            subd.adm_pressure = self.adm_pressure[subd.volumes]
            local_flux_prescription = self.global_molar_flux_prescription[:, subd.volumes]

            subd.local_rhs = Gips.mount_independent_term(
                self.params['Vbulk'][subd.volumes],
                self.params['porosity'][subd.volumes],
                self.params['Cf'],
                self.params['dVtdP'][subd.volumes],
                self.params['P'][subd.volumes],
                len(subd.volumes),
                self.params['n_components'],
                self.params['n_phases'],
                subd.map_volumes[subd.adj_intern_local_faces],
                self.params['dVtdk'][:, subd.volumes],
                self.params['z_centroids'][subd.volumes],
                self.params['xkj_internal_faces'][:, :, self.map_internal_faces[subd.intern_local_faces]],
                self.params['Csi_j_internal_faces'][:, :, self.map_internal_faces[subd.intern_local_faces]],
                self.params['mobilities_internal_faces'][: , :, self.map_internal_faces[subd.intern_local_faces]],
                self.params['pretransmissibility_internal_faces'][self.map_internal_faces[subd.intern_local_faces]],
                self.params['Pcap'][:, subd.volumes],
                self.params['Vp'][subd.volumes],
                self.params['Vt'][subd.volumes],
                subd.map_volumes[subd.ind_neum],
                local_flux_prescription[:, subd.map_volumes[subd.ind_neum]],
                self.params['delta_t'],
                self.params['g'],
                subd.map_volumes[subd.ind_diric],
                subd.adm_pressure[subd.map_volumes[subd.ind_diric]],
                subd.map_volumes[subd.ind_diric],
                self.params['rho_j'][:, :, subd.volumes],
                self.params['rho_j_internal_faces'][:, :, self.map_internal_faces[subd.intern_local_faces]]
            )
            # local_rhs += subd.flux_prescription
            # local_rhs[subd.map_volumes[subd.ind_diric]] = adm_pressure[subd.ind_diric]
            # subd.local_rhs = local_rhs
            subd.Tlocal = set_matrix_pressure_prescription(
                subd.Tlocal_no_bc,
                subd.map_volumes[subd.ind_diric]
            )
            
            resp = solver.direct_solver(subd.Tlocal, subd.local_rhs)
            self.queue.put([subd.volumes, resp])

        self.finish()
        
        