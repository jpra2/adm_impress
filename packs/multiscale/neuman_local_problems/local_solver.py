from packs.solvers.solvers_scipy.solver_sp import SolverSp
from multiprocessing import Queue
import os

def run_thing(local_solver_obj):
    local_solver_obj.run()

def run_thing2(local_solver_obj):
    local_solver_obj.run2()


class GlobalLocalSolver:
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

    def finish(self):
        print(f'\nProcess {self.id_process} finished \n')
        self.finished.value = True

    def initialize(self):
        print(f'\nProcess {self.id_process} start \n')


class LocalSolver1(GlobalLocalSolver):

    def run(self):
        self.initialize()
        solver = SolverSp()
        for subd in self.subdomains:
            resp = solver.direct_solver(subd.Tlocal, subd.local_rhs)
            self.queue.put([subd.volumes, resp])

        self.finish()
