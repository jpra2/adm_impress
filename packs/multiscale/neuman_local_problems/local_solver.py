from packs.solvers.solvers_scipy.solver_sp import SolverSp
from multiprocessing import Queue

def run_thing(local_solver_obj):
    local_solver_obj.run()


class GlobalLocalSolver:
    def __init__(self, subdomains, queue: Queue, comm, finished, id_process):
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
