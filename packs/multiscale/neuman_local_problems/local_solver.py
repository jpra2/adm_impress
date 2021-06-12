
def run_thing(local_solver_obj):
    local_solver_obj.run()


class GlobalLocalSolver:
    def __init__(self, subdomains, queue, comm, finished, id_process):
        """
        @param subdomains: list of subdomains
        @param queue: queue object
        @param comm: communicator
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
        self.finish()
