
def run_thing(local_solver_obj):
    local_solver_obj.run()


class GlobalLocalSolver:
    def __init__(self, subdomains, queue, comm, finished, id):
        """
        @param subdomains: list of subdomains
        @param queue: queue object
        @param comm: communicator
        """
        self.subdomains = subdomains
        self.queue = queue
        self.comm = comm
        self.finished = finished
        self.id = id

    def finish(self):
        print(f'\nProcess {self.id} finished \n')
        self.finished.value = True


class LocalSolver1(GlobalLocalSolver):

    def run(self):
        self.finish()
        pass
