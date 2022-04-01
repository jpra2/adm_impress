import multiprocessing as mp

from mpmath.libmp.libintmath import python_bitcount
from packs.solvers.solvers_scipy import solver_sp
import numpy as np
from multiprocessing import Queue, Value
from ctypes import c_bool
import queue
from packs.multiscale.neuman_local_problems.local_solver import LocalSolver1, run_thing, LocalSolver2
from packs.compositional.IMPEC.global_pressure_solver import GlobalIMPECPressureSolver as Gips

class CommonMasterMethods:
    
    def count_problems(self, n_problems: int, n_cpu: int):

        if n_problems < n_cpu:
            n_cpu = n_problems

        n_problems_per_cpu = n_problems // n_cpu
        resto = n_problems % n_cpu

        n_problems_per_cpu = np.repeat(n_problems_per_cpu, n_cpu)
        if resto != 0:
            n_problems_per_cpu[0:resto] += 1

        return n_problems_per_cpu

    
    def get_problems_per_cpu(self, n_problems_per_cpu, problems_list):
        
        problems_per_cpu = []
        summ = 0
        for i in n_problems_per_cpu:
            problems_per_cpu.append(problems_list[summ:summ+i])
            summ += i

        return problems_per_cpu
    

    def get_n_cpu(self):
        n_cpu = mp.cpu_count() - 2
        if n_cpu <= 0:
            n_cpu = 1
        
        return n_cpu
    
    
    def initialize_params(self, n_cpu):
        n = n_cpu
        master2worker = [mp.Pipe() for _ in range(n)]
        _queue = Queue()
        m2w, w2m = list(zip(*master2worker))
        values = [Value(c_bool, False) for _ in range(n)]
        
        return m2w, w2m, _queue, values
    
    def set_false_finished(self, finished_list):
        for val in finished_list:
            val.value = False
            
        
    def all_process_finished(self, finished_list):
        return np.all([i.value for i in finished_list])



class MasterLocalSolver(CommonMasterMethods):

    def __init__(self, problems_list, n_volumes):
        self.n_cpu = self.get_n_cpu() - 1
        if self.n_cpu <= 0:
            self.n_cpu = 1
        self.n_volumes = n_volumes
        self.problems_list = problems_list
        n_problems = len(problems_list)
        n_problems_per_cpu = self.count_problems(n_problems, self.n_cpu)
        self.problems_per_cpu = self.get_problems_per_cpu(n_problems_per_cpu, problems_list)
        self.m2w, self.w2m, self.queue, self.finished = self.initialize_params(len(self.problems_per_cpu))
        # self.procs, self.procs_args, self.procs_targets, self.procs_kwargs = self.init_subproblems(self.problems_per_cpu, self.w2m, self.finished, self.queue)


    def init_subproblems(self, problems_per_cpu, w2m, values, _queue):
        """

        @param problems_per_cpu: list of list local problems in cpu shape = (?, n_cpu)
        @return: m2w: master to worker list communicator
                 w2m: worker to master list communicator
                 procs: list of process
                 queue: queue shared between processes
        """
        n = len(problems_per_cpu)
        # master2worker = [mp.Pipe() for _ in range(n)]
        # _queue = Queue()
        # m2w, w2m = list(zip(*master2worker))
        # values = [Value(c_bool, False) for _ in range(n)]
        # m2w, w2m, _queue, values = self.initialize_params(n)
        procs = [mp.Process(target=run_thing, args=[LocalSolver1(subdomains, _queue, comm, finished, id_process)]) for subdomains, comm, finished, id_process in zip(problems_per_cpu, w2m, values, range(n))]
        process_args = [proc._args for proc in procs]
        process_targets = [proc._target for proc in procs]
        process_kwargs = [proc._kwargs for proc in procs]

        # return m2w, w2m, procs, _queue, values
        return procs, process_args, process_targets, process_kwargs

    def init_subproblems2(self, problems_per_cpu, w2m, values, _queue, local_params, **kwargs):
        """

        @param problems_per_cpu: list of list local problems in cpu shape = (?, n_cpu)
        @return: m2w: master to worker list communicator
                 w2m: worker to master list communicator
                 procs: list of process
                 queue: queue shared between processes
        """
        n = len(problems_per_cpu)
        # master2worker = [mp.Pipe() for _ in range(n)]
        # _queue = Queue()
        # m2w, w2m = list(zip(*master2worker))
        # values = [Value(c_bool, False) for _ in range(n)]
        # m2w, w2m, _queue, values = self.initialize_params(n)
        procs = [mp.Process(target=run_thing, args=[LocalSolver2(subdomains, _queue, comm, finished, id_process, local_params, **kwargs)]) for subdomains, comm, finished, id_process in zip(problems_per_cpu, w2m, values, range(n))]
        process_args = [proc._args for proc in procs]
        process_targets = [proc._target for proc in procs]
        process_kwargs = [proc._kwargs for proc in procs]

        # return m2w, w2m, procs, _queue, values
        return procs, process_args, process_targets, process_kwargs

    def run(self):
        
        solution = np.zeros(self.n_volumes)
        # procs = self.init_subproblems(self.problems_per_cpu, self.w2m, self.finished, self.queue)
        # procs = self.procs
        
        procs, process_args, process_targets, process_kwargs = self.init_subproblems(self.problems_per_cpu, self.w2m, self.finished, self.queue)
        
        # self.procs = procs

        for proc in procs:
            proc.start()

        while(not self.all_process_finished(self.finished)):
            try:
                resp = self.queue.get_nowait()
            except queue.Empty:
                # print('\nFila vazia\n')
                pass
            else:
                solution[resp[0]] = resp[1]

        # for i, proc in enumerate(self.procs):
        #     proc.join()
        
        # for i, proc in enumerate(self.procs):
        #     proc._popen = None
        #     proc._args = self.procs_args[i]
        #     proc._target = self.procs_targets[i]
        #     proc._kwargs = self.procs_kwargs[i]
        
        # import pdb; pdb.set_trace()

        while(not self.queue.empty()):
            resp = self.queue.get()
            solution[resp[0]] = resp[1]
        
        self.set_false_finished(self.finished)

        return solution

    def run_serial(self):
        solution = np.zeros(self.n_volumes)
        from packs.solvers.solvers_scipy.solver_sp import SolverSp
        solver = SolverSp()
        for subd in self.problems_list:
            resp = solver.direct_solver(subd.Tlocal, subd.local_rhs)
            solution[subd.volumes] = resp
        
        return solution
    
    def run2(self, neumann_subds, local_params, **kwargs):
        
        solution = np.zeros(self.n_volumes)
        
        procs, process_args, process_targets, process_kwargs = self.init_subproblems2(self.problems_per_cpu, self.w2m, self.finished, self.queue, local_params, **kwargs)
        
        for proc in procs:
            proc.start()

        while(not self.all_process_finished(self.finished)):
            try:
                resp = self.queue.get_nowait()
            except queue.Empty:
                # print('\nFila vazia\n')
                pass
            else:
                solution[resp[0]] = resp[1]

        while(not self.queue.empty()):
            resp = self.queue.get()
            solution[resp[0]] = resp[1]
        
        self.set_false_finished(self.finished)

        return solution