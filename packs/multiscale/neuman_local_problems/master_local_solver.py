import multiprocessing as mp
import numpy as np
from multiprocessing import Queue, Value
from ctypes import c_bool
from packs.multiscale.neuman_local_problems.local_solver import LocalSolver1, run_thing
import queue


class MasterLocalSolver:

    def __init__(self, problems_list, n_volumes):
        n_cpu = mp.cpu_count()
        self.n_cpu = n_cpu - 4
        self.n_volumes = n_volumes
        n_problems = len(problems_list)
        n_problems_per_cpu = MasterLocalSolver.count_problems(n_problems, self.n_cpu)
        problems_per_cpu = MasterLocalSolver.get_problems_per_cpu(n_problems_per_cpu, problems_list)
        self.m2w, self.w2m, self.procs, self.queue, self.finished = MasterLocalSolver.init_subproblems(problems_per_cpu)

    @staticmethod
    def count_problems(n_problems: int, n_cpu: int):

        if n_problems < n_cpu:
            n_cpu = n_problems

        n_problems_per_cpu = n_problems // n_cpu
        resto = n_problems % n_cpu

        n_problems_per_cpu = np.repeat(n_problems_per_cpu, n_cpu)
        if resto != 0:
            n_problems_per_cpu[0:resto] += 1

        return n_problems_per_cpu

    @staticmethod
    def get_problems_per_cpu(n_problems_per_cpu, problems_list):
        problems_per_cpu = []
        summ = 0
        for i in n_problems_per_cpu:
            problems_per_cpu.append(problems_list[summ:summ+i])
            summ += i

        return problems_per_cpu

    @staticmethod
    def init_subproblems(problems_per_cpu):
        """

        @param problems_per_cpu: list of list local problems in cpu shape = (?, n_cpu)
        @return: m2w: master to worker list communicator
                 w2m: worker to master list communicator
                 procs: list of process
                 queue: queue shared between processes
        """
        n = len(problems_per_cpu)
        master2worker = [mp.Pipe() for _ in range(n)]
        queue = Queue()
        m2w, w2m = list(zip(*master2worker))
        values = [Value(c_bool, False) for _ in range(n)]
        procs = [mp.Process(target=run_thing, args=[LocalSolver1(subdomains, queue, comm, finished, id_process)]) for subdomains, comm, finished, id_process in zip(problems_per_cpu, w2m, values, range(n))]

        return m2w, w2m, procs, queue, values

    def run(self):
        solution = np.zeros(self.n_volumes)

        for proc in self.procs:
            proc.start()

        while(not self.all_process_finished()):
            try:
                resp = self.queue.get_nowait()
            except queue.Empty:
                # print('\nFila vazia\n')
                pass
            else:
                solution[resp[0]] = resp[1]

        for proc in self.procs:
            proc.join()

        while(not self.queue.empty()):
            resp = self.queue.get()
            solution[resp[0]] = resp[1]

        return solution


    def all_process_finished(self):
        return np.all([i.value for i in self.finished])



