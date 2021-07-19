import multiprocessing as mp
import numpy as np
from multiprocessing import Queue, Value
from ctypes import c_bool
import queue
from packs.multiscale.neuman_local_problems.master_local_solver import CommonMasterMethods, run_thing
from packs.multiscale.preprocess.dual_domains import get_bool_updated_dual_subdomains, DualSubdomain
from collections.abc import Sequence
from packs.multiscale.operators.prolongation.AMS.paralell2.local_operator import LocalOperator

class MasterLocalOperator(CommonMasterMethods):
    
    def __init__(self, problems_list: Sequence[DualSubdomain], n_volumes, diagonal_term, T_global):
        n_cpu = self.get_n_cpu()
        self.n_cpu = n_cpu - 4
        self.n_volumes = n_volumes
        bool_updated_subdomains = get_bool_updated_dual_subdomains(problems_list)
        n_problems = bool_updated_subdomains.sum()
        subdomains_to_update = problems_list[bool_updated_subdomains]
        for dual in subdomains_to_update:
            # dual.update_local_source_term(source_term)
            dual.update_t_local(T_global, diagonal_term[dual.gids])
            dual.update_as(dual.Tlocal)
 
        n_problems_per_cpu = self.count_problems(n_problems, self.n_cpu)
        problems_per_cpu = self.get_problems_per_cpu(n_problems_per_cpu, problems_list[bool_updated_subdomains])
        self.m2w, self.w2m, self.procs, self.queue, self.finished = self.init_subproblems(problems_per_cpu)

    def init_subproblems(problems_per_cpu):
        
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
        m2w, w2m, _queue, values = self.initialize_params(n)
        procs = [mp.Process(target=run_thing, args=[LocalOperator(subdomains, _queue, comm, finished, id_process)]) for subdomains, comm, finished, id_process in zip(problems_per_cpu, w2m, values, range(n))]

        return m2w, w2m, procs, _queue, values
    
    def all_process_finished(self):
        return np.all([i.value for i in self.finished])
    
    def run(self):
        pass