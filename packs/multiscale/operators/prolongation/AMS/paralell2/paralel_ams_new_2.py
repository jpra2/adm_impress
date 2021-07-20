import multiprocessing as mp
from packs.compositional.prep_FR import correction_function
import numpy as np
from multiprocessing import Queue, Value
from ctypes import c_bool
import queue
from packs.multiscale.neuman_local_problems.master_local_solver import CommonMasterMethods, run_thing
from packs.multiscale.preprocess.dual_domains import DualSubdomainMethods, DualSubdomain
from collections.abc import Sequence
from packs.multiscale.operators.prolongation.AMS.paralell2.local_operator import LocalOperator

class MasterLocalOperator(CommonMasterMethods):
    
    def __init__(self, problems_list: Sequence[DualSubdomain], n_volumes, diagonal_term, T_global):
        n_cpu = self.get_n_cpu()
        self.n_cpu = n_cpu - 4
        self.n_volumes = n_volumes
        n_problems = len(problems_list)
        DualSubdomainMethods.update_matrices_dual_subdomains(problems_list, T_global, diagonal_term)
 
        n_problems_per_cpu = self.count_problems(n_problems, self.n_cpu)
        problems_per_cpu = self.get_problems_per_cpu(n_problems_per_cpu, problems_list)
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
        m2w, w2m, _queue, values = self.initialize_params(n)
        procs = [mp.Process(target=run_thing, args=[LocalOperator(subdomains, _queue, comm, finished, id_process)]) for subdomains, comm, finished, id_process in zip(problems_per_cpu, w2m, values, range(n))]

        return m2w, w2m, procs, _queue, values
    
    def all_process_finished(self):
        return np.all([i.value for i in self.finished])
    
    def run(self, OP_AMS, n_volumes):
        
        correction_function = np.zeros(n_volumes)
        
        for proc in self.procs:
            proc.start()

        while(not self.all_process_finished()):
            try:
                resp = self.queue.get_nowait()
            except queue.Empty:
                # print('\nFila vazia\n')
                pass
            else:
                set_data_to_op(OP_AMS, resp)
                set_data_to_cf(correction_function, resp)

        for proc in self.procs:
            proc.join()

        while(not self.queue.empty()):
            resp = self.queue.get()
            set_data_to_op(OP_AMS, resp)
            set_data_to_cf(correction_function, resp)

        return OP_AMS, correction_function
        


def set_data_to_op(OP_AMS, resp):
    
    if len(resp[0]) == 0:
        pass
    else:
        OP_AMS[resp['op_lines'], resp['op_cols']] = resp['op_data']

def set_data_to_cf(resp, correction_function):
    correction_function[resp['gids']] = resp['cf']