import multiprocessing as mp
from packs.compositional.prep_FR import correction_function
import numpy as np
from multiprocessing import Queue, Value
from ctypes import c_bool
import queue
from packs.multiscale.neuman_local_problems.master_local_solver import CommonMasterMethods
from packs.multiscale.neuman_local_problems.local_solver import run_thing, run_thing2
from packs.multiscale.preprocess.dual_domains import DualSubdomainMethods, DualSubdomain
from collections.abc import Sequence
from packs.multiscale.operators.prolongation.AMS.paralell2.local_operator import LocalOperator
import scipy.sparse as sp

class MasterLocalOperator(CommonMasterMethods):

    def __init__(self, problems_list, n_volumes, update_FC=False, **kwargs):
        self.update_FC = update_FC
        problems_list: Sequence[DualSubdomain]
        self.n_cpu = self.get_n_cpu()
        # self.n_cpu = 1
        # self.n_cpu = n_cpu - 1
        self.n_volumes = n_volumes
        n_problems = len(problems_list)

        n_problems_per_cpu = self.count_problems(n_problems, self.n_cpu)
        self.problems_per_cpu = self.get_problems_per_cpu(n_problems_per_cpu, problems_list)
        # self.m2w, self.w2m, self.procs, self.queue, self.finished = self.init_subproblems(problems_per_cpu)
        self.m2w, self.w2m, self.queue, self.finished = self.initialize_params(len(self.problems_per_cpu))
        self.procs, self.procs_args, self.procs_targets, self.procs_kwargs = self.init_subproblems(self.problems_per_cpu, self.w2m, self.finished, self.queue, self.update_FC, **kwargs)

    def init_subproblems(self, problems_per_cpu, w2m, values, _queue, update_FC, **kwargs):

        """

        @param problems_per_cpu: list of list local problems in cpu shape = (?, n_cpu)
        @return: m2w: master to worker list communicator
                w2m: worker to master list communicator
                procs: list of process
                queue: queue shared between processes
        """
        n = len(problems_per_cpu)
        # m2w, w2m, _queue, values = self.initialize_params(n)
        procs = [mp.Process(target=run_thing, args=[LocalOperator(subdomains, _queue, comm, finished, id_process, update_FC=update_FC, **kwargs)]) for subdomains, comm, finished, id_process in zip(problems_per_cpu, w2m, values, range(n))]
        process_args = [proc._args for proc in procs]
        process_targets = [proc._target for proc in procs]
        process_kwargs = [proc._kwargs for proc in procs]

        # return m2w, w2m, procs, _queue, values
        return procs, process_args, process_targets, process_kwargs

    def init_subproblems2(self, problems_per_cpu, w2m, values, _queue, update_FC, **kwargs):
        #### run_thing2 for local subproblems

        """

        @param problems_per_cpu: list of list local problems in cpu shape = (?, n_cpu)
        @return: m2w: master to worker list communicator
                w2m: worker to master list communicator
                procs: list of process
                queue: queue shared between processes
        """
        n = len(problems_per_cpu)
        # m2w, w2m, _queue, values = self.initialize_params(n)
        procs = [mp.Process(target=run_thing2, args=[LocalOperator(subdomains, _queue, comm, finished, id_process, update_FC=update_FC, **kwargs)]) for subdomains, comm, finished, id_process in zip(problems_per_cpu, w2m, values, range(n))]
        process_args = [proc._args for proc in procs]
        process_targets = [proc._target for proc in procs]
        process_kwargs = [proc._kwargs for proc in procs]

        # return m2w, w2m, procs, _queue, values
        return procs, process_args, process_targets, process_kwargs

    def run(self, OP_AMS, dual_subdomains):

        correction_function = np.zeros(self.n_volumes)
        # procs = self.init_subproblems(self.problems_per_cpu, self.w2m, self.finished, self.queue)
        procs, procs_args, procs_targets, procs_kwargs = self.init_subproblems(self.problems_per_cpu, self.w2m, self.finished, self.queue, self.update_FC)

        for proc in procs:
            proc.start()

        while(not self.all_process_finished(self.finished)):
            try:
                resp = self.queue.get_nowait()
            except queue.Empty:
                # print('\nFila vazia\n')
                pass
            else:
                set_data_to_op(OP_AMS, resp[0])
                set_data_to_cf(correction_function, resp[1])

        # for i, proc in enumerate(self.procs):
        #     proc.join()
        #     proc._popen = None
        #     proc._args = self.procs_args[i]
        #     proc._target = self.procs_targets[i]
        #     proc._kwargs = self.procs_kwargs[i]

        while(not self.queue.empty()):
            resp = self.queue.get()
            set_data_to_op(OP_AMS, resp[0])
            set_data_to_cf(correction_function, resp[1])

        self.set_false_finished(self.finished)

        return OP_AMS, correction_function

    def run2(self, OP_AMS, dual_subdomains, **kwargs):

        correction_function = np.zeros(self.n_volumes)
        # procs = self.init_subproblems(self.problems_per_cpu, self.w2m, self.finished, self.queue)
        procs, procs_args, procs_targets, procs_kwargs = self.init_subproblems2(self.problems_per_cpu, self.w2m, self.finished, self.queue, self.update_FC, **kwargs)

        for proc in procs:
            proc.start()

        while(not self.all_process_finished(self.finished)):
            try:
                resp = self.queue.get_nowait()
            except queue.Empty:
                # print('\nFila vazia\n')
                pass
            else:
                set_data_to_op(OP_AMS, resp[0])
                set_data_to_cf(correction_function, resp[1])

        # for i, proc in enumerate(self.procs):
        #     proc.join()
        #     proc._popen = None
        #     proc._args = self.procs_args[i]
        #     proc._target = self.procs_targets[i]
        #     proc._kwargs = self.procs_kwargs[i]

        while(not self.queue.empty()):
            resp = self.queue.get()
            set_data_to_op(OP_AMS, resp[0])
            set_data_to_cf(correction_function, resp[1])

        for comm in self.m2w:
            as_list = comm.recv()
            for resp in as_list:
                dual_subdomains[resp[0]].As.update(resp[1])

        self.set_false_finished(self.finished)

        return OP_AMS, correction_function

    def run_serial(self, OP_AMS, dual_subdomains, **kwargs):

        if self.update_FC:
            return self.run_serial_update_FC(OP_AMS, dual_subdomains, **kwargs)
        else:
            correction_function = np.zeros(self.n_volumes)

            for dual in dual_subdomains:
                if dual.test_update():
                    # dual.update_lu_matrices()
                    local_op = sp.find(dual.ams_solver.get_OP_AMS_TPFA_by_AS(dual.As))
                    # local_op = dual.ams_solver.get_OP_AMS_TPFA_by_AS_and_local_lu(dual.As, dual.lu_matrices)
                    # local_op = sp.find(local_op)
                    local_op[0][:] = dual.gids[local_op[0]]
                    local_op[1][:] = dual.rmap_lcid_cid[local_op[1]]
                    OP_AMS[local_op[0], local_op[1]] = local_op[2]

            return OP_AMS, correction_function

    def run_serial_update_FC(self, OP_AMS, dual_subdomains, **kwargs):
        correction_function = np.zeros(self.n_volumes)

        for dual in dual_subdomains:
            if dual.test_update():
                # dual.update_lu_matrices()
                local_op = sp.find(dual.ams_solver.get_OP_AMS_TPFA_by_AS(dual.As))
                # local_op = dual.ams_solver.get_OP_AMS_TPFA_by_AS_and_local_lu(dual.As, dual.lu_matrices)
                # local_op = sp.find(local_op)
                local_op[0][:] = dual.gids[local_op[0]]
                local_op[1][:] = dual.rmap_lcid_cid[local_op[1]]
                OP_AMS[local_op[0], local_op[1]] = local_op[2]

                local_pcorr = dual.ams_solver.get_pcorr2(dual.As, dual.local_source_term)
                correction_function[dual.gids] = local_pcorr

        return OP_AMS, correction_function

def set_data_to_op_dep0(OP_AMS, resp):

    print('setting_data')

    if len(resp['op_lines']) == 0:
        pass
    else:
        OP_AMS[resp['op_lines'], resp['op_cols']] = resp['op_data']

def set_data_to_cf_dep0(resp, correction_function):

    print('setting cf')
    print(resp)
    correction_function[resp['gids']] = resp['cf']

def set_data_to_op(OP_AMS, op_resp):
    if len(op_resp['op_lines']) == 0:
        pass
    else:
        OP_AMS[op_resp['op_lines'], op_resp['op_cols']] = op_resp['op_data']

def set_data_to_cf(correction_function, cf_resp):
    correction_function[cf_resp['gids']] = cf_resp['cf']

def set_local_lu_matrices(resp_lu_matrices, dual_subdomains):
    if len(resp_lu_matrices) > 0:
        dual_subdomains[resp_lu_matrices[0]].lu_matrices.update(resp_lu_matrices[1])