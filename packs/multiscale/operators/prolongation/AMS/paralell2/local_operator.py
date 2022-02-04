from packs.multiscale.neuman_local_problems.local_solver import GlobalLocalSolver
from packs.multiscale.preprocess.dual_domains import DualSubdomain
import numpy as np
import scipy.sparse as sp

class LocalOperator(GlobalLocalSolver):

    def run(self):
        dtype_op = [('op_lines', np.int), ('op_cols', np.int), ('op_data', np.float)]
        dtype_cf = [('gids', np.int), ('cf', np.float)]

        self.initialize()
        for dual in self.subdomains:
            dual: DualSubdomain
            if dual.test_update():
                # dual.update_lu_matrices()
                local_op = dual.ams_solver.get_OP_AMS_TPFA_by_AS(dual.As)
                # local_op = dual.ams_solver.get_OP_AMS_TPFA_by_AS_and_local_lu(dual.As, dual.lu_matrices)
                local_op = sp.find(local_op)
                local_op[0][:] = dual.gids[local_op[0]]
                local_op[1][:] = dual.rmap_lcid_cid[local_op[1]]
            else:
                local_op = np.array([[], [], []])

            if self.update_FC is True:
                local_pcorr = dual.ams_solver.get_pcorr2(dual.As, dual.local_source_term)
            else:
                local_pcorr = np.zeros(len(dual.local_source_term))
            # local_pcorr = dual.ams_solver.get_pcorr3(dual.As, dual.local_source_term, dual.lu_matrices)

            op_resp = np.zeros(len(local_op[0]), dtype=dtype_op)
            cf_resp = np.zeros(len(local_pcorr), dtype=dtype_cf)

            op_resp['op_lines'][:] = local_op[0]
            op_resp['op_cols'][:] = local_op[1]
            op_resp['op_data'][:] = local_op[2]
            cf_resp['gids'][:] = dual.gids
            cf_resp['cf'][:] = local_pcorr

            resp = np.array([
                op_resp,
                cf_resp
            ], dtype='O')


            self.queue.put(resp)

        self.finish()
