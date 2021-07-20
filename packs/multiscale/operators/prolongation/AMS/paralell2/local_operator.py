from packs.multiscale.neuman_local_problems.local_solver import GlobalLocalSolver
from packs.multiscale.preprocess.dual_domains import DualSubdomain
import numpy as np
import scipy.sparse as sp

class LocalOperator(GlobalLocalSolver):
    
    def run(self):
        self.initialize()
        for dual in self.subdomains:
            dual: DualSubdomain
            if dual.test_update():
                local_op = dual.ams_solver.get_OP_AMS_TPFA_by_AS(dual.As)
                local_op = sp.find(local_op)
                local_op[1][:] = dual.rmap_lcid_cid[local_op[1]]
            else:
                local_op = np.array([[], [], []])
            
            local_pcorr = dual.ams_solver.get_pcorr2(dual.As, dual.local_source_term)
        
            resp = {
                'op_lines': local_op[0],
                'op_cols':  local_op[1],
                'op_data':  local_op[2],
                'gids':     dual.gids,
                'cf':       local_pcorr
            }
            
            dual.reinitialize_local_update()
            self.queue.put(resp)
            
        self.finish()
    
