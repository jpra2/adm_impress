from packs.multiscale.neuman_local_problems.local_solver import GlobalLocalSolver
from packs.multiscale.preprocess.dual_domains import DualSubdomain

class LocalOperator(GlobalLocalSolver):
    
    def run(self):
        self.initialize()
        for dual in self.subdomains:
            dual: DualSubdomain
            local_op = dual.ams_solver.get_OP_AMS_TPFA_by_AS(dual.As)
            local_pcorr = dual.ams_solver.get_pcorr2(dual.As, dual.local_source_term)
            
            
            
            pass
            
        self.finish()
    
