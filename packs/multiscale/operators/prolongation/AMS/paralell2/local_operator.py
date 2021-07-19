from packs.multiscale.operators.prolongation.AMS.ams_tpfa import AMSTpfa
from packs.multiscale.neuman_local_problems.local_solver import GlobalLocalSolver
from packs.solvers.solvers_scipy.solver_sp import SolverSp
# from packs.multiscale.preprocess.dual_domains import DualSubdomain

class LocalOperator(GlobalLocalSolver):
    
    def run(self):
        self.initialize()
        for dual in self.subdomains:
            pass
            
        self.finish()
    
