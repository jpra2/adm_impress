import pyamg
import scipy.sparse as sp

class SolverAMG:

    def __init__(self):
        pass

    def smoothed_aggregation_solver(self, A: sp.csc_matrix, b, TOL=1e-10):
        
        ml = pyamg.smoothed_aggregation_solver(A)
        solution = ml.solve(b, tol=TOL)
        return solution
    
    def multigrid_solver(self, A: sp.csc_matrix, b, TOL=1e-10, max_coarse=10, max_levels=10, x0=None, accel=None, residuals=None):
        """
        Multigrid solver
        accel_list = ['gmres', 'cg']
        """
        ml = pyamg.smoothed_aggregation_solver(A, max_coarse=max_coarse, max_levels=max_levels)
        # residuals = []
        solution = ml.solve(b, tol=TOL, accel=accel, residuals=residuals, x0=x0)
        return solution, residuals
    
    def multigrid_as_precond(A, max_coarse=10, max_levels=10):
        ml = pyamg.smoothed_aggregation_solver(A, max_coarse=max_coarse, max_levels=max_levels)
        M = ml.aspreconditioner()
        return M


