from operator import inv
import pdb
from scipy.sparse import csc_matrix
import numpy as np
from scipy.sparse.linalg import gmres, lgmres, spsolve, cg, inv

class TamsSolverFV:
    """
        Finite Volume TAMS solver
        Hui Zhou and Hamdi A. Tchelepi
        Two-Stage Algebraic Multiscale Linear Solver for Highly Heterogeneous Reservoir Models
        doi - https://doi.org/10.2118/141473-pa
    """
    
    @staticmethod
    def richardson_solver(A: csc_matrix, b: np.ndarray, x0: np.ndarray, OR: csc_matrix, OP: csc_matrix, res_tol: float = 1e-10, x_tol: float = 1e-10, max_it: int = np.inf, **kwargs) -> np.ndarray:
        """Richardson iterative solver

        Args:
            A (csc): left side matrix
            b (np.ndarray): rhs term
            x0 (np.ndarray): initial guess
            OR (csc): restriction operator
            OP (csc): prolongation operator
            res_tol (float): tolerance for L2 residual error - np.linalg.norm(rk) 
            x_tol (float): L_inf error tolerance = np.absolute(x).max()
            max_it (int): max iteration counter, if max_it = np.inf (default), the code runs until b_tol and res_tol. 

        Returns:
            x (np.ndarray): answer
        """
        
        it_counter = 0
        eps = 1e100
        # initial step is to use OR = OP.T
        R = OP.transpose().copy()
        Ac_it = A*OP
        Ac_it = R*Ac_it
        # Tc_inv = inv(Ac_it)
        x = x0.copy()
        x0_in = x0.copy()
        res_f = b-A*x0
        res_c = OR*res_f
        
        while it_counter < max_it and eps > x_tol:
            # x += OP*spsolve(Ac_it, R*(b-A*x))
            # resp, exitcode = gmres(A, b-A*x, x0=np.zeros(A.shape[0]), tol=res_tol)
            # resp, exitcode = gmres(A, b-A*x, x0=np.zeros(A.shape[0]), tol=res_tol)
            # x += resp
            # x += spsolve(A, b-A*x)
            # res = OP*Tc_inv*OR*(b - A*x) 
            # res_c[:], exitcode = gmres(Ac_it, R*(b-A*x), x0=res_c, maxiter=2, tol=res_tol)
            res_c[:] = spsolve(Ac_it, R*(b-A*x))
            res_f[:] = OP*res_c
            x += res_f
            res_f[:], exitcode = cg(A, b-A*x, maxiter=5, x0=res_f ,tol=res_tol)
            x += res_f
            # eps = np.absolute((x - x0_in)/x).max()
            eps = np.absolute(res_f).max()/np.absolute(x).max()
            print(f'eps: {eps}')
            x0_in[:] = x.copy()
            it_counter += 1
        
        # x += OP*spsolve(OR*A*OP, OR*(b-A*x))
        # res_c[:], exitcode = gmres(OR*A*OP, OR*(b-A*x), x0=res_c, maxiter=2, tol=res_tol)
        res_c[:] = spsolve(OR*A*OP, OR*(b-A*x))
        res_f[:] = OP*res_c
        x += res_f
        
        return x