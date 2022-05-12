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
    def richardson_solver(
        A: csc_matrix, 
        b: np.ndarray, 
        x0: np.ndarray, 
        OR: csc_matrix, 
        OP: csc_matrix, 
        res_tol: float=1e-10, 
        x_tol: float=1e-10, 
        max_it: int=np.inf,
        pcorr: np.ndarray=None,
        **kwargs) -> np.ndarray:
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
            pcorr(np.ndarray): correction functions

        Returns:
            x (np.ndarray): answer
        """
        wells_producer = kwargs.get('wells_producer')
        if wells_producer:
            pass
        else:
            wells_producer = np.array([])
        # ids = np.setdiff1d(np.arange(len(x0)), wells_producer)

        it_counter = 0
        eps = 1e100
        # initial step is to use OR = OP.T
        R = OP.transpose().copy()
        Ac_it = A*OP
        Ac_it = R*Ac_it
        # Tc_inv = inv(Ac_it)
        x = x0.copy()
        x0_in = x0.copy()
        if pcorr:
            assert pcorr.shape == b.shape
        else:
            pcorr = np.zeros_like(b)
        # res_f = b-A*x
        res_f = b-A*(x-pcorr)
        res_c = R*res_f

        while it_counter < max_it and eps > x_tol:
            # res_c[:] = spsolve(Ac_it, R*(b-A*x))
            res_c[:] = spsolve(Ac_it, res_c)
            res_f[:] = OP*res_c + pcorr
            x += res_f
            res_f[:], exitcode = cg(A, b-A*x, maxiter=20, x0=res_f, tol=res_tol)
            x += res_f
            eps = np.absolute(res_f).max()/np.absolute(x).max()
            print(f'eps: {eps}')
            x0_in[:] = x.copy()
            it_counter += 1
            res_f[:] = b-A*(x-pcorr)
            res_c[:] = R*res_f

        res_c[:] = spsolve(OR*A*OP, OR*(b-A*(x-pcorr)))
        res_f[:] = OP*res_c + pcorr
        x += res_f

        # x=spsolve(A,b)

        return x, eps, it_counter
