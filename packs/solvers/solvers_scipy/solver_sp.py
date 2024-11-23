import scipy.sparse as sp
from scipy.sparse import linalg
import numpy as np


norm_rk = []
niter = 0

class ScipyCounter(object):
    def __init__(self, disp=False):
        self._disp = disp
        self.norm_rk = []
        self.niter = 0
    def __call__(self, rk=None):
        
        self.niter += 1
        self.norm_rk.append(np.linalg.norm(rk))
        if self._disp is True:
            print(self.niter)


class SolverSp:

    def __init__(self):
        pass

    def direct_solver(self, A, b):

        # print('\nSolving direct solver spsolve\n')

        # A2 = A.tocsc().copy()

        solution = linalg.spsolve(A.tocsc(),b)

        return solution

    def lu_solver(self, A, b):

        LU = linalg.splu(A)
        solution = LU.solve(b)

        return solution

    def gmres_solver(self, A, b, x0=None, tol=1e-5, M=None, maxiter=None):

        counter_callback = ScipyCounter(disp=False)

        # x, exitcode = linalg.gmres(A, b, x0=x0, tol=tol, M=M)
        x, exitcode = linalg.gmres(A, b, x0=x0, tol=tol, M=M, callback=counter_callback, maxiter=maxiter)
        ## exitcode = 0: indicates successful convergence

        return x

    def conjugate_gradient_solver(self, A, b, x0=None, tol=1e-5, M=None, maxiter=None):

        
        counter_callback = ScipyCounter(disp=False)

        x, exitcode = linalg.cg(A, b, x0=x0, tol=tol, M=M, maxiter=maxiter, callback=counter_callback)

        return x

    def LinearCG(self, A, b, x0, tol=1e-5, maxiter=100):
        xk = x0.copy()
        rk = A*xk - b
        pk = -rk
        rk_norm = np.linalg.norm(rk)
        
        # print()
        # print('Linear CG solve...')
        # print()
        
        num_iter = 0
        # curve_x = [xk]
        curve_x = xk.copy()
        while rk_norm > tol and num_iter < maxiter:
            apk = A*pk
            rkrk = np.dot(rk, rk)
            
            alpha = rkrk / np.dot(pk, apk)
            xk = xk + alpha * pk
            rk = rk + alpha * apk
            beta = np.dot(rk, rk) / rkrk
            pk = -rk + beta * pk
            
            num_iter += 1
            # curve_x.append(xk)
            curve_x = xk.copy()
            rk_norm = np.linalg.norm(rk)
            # print('Iteration: {} \t x = {} \t residual = {:.4f}'.format(num_iter, xk, rk_norm))
            # print('Iteration: {} \t x = {} \t residual = {}'.format(num_iter, xk, rk_norm))
        
        # print('\nSolution: \t x = {}'.format(xk))
            
        # return np.array(curve_x[-1])
        return np.array(curve_x)
    
    def get_spilu_precond(self, A, fill_factor=None):
        B = linalg.spilu(A, fill_factor=fill_factor)
        Mx = lambda x: B.solve(x)
        M = linalg.LinearOperator(A.shape, Mx)
        return M
    
    def bicgstab(self, A, b, x0=None, tol=1e-5, M=None, maxiter=None):
        counter_callback = ScipyCounter(disp=False)

        # x, exitcode = linalg.gmres(A, b, x0=x0, tol=tol, M=M)
        x, exitcode = linalg.bicgstab(A, b, x0=x0, tol=tol, M=M, callback=counter_callback, maxiter=maxiter)
        ## exitcode = 0: indicates successful convergence

        return x