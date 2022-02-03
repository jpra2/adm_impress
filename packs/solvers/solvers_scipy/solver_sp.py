import scipy.sparse as sp
from scipy.sparse import linalg
import numpy as np


norm_rk = []
niter = 0

class ScipyCounter(object):
    def __init__(self, disp=True):
        self._disp = disp
        global norm_rk
        global niter
        norm_rk = []
        niter = 0
    def __call__(self, rk=None):
        global norm_rk
        global niter
        niter += 1
        norm_rk.append(np.linalg.norm(rk))
        # if self._disp:
        #     print('iter %3i\trk = %s' % (self.niter, str(rk)))


class SolverSp:

    def __init__(self):
        pass

    def direct_solver(self, A, b):

        # print('\nSolving direct solver spsolve\n')

        # A2 = A.tocsc().copy()

        solution = linalg.spsolve(A.tocsc(),b)

        return solution

    def lu_solver(self, A, b):

        print('\nSolving direct solver lu_solver\n')

        A2 = A.tocsc().copy()

        LU = linalg.splu(A2)
        solution = LU.solve(b)

        return solution

    def gmres_solver(self, A, b, x0=None, tol=1e-5, precond=None):

        print('\nSolving gmres solver\n')
        counter_callback = ScipyCounter(disp=False)

        n = A.shape[0]
        if precond:
            # M1 = linalg.spilu(A)
            # M_x = lambda x: M1.solve(x)

            M_x = lambda x: linalg.spsolve(A, x)
            M = linalg.LinearOperator((n, n), M_x)
        else:
            M = None

        # x, exitcode = linalg.gmres(A, b, x0=x0, tol=tol, M=M)
        x, exitcode = linalg.gmres(A, b, x0=x0, tol=tol, M=M, callback=counter_callback)
        ## exitcode = 0: indicates successful convergence

        return x

    def conjugate_gradient_solver(self, A, b, x0=None, tol=1e-5, precond=None):

        print('\nSolving conjugate gradient solver\n')

        n = A.shape[0]
        if precond:
            # M1 = linalg.spilu(A)
            # M_x = lambda x: M1.solve(x)

            M_x = lambda x: linalg.spsolve(A, x)
            M = linalg.LinearOperator((n, n), M_x)
        else:
            M = None

        x, exitcode = linalg.cg(A, b, x0=x0, tol=tol, M=M)

        return x
