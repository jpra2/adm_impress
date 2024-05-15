import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import gmres, cg, bicgstab, spilu, splu, LinearOperator
from packs.solvers.solvers_scipy.solver_sp import SolverSp

def ms_solve_it(
        A: sp.csc_matrix, 
        OP: sp.csc_matrix, 
        OR: sp.csc_matrix,
        b: np.ndarray,
        x0: np.ndarray,
        epsilon=0.001,
        maxiter = 10,
        internal_loop_epsilon=0.001,
        internal_loop_maxiter=10
    ):

    '''
    multiscale iterative solver  
    '''

    R = OP.transpose(s)
    LU = splu(R*(A*OP))

    pn = x0.copy()
    rn = b - A*pn


    # B = A == A.transpose()
    # import pdb; pdb.set_trace()

    scipy_solver = SolverSp()

    x = x0.copy()

    # ilu0 = spilu(A)
    # Mx = lambda x: ilu0.solve(x)
    # M = LinearOperator(A.shape, Mx)

    # res_f = b-A*x
    # res_c = R*res_f

    rn_norm = 1e5
    count = 1
    while count <= maxiter and rn_norm > epsilon:
        # solution1[:], exit_code = splinalg.gmres(Ac, OR*rn, tol=internal_loop_epsilon, maxiter=internal_loop_maxiter)
        # res_c[:] = LU.solve(R*(b-A*x))
        # res_f[:] = OP*res_c
        # x += res_f
        # # res_f[:] = scipy_solver.LinearCG(A, b-A*x, x0=res_f, tol=internal_loop_epsilon, maxiter=100)
        # res_f[:] = scipy_solver.gmres_solver(A, b-A*x, x0=res_f, tol=internal_loop_epsilon, maxiter=100)
        # x += res_f
        
        delta_pms = OP*(LU.solve(OR*rn))
        rn2 = rn - A*delta_pms
        delta_p, exitcode = gmres(A, rn2, x0=delta_pms, tol=1e-5, maxiter=internal_loop_maxiter)
        pn += delta_pms + delta_p
        rn = b - A*pn
        x = pn

        # solution1[:] = LU.solve(OR*rn)
        # delta_p[:] = OP*solution1
        # rn[:] = rn - A*delta_p
        # pn[:] = pn + delta_p
        # delta_p[:], exit_code = gmres(A, rn, tol=internal_loop_epsilon, maxiter=internal_loop_maxiter, x0=delta_p[:])
        # pn[:] = pn + delta_p
        # rn[:] = b - A*pn

        rn_norm = np.linalg.norm(b-A*x)
        count += 1
        print(f'Count: {count}')
        print(f'Norm: {rn_norm} \n')
    
    print(f'Solved with {count} iterations \n')
    x = pn
    return x
        








    