import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import gmres, cg, bicgstab, spilu, splu, LinearOperator, spsolve
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
        

def iterative_ms_ilu0_bicgstab(
    A: sp.csc_matrix,
    b: np.ndarray,
    OP: sp.csc_matrix,
    OR: sp.csc_matrix,
    epsilon=1e-13,
    maxit=100
):
       
    rn = b.copy()
    rn2 = rn.copy()
    pn = b.copy()
    dp1 = b.copy()
    dp2 = b.copy()
    
    R = OR
    # R_ADM = OR_ADM
    
    ilu0 = spilu(A)
    Mx = lambda x: ilu0.solve(x)
    M = LinearOperator(A.shape, Mx)
    
    LU = splu((R*(A*OP)).tocsc())
    
    pn[:] = OP*LU.solve(R*b)
    rn[:] = b - A*pn
    
    err = 1e5
    
    it = 1

    # ##############################
    # ### metodo iterativo artur
    # while err > epsilon and it < maxit:
    #     dp1[:], exitcode = bicgstab(A, rn, M=M, tol=epsilon, maxiter=1)
    #     rn2[:] = rn - A*dp1
    #     # rc[:] = LU.solve(rn2)
    #     # dp2[:] = OP*rc
    #     dp2[:] = OP*LU.solve(R*rn2)
        
    #     pn[:] += dp1 + dp2
    #     rn[:] = b - A*pn
    #     err = np.linalg.norm(rn)
    #     it += 1
        
    #     # if err2 > err:
    #     #     import pdb; pdb.set_trace()
    #     # else:
    #     #     err = err2
    
    # dp1[:], exitcode = bicgstab(A, rn, M=M, tol=epsilon, maxiter=1)
    # pn[:] += dp1
    # rn[:] = b - A*pn
    # err = np.linalg.norm(rn)
    # #########################################
    
    ########################################
    ## metodo iterativo Filipe
    while err > epsilon and it < maxit:
        dp1[:] = OP*LU.solve(R*rn)
        rn2[:] = rn - A*dp1
        dp2[:], exitcode = bicgstab(A, rn2, M=M, tol=epsilon)
        pn[:] += dp1 + dp2
        rn[:] = b - A*pn
        err = np.linalg.norm(rn)
        # print(f'err: {err} \n')
        it += 1
    #######################################
    
    return pn, it, err
    
        
        
    
    
    
    






    