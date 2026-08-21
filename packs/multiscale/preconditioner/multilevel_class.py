import pyamg
import numpy as np
from pyvista import Line
from scipy.sparse.linalg import LinearOperator, factorized, spsolve_triangular
import scipy.sparse as sp
import ilupp
import pyamg
from pyamg.multilevel import multilevel_solver, MultilevelSolver
import copy

from packs.multiscale.transmissibility_correction.algorithimic_monotone import AlgorithimicMonotone


class Ilu0Precond(LinearOperator):
    def __init__(self, A: sp.csc_matrix, **kwargs):
        ilu = ilupp.ILU0Preconditioner(A)
        self.L = copy.deepcopy(ilu.factors()[0])
        self.U = copy.deepcopy(ilu.factors()[1])

    def aapply(self, v: np.ndarray):
        y = spsolve_triangular(self.L, v, lower=True)
        x = spsolve_triangular(self.U, y, lower=False)
        return x

    def _matvec(self, v: np.ndarray):
        y = spsolve_triangular(self.L, v, lower=True)
        x = spsolve_triangular(self.U, y, lower=False)
        return x
    
    def get_copy(self):
        return copy.deepcopy(self)
    
    @property
    def shape(self):
        return self.L.shape
    
    @property
    def dtype(self):
        return self.L.dtype

class Ilu1Precond(Ilu0Precond):
    def __init__(self, A: sp.csc_matrix, threshold: float=1e-16, **kwargs):
        nnz_por_linha = A.nnz // A.shape[0]
        ilu = ilupp.ILUTPreconditioner(A, fill_in=int(nnz_por_linha*1.5), threshold=threshold)
        # ilu = ilupp.ILUTPreconditioner(A, fill_in=2, threshold=threshold)
        self.L = ilu.factors()[0]
        self.U = ilu.factors()[1]

class MultiScaleIlu0Smoother2OP(LinearOperator):
    def __init__(self, A0, P1, R1, P2, R2, npre=0, npost=1, **kwargs):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        A1_1 = R1 @ A0 @ P1
        A1_2 = R2 @ A0 @ P2

        # 2. Solvers Diretos para as Malhas Grossas
        # factorized() é um wrapper rápido para o splu do scipy
        self.solve_grossa1 = factorized(A1_1.tocsc())
        self.solve_grossa2 = factorized(A1_2.tocsc())

        self.P1, self.R1 = P1, R1
        self.P2, self.R2 = P2, R2
        self.x = np.zeros(A0.shape[0])
        self.r = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])
        self.d2 = np.zeros(R2.shape[0])
        self.fs_smoother = self.create_fs_smoother(A0, **kwargs)
        
    
    def create_fs_smoother(self, A0, **kwargs) -> LinearOperator:
        return  Ilu0Precond(A0, **kwargs) 

    def _matvec_v0(self, b):
        
        self.x[:] = self.fs_smoother.matvec(b)

        # r1 = b - self.A0 @ x
        # d1 = self.R1 @ r1
        # c1 = self.solve_grossa1(d1)
        # x = x + (self.P1 @ c1)
        self.r[:] = b - self.A0 @ self.x
        self.d1[:] = self.R1 @ self.r
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + (self.P1 @ self.d1)

        # r2 = b - self.A0 @ x
        # d2 = self.R2 @ r2
        # c2 = self.solve_grossa2(d2)
        # x[:] = x + (self.P2 @ c2)
        self.r[:] = b - self.A0 @ self.x
        self.d2[:] = self.R2 @ self.r
        self.d2[:] = self.solve_grossa1(self.d2)
        self.x[:] = self.x + (self.P2 @ self.d2)

        return self.x.copy()
    
    def _matvec_v1(self, b):
        
        self.d1[:] = self.R1 @ b
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.P1 @ self.d1
        
        self.r[:] = b - self.A0 @ self.x
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        
        self.r[:] = b - self.A0 @ self.x
        self.d2[:] = self.R2 @ self.r
        self.d2[:] = self.solve_grossa2(self.d2)
        self.x[:] = self.x + (self.P2 @ self.d2)
        
        self.r[:] = b - self.A0 @ self.x
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)

        return self.x.copy()
    
    def _matvec_v2(self, b):
        self.x[:] = 0
        self.r[:] = b
        
        # for _ in range(self.npre):
        #     self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        #     self.r[:] = b - self.A0 @ self.x
        
        self.d1[:]  = self.R1 @ self.r
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + self.P1 @ self.d1
        self.r[:] = b - self.A0 @ self.x    
    
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        self.r[:] = b - self.A0 @ self.x
        
        # for _ in range(self.npost):
        #     self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        #     self.r[:] = b - self.A0 @ self.x
        
        self.d2[:] = self.R2 @ self.r
        self.d2[:] = self.solve_grossa2(self.d2)
        self.x[:] = self.x + self.P2 @ self.d2
        self.r[:] = b - self.A0 @ self.x
        
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        self.r[:] = b - self.A0 @ self.x
        
        # for _ in range(self.npost):
        #     self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        #     self.r[:] = b - self.A0 @ self.x
        
        return self.x.copy()
        
        
        
    
    def _matvec(self, b):
        # return self._matvec_v1(b)
        return self._matvec_v2(b)


class MultiScaleIlu1Smoother2OP(MultiScaleIlu0Smoother2OP):
    def __init__(self, A0, P1, R1, P2, R2, npre=0, npost=1, threshold: float=1e-16, ilu_factor: float=1.5, **kwargs):
        
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        A1_1 = R1 @ A0 @ P1
        A1_2 = R2 @ A0 @ P2

        # 2. Solvers Diretos para as Malhas Grossas
        # factorized() é um wrapper rápido para o splu do scipy
        self.solve_grossa1 = factorized(A1_1.tocsc())
        self.solve_grossa2 = factorized(A1_2.tocsc())

        self.P1, self.R1 = P1, R1
        self.P2, self.R2 = P2, R2
        self.x = np.zeros(A0.shape[0])
        self.r = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])
        self.d2 = np.zeros(R2.shape[0])
        self.fs_smoother = self.create_fs_smoother(A0, threshold=threshold, ilu_factor=ilu_factor, **kwargs)
        
    def create_fs_smoother(self, A0, threshold: float=1e-16, ilu_factor: float=1.5, **kwargs) -> LinearOperator:
        return Ilu1Precond(A0, threshold=threshold)


class MultiScaleIlu0Smoother2OP_new1(LinearOperator):
    def __init__(self, A0, P1, R1, P2, R2, x0, npre=0, npos=1):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npos = npos

        self.fs_smoother = Ilu0Precond(A0)

        A1_1 = R1 @ A0 @ P1
        A1_2 = R2 @ A0 @ P2

        # 2. Solvers Diretos para as Malhas Grossas
        # factorized() é um wrapper rápido para o splu do scipy
        self.solve_grossa1 = factorized(A1_1.tocsc())
        self.solve_grossa2 = factorized(A1_2.tocsc())

        self.P1, self.R1 = P1, R1
        self.P2, self.R2 = P2, R2
        if x0:
            self.x = x0.copy()
        else:
            self.x = np.zeros(A0.shape[0])
        self.r = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])
        self.d2 = np.zeros(R2.shape[0])

    def _matvec_v0(self, b):
        
        self.x[:] = self.fs_smoother.apply(b)

        # r1 = b - self.A0 @ x
        # d1 = self.R1 @ r1
        # c1 = self.solve_grossa1(d1)
        # x = x + (self.P1 @ c1)
        self.r[:] = b - self.A0 @ self.x
        self.d1[:] = self.R1 @ self.r
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + (self.P1 @ self.d1)

        # r2 = b - self.A0 @ x
        # d2 = self.R2 @ r2
        # c2 = self.solve_grossa2(d2)
        # x[:] = x + (self.P2 @ c2)
        self.r[:] = b - self.A0 @ self.x
        self.d2[:] = self.R2 @ self.r
        self.d2[:] = self.solve_grossa1(self.d2)
        self.x[:] = self.x + (self.P2 @ self.d2)

        return self.x.copy()
    
    def _matvec_v1(self, b):
        
        for _ in range(self.npre):
            self.x[:] = self.x + self.fs_smoother.matvec(b - self.A0 @ self.x)
        
        self.r[:] = b - self.A0 @ self.x
        self.d1[:] = self.R1 @ (self.r)
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + self.P1 @ self.d1
        
        for _ in range(self.npos):
            self.x[:] = self.x + self.fs_smoother.matvec(b - self.A0 @ self.x)
        
        self.r[:] = b - self.A0 @ self.x
        self.d2[:] = self.R2 @ self.r
        self.d2[:] = self.solve_grossa2(self.d2)
        self.x[:] = self.x + (self.P2 @ self.d2)
        
        for _ in range(self.npos):
            self.x[:] = self.x + self.fs_smoother.matvec(b - self.A0 @ self.x)

        return self.x.copy()
    
    def _matvec(self, b):
        return self._matvec_v1(b)

class MultiscaleIlu0Smoother2Levels(LinearOperator):
    def __init__(self, A0, P_lv1, R_lv1, P_lv2, R_lv2, npre=1, npost=1):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        self.fs_smoother = Ilu0Precond(A0.tocsr())
        A1 = R_lv1 @ A0 @ P_lv1
        self.A1 = A1
        self.smoother_lv1 = Ilu0Precond(A1.tocsr())
        
        A2 = R_lv2 @ A1 @ P_lv2

        # 2. Solvers Diretos para as Malhas Grossas
        self.solve_grossa2 = factorized(A2.tocsc())

        self.P1, self.R1 = P_lv1, R_lv1
        self.P2, self.R2 = P_lv2, R_lv2
        self.x = np.zeros(A0.shape[0])
        self.x1 = np.zeros(R_lv1.shape[0])
        self.x2 = np.zeros(R_lv2.shape[0])
        
        self.r0 = np.zeros(A0.shape[0])
        self.r1 = np.zeros(R_lv1.shape[0])
        self.b1 = np.zeros(R_lv1.shape[0])
        
        self.b2 = np.zeros(R_lv2.shape[0])
    
    def _matvec_v1(self, b):
        self.x[:] = 0
        self.r[:] = b
        self.r1[:] = 0
        self.d2[:] = 0
        
        for _ in range(self.npre):
            self.x[:] = self.x + self.fs_smoother.matvec(self.r)
            self.r[:] = b - self.A0 @ self.x
        
        
        self.b1[:] = self.R1 @ (self.r)
        self.r1[:] = self.b1
        
        for _ in range(self.npre):
            self.x1[:] = self.x1 + self.smoother_lv1.matvec(self.r1)
            self.r1[:] = self.b1 - self.A1 @ self.x1
        
        self.b2[:] = self.R2 @ self.r1
        
        self.x2[:] = self.solve_grossa2(self.b2)
        
        self.x1[:] = self.x1 + self.P2 @ self.x2
        
        self.r1[:] = self.b1 - self.A1 @ self.x1
        
        for _ in range(self.npost):
            self.x1[:] = self.x1 + self.smoother_lv1.matvec(self.r1)
            self.r1[:] = self.b1 - self.A1 @ self.x1
        
        self.x[:] = self.x + self.P1 @ self.x1
        self.r[:] = b - self.A0 @ self.x
        
        for _ in range(self.npost):
            self.x[:] = self.x + self.fs_smoother.matvec(self.r)
            self.r[:] = b - self.A0 @ self.x
        
        return self.x.copy()
    
    def _matvec(self, v):
        return self._matvec_v1(v)
        

class MultiScaleIlu0Smoother(LinearOperator):
    def __init__(self, A0, P1, R1, npre=0, npost=1):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        self.fs_smoother = Ilu0Precond(A0)

        A1_1 = R1 @ A0 @ P1

        # 2. Solvers Diretos para as Malhas Grossas
        # factorized() é um wrapper rápido para o splu do scipy
        self.solve_grossa1 = factorized(A1_1.tocsc())

        self.P1, self.R1 = P1, R1

        self.x = np.zeros(A0.shape[0])
        self.r0 = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])
    
    def _matvec_v1(self, b):
        
        for _ in range(self.npre):
            self.x[:] = self.x + self.fs_smoother.matvec(b - self.A0 @ self.x)
        
        self.d1[:] = self.R1 @ (b - self.A0 @ self.x)
        self.d1[:] = self.R1 @ b
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + self.P1 @ self.d1
        
        for _ in range(self.npost):
            self.x[:] = self.x + self.fs_smoother.matvec(b - self.A0 @ self.x)

        return self.x.copy()
    
    def _matvec_v2(self, b):
        self.x[:] = 0
        self.r0[:] = b
        
        for _ in range(self.npre):
            self.x[:] = self.x + self.fs_smoother.matvec(self.r0)
            self.r0[:] = b - self.A0 @ self.x
        
        
        self.d1[:]  = self.R1 @ self.r0
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + self.P1 @ self.d1
        self.r0[:] = b - self.A0 @ self.x
        
        for _ in range(self.npost):
            self.x[:] = self.x + self.fs_smoother.matvec(self.r0)
            self.r0[:] = b - self.A0 @ self.x
        
        return self.x.copy()
        
        
    
    def _matvec(self, b):
        # return self._matvec_v1(b)
        return self._matvec_v2(b)

class MultiscaleIlu1Smoother(MultiScaleIlu0Smoother):
    def __init__(self, A0, P1, R1, npre=0, npost=1, threshold=1e-16):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        self.fs_smoother = Ilu1Precond(A0, threshold=threshold)

        A1_1 = R1 @ A0 @ P1

        # 2. Solvers Diretos para as Malhas Grossas
        # factorized() é um wrapper rápido para o splu do scipy
        self.solve_grossa1 = factorized(A1_1.tocsc())

        self.P1, self.R1 = P1, R1

        self.x = np.zeros(A0.shape[0])
        self.r0 = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])


class MultiScaleRootnodeSmoother2OP(LinearOperator):
    def __init__(self, A0, P1, R1, P2, R2):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype

        self.fs_smoother = pyamg.rootnode_solver(A0).aspreconditioner()

        A1_1 = R1 @ A0 @ P1
        A1_2 = R2 @ A0 @ P2

        # 2. Solvers Diretos para as Malhas Grossas
        # factorized() é um wrapper rápido para o splu do scipy
        self.solve_grossa1 = factorized(A1_1.tocsc())
        self.solve_grossa2 = factorized(A1_2.tocsc())

        self.P1, self.R1 = P1, R1
        self.P2, self.R2 = P2, R2

        self.x = np.zeros(A0.shape[0])
        self.r = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])
        self.d2 = np.zeros(R2.shape[0])
    
    def _matvec_v1(self, b):
        
        self.d1[:] = self.R1 @ b
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.P1 @ self.d1
        
        self.r[:] = b - self.A0 @ self.x
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        
        self.r[:] = b - self.A0 @ self.x
        self.d2[:] = self.R2 @ self.r
        self.d2[:] = self.solve_grossa2(self.d2)
        self.x[:] = self.x + (self.P2 @ self.d2)
        
        self.r[:] = b - self.A0 @ self.x
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)

        return self.x.copy()
    
    def _matvec(self, b):
        return self._matvec_v1(b)

class MultilevelSolverMultiplicativePyamg(LinearOperator):
    
    def __init__(self, solver_1: MultilevelSolver, solver_2: MultilevelSolver, A0: sp.csr_matrix):
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.solver_1 = solver_1.aspreconditioner(cycle='V')
        self.solver_2 = solver_2.aspreconditioner(cycle='V')
        self.A0 = A0
        
        self.x = np.zeros(A0.shape[0])
        self.r0 = np.zeros(A0.shape[0])
        self.xf1 = np.zeros(A0.shape[0])
        self.xf2 = np.zeros(A0.shape[0])
    
    def _matvec_v1(self, b):
        self.x[:] = 0
        self.xf1[:] = 0
        self.xf2[:] = 0
        self.r0[:] = b
        
        self.xf1[:] = self.solver_1.matvec(self.r0)
        self.x[:] = self.x + self.xf1
        self.r0[:] = b - self.A0 @ self.x
        
        self.xf2[:] = self.solver_2.matvec(self.r0)
        self.x[:] = self.x + self.xf2
        # self.r0[:] = b - self.A0 @ self.x
        
        return self.x.copy()
    
    def _matvec(self, b):
        return self._matvec_v1(b)

class MultiscaleSolverPyamg(LinearOperator):
    
    def __init__(self, A0: sp.csr_matrix, OP1: sp.csr_matrix, OR1: sp.csr_matrix, finescale_smoother: LinearOperator, **kwargs):
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.A0 = A0
        self.fs_smoother = finescale_smoother
        
        self.x = np.zeros(A0.shape[0])
        self.r0 = np.zeros(A0.shape[0])
        self.xf1 = np.zeros(A0.shape[0])
        # self.preprocess()
        
        levels = []
        
        ## Nivel 0
        levels.append(pyamg.multilevel.multilevel_solver.level())
        levels[0].A = A0.tocsr()
        levels[0].P = OP1.tocsr()
        levels[0].R = OR1.tocsr()
        
        ## Nivel 1
        levels.append(pyamg.multilevel.multilevel_solver.level())
        levels[1].A = levels[0].R @ levels[0].A @ levels[0].P
        
        # ml = multilevel_solver(levels, coarse_solver='splu')
        ml1 = MultilevelSolver(levels, coarse_solver='splu')
        
        def no_op(A, x, b):
            pass
        
        def ilu_smoother(A: sp.csr_matrix, x: np.ndarray, b: np.ndarray):
            niter = 1
            for _ in range(niter):
                x[:] = x + self.fs_smoother.matvec(b - A @ x)
            return x
        
        ml1.levels[0].presmoother = no_op
        ml1.levels[0].postsmoother = ilu_smoother
        
        self.solver_1 = ml1.aspreconditioner()
        
    
    def preprocess(self):
        
        levels = []
        
        ## Nivel 0
        levels.append(pyamg.multilevel.multilevel_solver.level())
        levels[0].A = self.A0.tocsr()
        levels[0].P = self.OP1.tocsr()
        levels[0].R = self.OR1.tocsr()
        
        ## Nivel 1
        levels.append(pyamg.multilevel.multilevel_solver.level())
        levels[1].A = levels[0].R @ levels[0].A @ levels[0].P
        
        # ml = multilevel_solver(levels, coarse_solver='splu')
        ml1 = MultilevelSolver(levels, coarse_solver='splu')
        
        def no_op(A, x, b):
            pass
        
        def ilu_smoother(A: sp.csr_matrix, x: np.ndarray, b: np.ndarray):
            niter = 1
            for _ in range(niter):
                x[:] = x + self.fs_smoother.matvec(b - A @ x)
            return x
        
        ml1.levels[0].presmoother = no_op
        ml1.levels[0].postsmoother = ilu_smoother
        
        self.solver_1 = ml1.aspreconditioner()
    
    def _matvec(self, b):
        return self.solver_1.matvec(b)


class MultiscaleSolverPyamgIlu0(MultiscaleSolverPyamg):
    
    def __init__(self, A0: sp.csr_matrix, OP1: sp.csr_matrix, OR1: sp.csr_matrix, **kwargs):
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.A0 = A0
        self.fs_smoother = Ilu0Precond(A0)
        self.OP1 = OP1
        self.OR1 = OR1
        
        self.x = np.zeros(A0.shape[0])
        self.r0 = np.zeros(A0.shape[0])
        self.xf1 = np.zeros(A0.shape[0])
        self.preprocess()

class MultiscaleSolverPyamgIlu1(MultiscaleSolverPyamg):
    
    def __init__(self, A0: sp.csr_matrix, OP1: sp.csr_matrix, OR1: sp.csr_matrix, threshold=1e-16, **kwargs):
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.A0 = A0
        self.fs_smoother = Ilu1Precond(A0, threshold=threshold)
        self.OP1 = OP1
        self.OR1 = OR1
        
        self.x = np.zeros(A0.shape[0])
        self.r0 = np.zeros(A0.shape[0])
        self.xf1 = np.zeros(A0.shape[0])
        self.preprocess()

    
class MultilevelSolverMultiplicativePyamg_2OP(LinearOperator):
    
    def __init__(self, A0: sp.csr_matrix, OP1: sp.csr_matrix, OR1: sp.csr_matrix, OP2: sp.csr_matrix, OR2: sp.csr_matrix, finescale_smoother: LinearOperator, **kwargs):
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.A0 = A0
        
        self.solver_1 = MultiscaleSolverPyamg(A0, OP1, OR1, finescale_smoother)
        self.solver_2 = MultiscaleSolverPyamg(A0, OP2, OR2, finescale_smoother)
        
        self.x = np.zeros(A0.shape[0])
        self.r0 = np.zeros(A0.shape[0])       
    
    def _matvec_v1(self, b):
        self.x[:] = 0
        self.r0[:] = b
        
        self.x[:] = self.x + self.solver_1.matvec(self.r0)
        self.r0[:] = b - self.A0 @ self.x
        
        self.x[:] = self.x + self.solver_2.matvec(self.r0)
        
        return self.x.copy()
    
    def _matvec(self, b):
        return self._matvec_v1(b)

class MultilevelSolverMultiplicativePyamg_2OPIlu0(MultilevelSolverMultiplicativePyamg_2OP):
    
    def __init__(self, A0: sp.csr_matrix, OP1: sp.csr_matrix, OR1: sp.csr_matrix, OP2: sp.csr_matrix, OR2: sp.csr_matrix, **kwargs):
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.A0 = A0
        
        
        self.solver_1 = MultiscaleSolverPyamgIlu0(A0, OP1, OR1)
        self.solver_2 = MultiscaleSolverPyamgIlu0(A0, OP2, OR2)
        
        self.x = np.zeros(A0.shape[0])
        self.r0 = np.zeros(A0.shape[0])
        

class MultilevelSolverMultiplicativePyamg_2OPIlu1(MultilevelSolverMultiplicativePyamg_2OP):
    
    def __init__(self, A0: sp.csr_matrix, OP1: sp.csr_matrix, OR1: sp.csr_matrix, OP2: sp.csr_matrix, OR2: sp.csr_matrix, threshold=1e-16, **kwargs):
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.A0 = A0
        
        self.solver_1 = MultiscaleSolverPyamgIlu1(A0, OP1, OR1)
        self.solver_2 = MultiscaleSolverPyamgIlu1(A0, OP2, OR2)
        
        self.x = np.zeros(A0.shape[0])
        self.r0 = np.zeros(A0.shape[0])

class MultiscaleSolverPyamg2levels(LinearOperator):
    def __init__(self, A0, P_lv0, R_lv0, P_lv1, R_lv1, finescale_smoother: LinearOperator, lv1_smoother: LinearOperator, **kwargs):
        
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.A0 = A0
        self.x = np.zeros(A0.shape[0])
        self.r = np.zeros(A0.shape[0])
        
        
        levels = []
        
        ## Nivel 0
        levels.append(pyamg.multilevel.multilevel_solver.level())
        levels[0].A = self.A0.tocsr()
        levels[0].P = P_lv0.tocsr()
        levels[0].R = R_lv0.tocsr()
        
        ## Nivel 1
        levels.append(pyamg.multilevel.multilevel_solver.level())
        levels[1].A = levels[0].R @ levels[0].A @ levels[0].P
        levels[1].P = P_lv1.tocsr()
        levels[1].R = R_lv1.tocsr()
        
        ## Nivel 2
        levels.append(pyamg.multilevel.multilevel_solver.level())
        levels[2].A = levels[1].R @ levels[1].A @ levels[1].P
        
        def smoother_lv0(A: sp.csr_matrix, x: np.ndarray, b: np.ndarray):
            niter = 1
            for _ in range(niter):
                r = b-A@x
                x[:] = x + finescale_smoother.matvec(r)
            return x
        
        def smoother_lv1(A: sp.csr_matrix, x: np.ndarray, b: np.ndarray):
            niter = 1
            for _ in range(niter):
                r = b-A@x
                x[:] = x + lv1_smoother.matvec(r)
            return x
        
        levels[0].presmoother = smoother_lv0
        levels[0].postsmoother = smoother_lv0
        levels[1].presmoother = smoother_lv1
        levels[1].postsmoother = smoother_lv1
        
        ml = MultilevelSolver(levels, coarse_solver='splu')
        
        self.solver = ml.aspreconditioner()
        
    def _matvec(self, b):
        return self.solver.matvec(b)
        

class MultiscaleIlu0SmootherAlgorithmic(MultiScaleIlu0Smoother):
     def __init__(self, A0, P1, R1, npre=0, npost=1, zeta=0.0001, weight=1.0):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        self.fs_smoother = Ilu0Precond(A0)

        A1_1 = R1 @ A0 @ P1
        ag = AlgorithimicMonotone()
        A1 = ag.get_monotone_matrix(A1_1, epsilon=zeta, w=weight)
        

        # 2. Solvers Diretos para as Malhas Grossas
        # factorized() é um wrapper rápido para o splu do scipy
        self.solve_grossa1 = factorized(A1.tocsc())

        self.P1, self.R1 = P1, R1

        self.x = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])
        
class MultiscaleIlu1SmootherAlgorithmic(MultiScaleIlu0Smoother):
     def __init__(self, A0, P1, R1, npre=0, npost=1, zeta=0.0001, weight=1.0, threshold=1e-16, ilu_factor=1.5):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        self.fs_smoother = Ilu1Precond(A0, threshold=threshold)

        A1_1 = R1 @ A0 @ P1
        ag = AlgorithimicMonotone()
        A1 = ag.get_monotone_matrix(A1_1, epsilon=zeta, w=weight)
        

        # 2. Solvers Diretos para as Malhas Grossas
        # factorized() é um wrapper rápido para o splu do scipy
        self.solve_grossa1 = factorized(A1.tocsc())

        self.P1, self.R1 = P1, R1

        self.x = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])   

   
class MultiScaleSmoother(LinearOperator):
    def __init__(self, A0, P1, R1, finescale_smoother: LinearOperator, npre=0, npost=1, *args, **kwargs):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        self.fs_smoother = finescale_smoother

        A1_1 = R1 @ A0 @ P1
        self.solve_grossa1 = factorized(A1_1.tocsc())

        self.P1, self.R1 = P1, R1

        self.x = np.zeros(A0.shape[0])
        self.r0 = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])
    
    def _matvec_v2(self, b):
        self.x[:] = 0
        self.r0[:] = b
        
        for _ in range(self.npre):
            self.x[:] = self.x + self.fs_smoother.matvec(self.r0)
            self.r0[:] = b - self.A0 @ self.x
        
        
        self.d1[:]  = self.R1 @ self.r0
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + self.P1 @ self.d1
        self.r0[:] = b - self.A0 @ self.x
        
        for _ in range(self.npost):
            self.x[:] = self.x + self.fs_smoother.matvec(self.r0)
            self.r0[:] = b - self.A0 @ self.x
        
        return self.x.copy()
        
        
    
    def _matvec(self, b):
        # return self._matvec_v1(b)
        return self._matvec_v2(b)

class MultiScaleIlu0Smoother2OP_v2(LinearOperator):
    def __init__(self, A0, P1, R1, P2, R2, npre=0, npost=1, **kwargs):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        A1_1 = R1 @ A0 @ P1
        A1_2 = R2 @ A0 @ P2

        self.solve_grossa1 = factorized(A1_1.tocsc())
        self.solve_grossa2 = factorized(A1_2.tocsc())

        self.P1, self.R1 = P1, R1
        self.P2, self.R2 = P2, R2
        self.x = np.zeros(A0.shape[0])
        self.r = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])
        self.d2 = np.zeros(R2.shape[0])
        self.fs_smoother = self.create_fs_smoother(A0, **kwargs)
        self.count = 0
        
    
    def create_fs_smoother(self, A0, **kwargs) -> LinearOperator:
        return  Ilu0Precond(A0, **kwargs) 

    def _matvec_v0(self, b):
        
        self.x[:] = self.fs_smoother.matvec(b)

        # r1 = b - self.A0 @ x
        # d1 = self.R1 @ r1
        # c1 = self.solve_grossa1(d1)
        # x = x + (self.P1 @ c1)
        self.r[:] = b - self.A0 @ self.x
        self.d1[:] = self.R1 @ self.r
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + (self.P1 @ self.d1)

        # r2 = b - self.A0 @ x
        # d2 = self.R2 @ r2
        # c2 = self.solve_grossa2(d2)
        # x[:] = x + (self.P2 @ c2)
        self.r[:] = b - self.A0 @ self.x
        self.d2[:] = self.R2 @ self.r
        self.d2[:] = self.solve_grossa1(self.d2)
        self.x[:] = self.x + (self.P2 @ self.d2)

        return self.x.copy()
    
    def _matvec_v1(self, b):
        
        self.d1[:] = self.R1 @ b
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.P1 @ self.d1
        
        self.r[:] = b - self.A0 @ self.x
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        
        self.r[:] = b - self.A0 @ self.x
        self.d2[:] = self.R2 @ self.r
        self.d2[:] = self.solve_grossa2(self.d2)
        self.x[:] = self.x + (self.P2 @ self.d2)
        
        self.r[:] = b - self.A0 @ self.x
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)

        return self.x.copy()
    
    def _matvec_v2(self, b):
        self.x[:] = 0
        self.r[:] = b
        
        # for _ in range(self.npre):
        #     self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        #     self.r[:] = b - self.A0 @ self.x
        
        self.d1[:]  = self.R1 @ self.r
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + self.P1 @ self.d1
        self.r[:] = b - self.A0 @ self.x    
    
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        self.r[:] = b - self.A0 @ self.x
        
        # for _ in range(self.npost):
        #     self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        #     self.r[:] = b - self.A0 @ self.x
        
        self.d2[:] = self.R2 @ self.r
        self.d2[:] = self.solve_grossa2(self.d2)
        self.x[:] = self.x + self.P2 @ self.d2
        self.r[:] = b - self.A0 @ self.x
        
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        self.r[:] = b - self.A0 @ self.x
        
        # for _ in range(self.npost):
        #     self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        #     self.r[:] = b - self.A0 @ self.x
        
        # print(f'Count: {self.count}')
        # self.count += 1
        
        return self.x.copy()
        
        
        
    
    def _matvec(self, b):
        # return self._matvec_v1(b)
        return self._matvec_v2(b)
    

class MultiScaleIlu1Smoother2OP_v2(MultiScaleIlu0Smoother2OP_v2):
    def __init__(self, A0, P1, R1, P2, R2, npre=0, npost=1, **kwargs):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        A1_1 = R1 @ A0 @ P1
        A1_2 = R2 @ A0 @ P2

        self.solve_grossa1 = factorized(A1_1.tocsc())
        self.solve_grossa2 = factorized(A1_2.tocsc())

        self.P1, self.R1 = P1, R1
        self.P2, self.R2 = P2, R2
        self.x = np.zeros(A0.shape[0])
        self.r = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])
        self.d2 = np.zeros(R2.shape[0])
        self.fs_smoother = self.create_fs_smoother(A0, **kwargs)
        
    
    def create_fs_smoother(self, A0, **kwargs) -> LinearOperator:
        return  Ilu1Precond(A0, **kwargs)



class MultiScaleSmoother2OP(LinearOperator):
    def __init__(self, A0, P1, R1, P2, R2, finescale_smoother: LinearOperator, npre=0, npost=1, **kwargs):
        
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost
        
        A1_1 = R1 @ A0 @ P1
        A1_2 = R2 @ A0 @ P2
        self.solve_grossa1 = factorized(A1_1.tocsc())
        self.solve_grossa2 = factorized(A1_2.tocsc())
        
        self.P1, self.R1 = P1, R1
        self.P2, self.R2 = P2, R2
        
        self.x = np.zeros(A0.shape[0])
        self.r = np.zeros(A0.shape[0])        
        self.d1 = np.zeros(R1.shape[0])
        self.d2 = np.zeros(R2.shape[0])
        self.fs_smoother = finescale_smoother
        self.count = 0
        
        
    
    def _matvec_v2(self, b):
        self.x[:] = 0
        self.r[:] = b
        
        # for _ in range(self.npre):
        #     self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        #     self.r[:] = b - self.A0 @ self.x
        
        self.d1[:]  = self.R1 @ self.r
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + self.P1 @ self.d1
        self.r[:] = b - self.A0 @ self.x    
    
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        # self.x[:] = self.x + self.fs_smoother @ self.r
        self.r[:] = b - self.A0 @ self.x
        
        # for _ in range(self.npost):
        #     self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        #     self.r[:] = b - self.A0 @ self.x
        
        self.d2[:] = self.R2 @ self.r
        self.d2[:] = self.solve_grossa2(self.d2)
        self.x[:] = self.x + self.P2 @ self.d2
        self.r[:] = b - self.A0 @ self.x
        
        self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        # self.x[:] = self.x + self.fs_smoother @ self.r
        self.r[:] = b - self.A0 @ self.x
        
        # for _ in range(self.npost):
        #     self.x[:] = self.x + self.fs_smoother.matvec(self.r)
        #     self.r[:] = b - self.A0 @ self.x
        
        # print(f'Count: {self.count}')
        # self.count += 1
        
        return self.x.copy()        
        
    
    def _matvec(self, b):
        self._matvec_v2(b)
        # self._matvec_v3(b)
        
class MultiscaleSmoother2Levels(LinearOperator):
    def __init__(self, A0, P_lv1, R_lv1, P_lv2, R_lv2, finescale_smoother: LinearOperator, lv1_smoother: LinearOperator, npre=1, npost=1, **kwargs):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        self.fs_smoother = finescale_smoother
        self.smoother_lv1 = lv1_smoother
        A1 = R_lv1 @ A0 @ P_lv1
        A2 = R_lv2 @ A1 @ P_lv2
        self.A1 = A1
        # 2. Solvers Diretos para as Malhas Grossas
        self.solve_grossa2 = factorized(A2.tocsc())

        self.P1, self.R1 = P_lv1, R_lv1
        self.P2, self.R2 = P_lv2, R_lv2
        self.x = np.zeros(A0.shape[0])
        self.x1 = np.zeros(R_lv1.shape[0])
        self.x2 = np.zeros(R_lv2.shape[0])
        
        self.r0 = np.zeros(A0.shape[0])
        self.r1 = np.zeros(R_lv1.shape[0])
        self.b1 = np.zeros(R_lv1.shape[0])
        self.b2 = np.zeros(R_lv2.shape[0])
    
    def _matvec_v1(self, b):
        self.x[:] = 0
        self.r0[:] = b
        self.r1[:] = 0
        
        for _ in range(self.npre):
            self.x[:] = self.x + self.fs_smoother.matvec(self.r0)
            self.r0[:] = b - self.A0 @ self.x
        
        
        self.b1[:] = self.R1 @ (self.r0)
        self.r1[:] = self.b1
        
        for _ in range(self.npre):
            self.x1[:] = self.x1 + self.smoother_lv1.matvec(self.r1)
            self.r1[:] = self.b1 - self.A1 @ self.x1
        
        self.b2[:] = self.R2 @ self.r1
        
        self.x2[:] = self.solve_grossa2(self.b2)
        
        self.x1[:] = self.x1 + self.P2 @ self.x2
        
        self.r1[:] = self.b1 - self.A1 @ self.x1
        
        for _ in range(self.npost):
            self.x1[:] = self.x1 + self.smoother_lv1.matvec(self.r1)
            self.r1[:] = self.b1 - self.A1 @ self.x1
        
        self.x[:] = self.x + self.P1 @ self.x1
        self.r0[:] = b - self.A0 @ self.x
        
        for _ in range(self.npost):
            self.x[:] = self.x + self.fs_smoother.matvec(self.r0)
            self.r0[:] = b - self.A0 @ self.x
        
        return self.x.copy()
    
    def _matvec(self, v):
        return self._matvec_v1(v)

    

    
    
class MultiScaleRootnodeSmoother(LinearOperator):
    def __init__(self, A0, P1, R1, npre=0, npost=1):
        self.A0 = A0
        self.shape = A0.shape
        self.dtype = A0.dtype
        self.npre = npre
        self.npost = npost

        self.fs_smoother = pyamg.rootnode_solver(A0).aspreconditioner()

        A1_1 = R1 @ A0 @ P1

        # 2. Solvers Diretos para as Malhas Grossas
        # factorized() é um wrapper rápido para o splu do scipy
        self.solve_grossa1 = factorized(A1_1.tocsc())

        self.P1, self.R1 = P1, R1

        self.x = np.zeros(A0.shape[0])
        self.r = np.zeros(A0.shape[0])
        self.d1 = np.zeros(R1.shape[0])

    def _matvec_v0(self, b):
        
        self.x[:] = self.ilu.apply(b)

        # r1 = b - self.A0 @ x
        # d1 = self.R1 @ r1
        # c1 = self.solve_grossa1(d1)
        # x = x + (self.P1 @ c1)
        self.r[:] = b - self.A0 @ self.x
        self.d1[:] = self.R1 @ self.r
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + (self.P1 @ self.d1)

        return self.x.copy()
    
    def _matvec_v1(self, b):
        
        for _ in range(self.npre):
            self.x[:] = self.x + self.fs_smoother.matvec(b - self.A0 @ self.x)
        
        self.d1[:] = self.R1 @ (b - self.A0 @ self.x)
        self.d1[:] = self.solve_grossa1(self.d1)
        self.x[:] = self.x + self.P1 @ self.d1
        
        for _ in range(self.npost):
            self.x[:] = self.x + self.fs_smoother.matvec(b - self.A0 @ self.x)

        return self.x.copy()
    
    def _matvec(self, b):
        return self._matvec_v1(b)