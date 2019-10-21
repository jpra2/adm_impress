import scipy.sparse as sp


class SolverSp:

    def __init__(self, A, b, x=None):
        self.A = A.tocsc()
        self.b = b

    def direct_solver(self):
        from scipy.sparse import linalg
        b = self.b
        A = self.A

        solution = linalg.spsolve(A,b)

        self.x = solution

    def lu_solver(self):
        from scipy.sparse import linalg
        b = self.b
        A = self.A

        LU = linalg.splu(A)
        solution = LU.solve(b)

        self.LU = LU
        self.x = solution
