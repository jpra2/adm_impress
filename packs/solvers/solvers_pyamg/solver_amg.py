import pyamg

class SolverAMG:

    def __init__(self):
        pass

    def smoothed_aggregation_solver(self, A, b, TOL=1e-10):
        ml = pyamg.smoothed_aggregation_solver(A)
        solution = ml.solve(b, tol=TOL)
        return solution
