from ..solvers.solvers_scipy.solver_sp import SolverSp
# from ..solvers.solvers_trilinos.solvers_tril import solverTril
from ..direct_solution.biphasic.biphasic1 import Biphasic
from .init_simulation import rodar

# __all__ = []

M = rodar.M

biphasic = Biphasic(M)
solver = SolverSp()

biphasic.get_transmissibility_matrix_without_contours()
biphasic.get_transmissibility_matrix()
biphasic.get_RHS_term()

# x = solver.direct_solver(biphasic.datas['T'], biphasic.datas['b'])
# biphasic.get_solution(x)
biphasic.get_solution(solver.direct_solver(biphasic.datas['T'], biphasic.datas['b']))
biphasic.get_flux_faces_and_volumes()

import pdb; pdb.set_trace()
