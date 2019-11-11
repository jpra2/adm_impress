from ..solvers.solvers_scipy.solver_sp import SolverSp
from ..direct_solution.biphasic.biphasic1 import Biphasic
from .init_simulation import rodar

M = rodar.M

biphasic = Biphasic(M)
