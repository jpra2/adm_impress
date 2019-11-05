from ..solvers.solvers_scipy.solver_sp import SolverSp
from ..direct_solution.monophasic.monophasic1 import Monophasic
import pdb
from .init_simulation import rodar

M = rodar.M

# prep1 = Preprocess1()
# prep1.set_saturation_regions(M)

m1 = Monophasic(M)
m1.get_transmissibility_matrix_without_contours()
m1.get_transmissibility_matrix()
m1.get_RHS_term()


solver = SolverSp()
x = solver.direct_solver(m1.datas['T'], m1.datas['b'])
m1.get_solution(x)
m1.get_flux_faces_and_volumes()

M.data.update_variables_to_mesh()
m1.export_datas_to_npz()
