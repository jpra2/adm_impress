from packs.running.run_simulation import RunSimulation
rodar = RunSimulation(state=5)
M = rodar.M

from packs.direct_solution.monophasic.monophasic1 import Monophasic

m1 = Monophasic(M)
m1.get_transmissibility_matrix_without_contours()
m1.get_transmissibility_matrix()
m1.get_RHS_term()

from packs.solvers.solvers_scipy.solver_sp import SolverSp
solver = SolverSp()
x = solver.direct_solver(m1.datas['T'], m1.datas['b'])
m1.get_solution(x)
m1.get_flux_faces_and_volumes()
import pdb
pdb.set_trace()
M.data.update_variables_to_mesh()




import pdb; pdb.set_trace()





import pdb; pdb.set_trace()








# import pickle
#
# from tcc.load_save_initialize.load_infos import LoadInfos
# from tcc.dual_mesh.create_dual_mesh import DualMesh1
# import numpy as np
# from . import directories
# import pdb; pdb.set_trace()
#
# # file_name = os.path.join(path_flying, 'mesh_obj.txt')
# # with open(file_name, 'wb') as handle:
# #     pickle.dump(M, handle)
#
#
#
# LoadInfos(M)
# DualMesh1(M)
#
# import pdb; pdb.set_trace()
