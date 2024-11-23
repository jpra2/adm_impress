import pyfvm
from pyfvm.form_language import *
import meshzoo
from scipy.sparse import linalg
import meshplex
import time
import pyamg
class Poisson(object):
    def apply(self, u):
        return integrate(lambda x: -n_dot_grad(u(x)), dS) \
             - integrate(lambda x: 1.0, dV)

    def dirichlet(self, u):
        return [(lambda x: u(x) - 0.0, Boundary())]

# Create mesh using meshzo
vertices, cells = meshzoo.rectangle(0.0, 2.0, 0.0, 1.0, 401, 201)
mesh = meshplex.MeshTri(vertices, cells)
t1=time.time()
matrix, rhs = pyfvm.discretize_linear(Poisson(), mesh)
u = linalg.spsolve(matrix, rhs)
print(time.time()-t1)
t2=time.time()

ml = pyamg.smoothed_aggregation_solver(matrix)
u2 = ml.solve(rhs, tol=1e-10)
import pdb; pdb.set_trace()
print(time.time()-t2)
#mesh.write('out.vtk', point_data={'u': u})
import pdb; pdb.set_trace()
