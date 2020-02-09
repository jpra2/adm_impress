from impress import FineScaleMesh as mesh
from packs.stokes_brinkman_3d import stokes_brinkman as stokes
# from implicit_impress.jacobian.impress_assembly import assembly
# from implicit_impress.newton import newton
M=mesh('./mesh/3x3x3.msh')
stokes.stokes_solver(M)

# newton(M)
# J=symbolic_J().J
# assembly(M)

import pdb; pdb.set_trace()
