from impress import FineScaleMesh as mesh
from implicit_impress.jacobian.symbolic_jacobian import symbolic_J
from implicit_impress.jacobian.impress_assembly import assembly
M=mesh('./mesh/5x5x5.msh')
J=symbolic_J().J
assembly(M)

import pdb; pdb.set_trace()
