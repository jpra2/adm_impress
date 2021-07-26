from impress import FineScaleMesh as mesh
from implicit_impress.jacobian.impress_assembly import assembly
from implicit_impress.newton import newton


M=mesh('./mesh/9x9x9.h5m')
newton(M)
# J=symbolic_J().J
# assembly(M)

import pdb; pdb.set_trace()
