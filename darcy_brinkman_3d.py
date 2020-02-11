from impress import FineScaleMesh as mesh
from packs.stokes_brinkman_3d import stokes_brinkman as stokes

M=mesh('./mesh/3x3x3.msh')
stokes.stokes_solver(M)

import pdb; pdb.set_trace()
