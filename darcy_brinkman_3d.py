from impress import FineScaleMesh as mesh
from packs.stokes_brinkman_3d import stokes_brinkman as stokes
# import update_inputs


M=mesh('./mesh/10x10x10.msh')
stokes.stokes_solver(M)

import pdb; pdb.set_trace()
