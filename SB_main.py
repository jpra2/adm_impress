from impress import FineScaleMesh as mesh
from packs.stokes_brinkman_3d import stokes_brinkman as stokes

file='./mesh/9x9x9.msh'
M = mesh(file)
import pdb; pdb.set_trace()
stokes.stokes_solver(M)
