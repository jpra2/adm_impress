from impress import FineScaleMesh as mesh
from packs.vugs_generator import  vugs_generator
import numpy as np
file="mesh/isotropic/90x90x1.msh"
M = mesh(file)
vugs_generator.GenarateVugs(M)
