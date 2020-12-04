from impress import FineScaleMesh as mesh
from packs.vugs_generator import  vugs_generator
import numpy as np
file="mesh/100x100x1_unst.msh"
M = mesh(file)
vugs_generator.GenarateVugs(M)
