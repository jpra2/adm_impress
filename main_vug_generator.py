from impress import FineScaleMesh as mesh
from packs.vugs_generator import  vugs_generator
import numpy as np
file="mesh/27x27x27.msh"
M = mesh(file)
vugs_generator.GenarateVugs(M)
