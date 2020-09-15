# from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as mesh
# from packs.vugs_generator.vugs_generator_old import GenarateVugs
from packs.vugs_generator.vugs_generator import VugGenerator
import numpy as np
input_file = "mesh/27x27x27.msh"
# M = mesh(file)
# vug_gen = GenarateVugs(M)
vug_gen = VugGenerator(input_file, (50.0, 100.0))
vug_gen.run()
vug_gen.write_file()
