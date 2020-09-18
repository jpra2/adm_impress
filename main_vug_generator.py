from packs.vugs_generator.vugs_generator import VugGenerator
import numpy as np
input_file = "mesh/100x100x1_unst.msh"
vug_gen = VugGenerator(input_file, (1.0, 5.0), num_ellipsoids=20, num_fractures=5)
vug_gen.run()
vug_gen.write_file()
