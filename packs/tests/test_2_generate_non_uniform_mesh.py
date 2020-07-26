import numpy as np
from packs.preprocess import createNonUniformMesh
from packs.tests.test_1_load_m_object import load_m_object
import pdb

mesh_name = 'mesh/10x10x4.h5m'
# x_points = np.linspace(0, 500, 501) #10 blocos em x e y
# y_points = x_points
# z_points = np.array([0.0, 15.24]) # 4 blocos em z
# generator = createNonUniformMesh()
# generator.create_mesh(x_points, y_points, z_points, mesh_name)

pdb.set_trace()
M = load_m_object(mesh_name)
