import numpy as np
from packs.preprocess import createNonUniformMesh
from packs.tests.test_1_load_m_object import load_m_object
import pdb

mesh_name = '10x10x4.h5m'
x_points = np.linspace(0, 804.6, 11) #10 blocos em x e y
y_points = x_points
z_points = np.array([0.0, 15.24, 30.48, 39.62, 48.76]) # 4 blocos em z
generator = createNonUniformMesh()
generator.create_mesh(x_points, y_points, z_points, mesh_name)

M = load_m_object(mesh_name)
