import os
from packs import defpaths
from packs.preprocess.create_mesh import createMesh
import numpy as np

mesh = createMesh()
mesh.create_fine_vertices()
mesh.create_elements()
mesh.export_mesh()

def test_weights():
    ## verify if exists the mesh test
    mesh_test_name = 'test_weights_2d_structured.h5m'
    mesh_path = os.path.join(defpaths.mesh, mesh_test_name)
    if os.path.exists(mesh_path):
        pass
    else:
        ## create mesh
        # starting point
        starting_point = np.array([0.0, 0.0, 0.0])

        # tamanho de cada bloco da malha fina
        block_size = np.array([20, 10, 2])

        # tamanho total da malha
        block_number = np.array([60, 220, 1])

        # mesh_size: [9.0, 1.0, 9.0]

        # nome da malha para salvar

        mesh_name = mesh_test_name


        mesh = createMesh()

        import pdb; pdb.set_trace()



        mesh.create_fine_vertices()
        mesh.create_elements()
        mesh.export_mesh()


        import pdb; pdb.set_trace()