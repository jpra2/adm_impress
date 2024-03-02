from packs.manager.meshio_wrapper import MeshioWrapper
import os
from packs import defpaths

def run():
    
    mesh_path = os.path.join(defpaths.mesh, defpaths.unstructured_coarse_test_mesh_folder, 'mesh0.msh')
    mesh_data = MeshioWrapper(mesh_path)
    print('oi')

