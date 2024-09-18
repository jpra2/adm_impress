from packs.manager.meshio_wrapper import MeshioWrapper
import os
from packs import defpaths

def run():
    
    # mesh_path = os.path.join(defpaths.mesh, defpaths.unstructured_coarse_test_mesh_folder, 'mesh0.msh')
    mesh_path = os.path.join(defpaths.mesh, 'quad_test.msh')
    mesh_data = MeshioWrapper(mesh_path)

    for tag in mesh_data.physical_tags:
        print(mesh_data.get_elements_by_physical_tag(tag))
    
    print('end test')

