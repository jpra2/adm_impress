from packs.mpfa_methods.mesh_preprocess import preprocess_mesh
from packs import defpaths
import os


def run():
    fine_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh0.msh')
    fine_mesh_properties_name = 'fine_properties_uns'
    coarse_mesh_path = os.path.join(defpaths.unstructured_coarse_test_mesh_folder, 'mesh1.msh')
    coarse_mesh_properties_name = 'coarse_properties_uns'

    fine_mesh_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)



    import pdb; pdb.set_trace()


