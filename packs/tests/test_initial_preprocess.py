import os
from packs import defpaths
from packs.mpfa_methods.mesh_preprocess import preprocess_mesh

def run():

    rel_path = os.path.join(
            defpaths.unstructured_coarse_test_mesh_folder,
            'brazil'
        )

    fine_mesh_path = os.path.join(rel_path, 'brazilf.msh')
    fine_mesh_properties_name = 'brazilf' 
    fine_mesh_path_v4 = os.path.join(rel_path, 'brazilf_v4.msh')

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)

    print('foi')