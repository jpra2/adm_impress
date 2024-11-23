import os
from packs import defpaths
from packs.mpfa_methods.mesh_preprocess import preprocess_mesh
from packs.manager import MeshData, MeshProperty

def run():

    # rel_path = os.path.join(
    #         'malhas_diss',
    #         'brazil'
    #     )
    
    rel_path = 'malhas_diss'

    # fine_mesh_path = os.path.join(rel_path, 'square_uns.msh')
    fine_mesh_path = os.path.join(rel_path, 'square_uns_coarse.msh')
    # fine_mesh_properties_name = 'square_uns'
    fine_mesh_properties_name = 'square_uns_coarse'
    # fine_mesh_path_v4 = os.path.join(rel_path, 'brazilf_v4.msh')

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)

    mesh_data = MeshData(mesh_path=fine_mesh_path)

    print('foi')