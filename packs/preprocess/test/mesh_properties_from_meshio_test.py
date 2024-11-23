from packs import defpaths
from packs.preprocess.create_mesh_properties_from_meshiowrapper import create_meshproperties_from_meshio, create_meshproperties_from_meshio_if_not_exists
import os

def run():

    mesh_path = os.path.join(
        defpaths.unstructured_coarse_test_mesh_folder,
        'mesh0.msh'
    )

    mesh_properties_name = 'mesh0_test'

    mesh_properties = create_meshproperties_from_meshio(mesh_path, mesh_properties_name)

    # mesh_properties = create_meshproperties_from_meshio_if_not_exists(mesh_path, mesh_properties_name)

    return mesh_properties
