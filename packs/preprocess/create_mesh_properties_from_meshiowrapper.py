from packs.manager.meshio_wrapper import MeshioWrapper
from pymoab import core, types, rng, topo_util
from packs.manager.meshmanager import MeshProperty
import numpy as np
from packs import defpaths
from packs.mpfa_methods.mesh_preprocess import preprocess_mesh
import os
from packs.utils.test_functions import test_mesh_path

def _create_vertices(points: np.ndarray, mb: core.Core) -> np.ndarray:
    verts = mb.create_vertices(points.flatten())
    return verts

def _create_triangles(verts: np.ndarray, triangle_points: np.ndarray, mb: core.Core):
    alltriangles = [verts[i] for i in triangle_points]
    triangles_moab = mb.create_elements(types.MBTRI, alltriangles)

def _create_flying_mesh(mesh_path):
    meshio_data = MeshioWrapper(mesh_path)
    flying_mesh_path = os.path.join(
        defpaths.remove_folder,
        'flying_mesh.h5m'
    )

    mb = core.Core()
    # root_set = mb.get_root_set()
    mtu = topo_util.MeshTopoUtil(mb)

    all_nodes = _create_vertices(meshio_data.points_centroids, mb)
    _create_triangles(all_nodes, meshio_data.triangles_points, mb)
    mtu.construct_aentities(all_nodes)
    mb.write_file(flying_mesh_path)
    return flying_mesh_path

def _verify_mesh_properties_exists(mesh_properties_name):
    mesh_properties = MeshProperty()
    mesh_properties.insert_mesh_name(mesh_properties_name)
    return mesh_properties.exists()

def create_meshproperties_from_meshio(mesh_path:str, mesh_properties_name: str):
    flying_mesh_path = _create_flying_mesh(mesh_path)
    mesh_properties = preprocess_mesh(flying_mesh_path, mesh_properties_name)
    
    return mesh_properties

def create_meshproperties_from_meshio_if_not_exists(mesh_path:str, mesh_properties_name: str):
    if _verify_mesh_properties_exists(mesh_properties_name):
        mesh_properties = MeshProperty()
        mesh_properties.insert_mesh_name(mesh_properties_name)
        mesh_properties.load_data()
    else:
        mesh_properties = create_meshproperties_from_meshio(mesh_path, mesh_properties_name)
    
    return mesh_properties
    



    



















