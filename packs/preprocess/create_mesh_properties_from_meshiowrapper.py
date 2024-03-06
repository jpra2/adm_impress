from packs.manager.meshio_wrapper import MeshioWrapper
from pymoab import core, types, rng, topo_util
from packs.manager.meshmanager import MeshProperty
import numpy as np
from packs import defpaths
from packs.mpfa_methods.mesh_preprocess import preprocess_mesh
import os

def create_vertices(points: np.ndarray, mb: core.Core) -> np.ndarray:
    verts = mb.create_vertices(points.flatten())
    return verts

def create_triangles(verts: np.ndarray, triangle_points: np.ndarray, mb: core.Core):
    alltriangles = [verts[i] for i in triangle_points]
    triangles_moab = mb.create_elements(types.MBTRI, alltriangles)

def create_flying_mesh(mesh_path):
    meshio_data = MeshioWrapper(mesh_path)
    flying_mesh_path = os.path.join(
        defpaths.flying,
        'flying_mesh.h5m'
    )

    mb = core.Core()
    # root_set = mb.get_root_set()
    mtu = topo_util.MeshTopoUtil(mb)

    all_nodes = create_vertices(meshio_data.points_centroids, mb)
    create_triangles(all_nodes, meshio_data.triangles_points, mb)
    mtu.construct_aentities(all_nodes)
    mb.write_file(flying_mesh_path)
    mb.write_file(os.path.join(defpaths.flying, 'teste_mesh0.vtk'))
    return flying_mesh_path

def create_meshproperties_from_meshio(mesh_path:str, mesh_properties_name: str):
    flying_mesh_path = create_flying_mesh(mesh_path)
    mesh_properties = preprocess_mesh(flying_mesh_path, mesh_properties_name)
    print('fim')
    
    return mesh_properties
    



    



















