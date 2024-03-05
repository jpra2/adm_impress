from packs.manager.meshio_wrapper import MeshioWrapper
from pymoab import core, types, rng, topo_util
from packs.manager.meshmanager import MeshProperty
import numpy as np

def create_vertices(points: np.ndarray, mb: core.Core) -> np.ndarray:
    verts = mb.create_vertices(points.flatten())
    return verts

def create_triangles(verts: np.ndarray, triangle_points: np.ndarray, mb: core.Core):
    alltriangles = [verts[i] for i in triangle_points]
    triangles_moab = mb.create_elements(types.MBTRI, alltriangles)

def create_meshproperties_from_meshio(mesh_path:str, mesh_properties_name: str):
    
    meshio_data = MeshioWrapper(mesh_path)

    mb = core.Core()
    root_set = mb.get_root_set()
    mtu = topo_util.MeshTopoUtil(mb)

    all_nodes = create_vertices(meshio_data.points_centroids, mb)
    create_triangles(all_nodes, meshio_data.triangles_points, mb)
    all_faces = mb.get_entities_by_dimension(root_set, 2)

    initial_edges = mb.get_entities_by_dimension(root_set, 1)
    mtu.construct_aentities(all_nodes)
    all_edges = mb.get_entities_by_dimension(root_set, 1)
    edges_centroids_meshio = meshio_data.lines_centroids

    all_nodes_coords = mb.get_coords(all_nodes).reshape(len(all_nodes), 3)



    print('fim')



















