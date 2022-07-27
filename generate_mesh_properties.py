from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh
from pymoab import core, types, rng, topo_util
# from pymoab import skinner as sk
import numpy as np
import pandas as pd
from numba import jit

from packs.data_class.data_manager import DataManager

def init_mesh(mesh_path: str):
    mb = core.Core()
    mtu = topo_util.MeshTopoUtil(mb)
    mb.load_file(mesh_path)
    root_set = mb.get_root_set()
    return mb, mtu, root_set

def init_mesh_entities(mb: core.Core, mtu: topo_util.MeshTopoUtil):
    all_volumes = mb.get_entities_by_dimension(0, 3)
    all_nodes = mb.get_entities_by_dimension(0, 0)
    mtu.construct_aentities(all_nodes)
    all_faces = mb.get_entities_by_dimension(0, 2)
    all_edges = mb.get_entities_by_dimension(0, 1)
    return all_volumes, all_faces, all_edges, all_nodes

def init_properties(faces: np.ndarray, volumes: np.ndarray, nodes: np.ndarray, all_volumes: rng, all_faces: rng, all_nodes: rng, mtu: topo_util.MeshTopoUtil):
    n_faces = len(faces)
    volumes_adj_by_faces = np.repeat(-1, n_faces*2).reshape((n_faces, 2)).astype(np.uint64)
    volumes_series = pd.DataFrame({
        'vol_ids': volumes
    }, index = np.array(all_volumes))

    nodes_series = pd.DataFrame({
        'nodes_ids': nodes
    }, index=all_nodes)

    nodes_of_faces = np.repeat(-1, n_faces*4).reshape((n_faces, 4)).astype(np.uint64)

    for i, face in enumerate(all_faces):
        volumes_adj_by_faces[i][:] = mtu.get_bridge_adjacencies(face, 2, 3)
        nodes_of_faces[i][:] = mtu.get_bridge_adjacencies(face, 2, 0)

    test = volumes_adj_by_faces[:, 0] == volumes_adj_by_faces[:, 1]
    volumes_adj_by_faces[:, 0] = volumes_series.loc[volumes_adj_by_faces[:, 0]].to_numpy().flatten()
    volumes_adj_by_faces[:, 1] = volumes_series.loc[volumes_adj_by_faces[:, 1]].to_numpy().flatten()
    volumes_adj_by_faces = volumes_adj_by_faces.astype(np.int64)
    volumes_adj_by_faces[test, 1] = -1

    for i in range(nodes_of_faces.shape[1]):
        nodes_of_faces[:, i] = nodes_series.loc[nodes_of_faces[:, i]].to_numpy().flatten()
    nodes_of_faces = nodes_of_faces.astype(np.int64)
    bool_internal_faces = test        

    return bool_internal_faces, volumes_adj_by_faces, nodes_of_faces

def get_mesh_properties(mesh_path: str):
    data = dict()
    mb, mtu, root_set = init_mesh(mesh_path)
    all_volumes, all_faces, all_edges, all_nodes = init_mesh_entities(mb, mtu)
    volumes = np.arange(len(all_volumes), dtype=np.uint64)
    faces = np.arange(len(all_faces), dtype=np.uint64)
    edges = np.arange(len(all_edges), dtype=np.uint64)
    nodes = np.arange(len(all_nodes), dtype=np.uint64)
    bool_internal_faces, volumes_adj_by_faces, nodes_of_faces = init_properties(faces, volumes, nodes, all_volumes, all_faces, all_nodes, mtu)
    
    nodes_centroids = np.array([mb.get_coords(node) for node in all_nodes])
    
    data['volumes'] = volumes
    data['faces'] = faces
    data['edges'] = edges
    data['nodes'] = nodes
    data['bool_internal_faces'] = bool_internal_faces
    data['volumes_adj_by_faces'] = volumes_adj_by_faces
    data['nodes_of_faces'] = nodes_of_faces
    data['nodes_centroids'] = nodes_centroids
    data['mesh_name'] = np.array([mesh_name])
    
    return data

@jit(nopython=True)
def get_faces_of_volumes(volumes_adj_by_faces, volumes, faces):
    
    faces_of_volumes = np.zeros((len(volumes), 6), dtype=np.int64)
    faces_of_volumes[:] = -1
    # test = np.full(volumes_adj_by_faces.shape[0], False, dtype=bool)
    iterator = np.array([1], dtype=np.int64)
    # i = iterator[0]
    for iterator[0] in volumes:
        # test[:] = (volumes_adj_by_faces[:, 0] == iterator[0]) | (volumes_adj_by_faces[:, 1] == iterator[0])
        test = (volumes_adj_by_faces[:, 0] == iterator[0]) | (volumes_adj_by_faces[:, 1] == iterator[0])
        faces_of_volumes[iterator[0], :] = faces[test]
        print(iterator[0])
    
    return faces_of_volumes
    
    

class MeshProperties():

    '''
        Create mesh properties using pymoab
    '''

    def __init__(self, mesh_path='', mesh_name=''):
        """MeshProperties creator

        Args:
            mesh_path (str, optional): Path of mesh. Defaults to ''.
            mesh_name (str, optional): Name to save. Defaults to ''.
        """

        self.mesh_path = mesh_path
        self.data = dict()
        self.name = mesh_name

    def create_properties(self):
        self.data.update(get_mesh_properties(self.mesh_path))

    @property
    def volumes(self):
        """The volumes property."""
        return self.data['volumes']

    @property
    def faces(self):
        """The faces property."""
        return self.data['faces']
    
    @property
    def nodes(self):
        """The nodes property."""
        return self.data['nodes']

    @property
    def nodes_of_faces(self):
        """The nodes_of_faces property."""
        return self.data['nodes_of_faces']
    
    @property
    def volumes_adjacencies_by_faces(self):
        """The volumes_adjacencies_by_faces property."""
        return self.data['volumes_adj_by_faces']

    @property
    def nodes_centroids(self):
        """The nodes_centroids property."""
        return self.data['nodes_centroids']

    @property
    def internal_faces(self):
        """The internal_faces property."""
        return self.faces[self.data['bool_internal_faces']]
    
    @property
    def boundary_faces(self):
        """The boundary_faces property."""
        bool_boundary_faces = ~self.data['bool_internal_faces']
        return self.faces[bool_boundary_faces]
    
    @property
    def faces_of_volumes(self):
        """The faces_of_volumes property."""
        try:
            return self.data['faces_of_volumes']
        except KeyError:
            self.data['faces_of_volumes'] = get_faces_of_volumes(
                self.volumes_adjacencies_by_faces,
                self.volumes,
                self.faces
            )
            self.export_data()
            return self.data['faces_of_volumes']
        except Exception as e:
            raise e
            
    def export_data(self):
        manager = DataManager(description=self.name)
        manager._data.update(self.data)
        manager.export_to_npz()

    def load_data(self):
        manager = DataManager(description=self.name)
        manager.load_from_npz()
        self.data.update(manager._data)

if __name__ == '__main__':
    
    mesh_path = 'mesh/80x80x1_ufce.h5m'
    mesh_name = '80x80x1_ufce'
    mesh_properties = MeshProperties(mesh_path=mesh_path, mesh_name=mesh_name)
    mesh_properties.create_properties()
    fvols = mesh_properties.faces_of_volumes
        
    mesh_properties.export_data()
    mesh2 = MeshProperties(mesh_path=mesh_path, mesh_name=mesh_name)
    mesh2.load_data()
    
    assert np.all(mesh_properties.faces == mesh2.faces)
