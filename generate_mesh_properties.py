# from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh
from pymoab import core, types, rng, topo_util
from pymoab import skinner as sk
import numpy as np
import pandas as pd
import yaml
import copy

from packs.data_class.data_manager import DataManager

class MeshProperties():

    '''
        Create mesh properties using pymoab
    '''

    def initialize(self, mesh_path='', mesh_name=''):

        self.mesh_path = mesh_path
        self.all_volumes = None
        self.all_faces = None
        self.all_nodes = None
        self.all_edges = None
        self.mb: core.Core = None
        self.mtu: topo_util.MeshTopoUtil = None
        self.root_set = None
        self.data = dict()
        self.mesh_name = ''

    def _init_mesh(self):
        mb = core.Core()
        mtu = topo_util.MeshTopoUtil(mb)
        mb.load_file(self.mesh_path)
        root_set = mb.get_root_set()
        return mb, mtu, root_set

    def init_mesh(self):
        self.mb, self.mtu, self.root_set = self._init_mesh()

    def init_mesh_entities(self):
        self.all_volumes = self.mb.get_entities_by_dimension(0, 3)
        self.all_nodes = self.mb.get_entities_by_dimension(0, 0)
        self.mtu.construct_aentities(self.all_nodes)
        self.all_faces = self.mb.get_entities_by_dimension(0, 2)
        self.all_edges = self.mb.get_entities_by_dimension(0, 1)

    def _init_properties(self, faces, volumes, nodes):
        n_faces = len(faces)
        volumes_adj_by_faces = np.repeat(-1, n_faces*2).reshape((n_faces, 2)).astype(np.uint64)
        volumes_series = pd.DataFrame({
            'vol_ids': volumes
        }, index = np.array(self.all_volumes))

        nodes_series = pd.DataFrame({
            'nodes_ids': nodes
        }, index=self.all_nodes)

        nodes_of_faces = np.repeat(-1, n_faces*4).reshape((n_faces, 4)).astype(np.uint64)

        for i, face in enumerate(self.all_faces):
            volumes_adj_by_faces[i][:] = self.mtu.get_bridge_adjacencies(face, 2, 3)
            nodes_of_faces[i][:] = self.mtu.get_bridge_adjacencies(face, 2, 0)

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

    def create_initial_array_properties(self):

        volumes = np.arange(len(self.all_volumes), dtype=int)
        faces = np.arange(len(self.all_faces), dtype=int)
        edges = np.arange(len(self.all_edges), dtype=int)
        nodes = np.arange(len(self.all_nodes), dtype=int)
        bool_internal_faces, volumes_adj_by_faces, nodes_of_faces =  self._init_properties(faces, volumes, nodes)

        nodes_centroids = np.array([self.mb.get_coords(node) for node in self.all_nodes])

        # av2 = np.array(self.all_volumes, np.uint64)
        #
        # vc = self.mtu.get_average_position(av2[0])

        self.data['volumes'] = volumes
        self.data['faces'] = faces
        self.data['edges'] = edges
        self.data['nodes'] = nodes
        self.data['bool_internal_faces'] = bool_internal_faces
        self.data['volumes_adj_by_faces'] = volumes_adj_by_faces
        self.data['nodes_of_faces'] = nodes_of_faces
        self.data['nodes_centroids'] = nodes_centroids
        self.data['mesh_name'] = np.array([self.mesh_name])

    def volumes():
        doc = "The volumes property."
        def fget(self):
            return self.data['volumes']
        return locals()
    volumes = property(**volumes())

    def faces():
        doc = "The faces property."
        def fget(self):
            return self.data['faces']
    faces = property(**faces())

    def nodes():
        doc = "The nodes property."
        def fget(self):
            return self.data['nodes']
        return locals()
    nodes = property(**nodes())

    def nodes_of_faces():
        doc = "The nodes_of_faces property."
        def fget(self):
            return self.data['nodes_of_faces']
        return locals()
    nodes_of_faces = property(**nodes_of_faces())

    def volumes_adjacencies_by_faces():
        doc = "The volumes_adjacencies_by_faces property."
        def fget(self):
            return self.data['volumes_adj_by_faces']
        return locals()
    volumes_adjacencies_by_faces = property(**volumes_adjacencies_by_faces())

    def nodes_centroids():
        doc = "The nodes_centroids property."
        def fget(self):
            return self.data['nodes_centroids']
        return locals()
    nodes_centroids = property(**nodes_centroids())

    def internal_faces():
        doc = "The internal_faces property."
        def fget(self):
            return self.faces[self.data['bool_internal_faces']]
        return locals()
    internal_faces = property(**internal_faces())

    def export_data(self):
        manager = DataManager(description=self.mesh_name)
        manager._data = self.data
        manager.export_to_npz()

    def load_data(self):
        manager = DataManager(description=self.mesh_name)
        manager.load_data()
        self.data = copy.deepcopy(manager._data)





if __name__ == '__main__':
    mesh_path = 'mesh/500x500_ufce.h5m'
    mesh_name = '500x500_ufce'
    mesh_properties = MeshProperties()
    mesh_properties.initialize(mesh_path=mesh_path, mesh_name=mesh_name)
    mesh_properties.init_mesh()
    mesh_properties.init_mesh_entities()
    mesh_properties.create_initial_array_properties()
