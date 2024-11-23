from .data_manager import DataManager
import numpy as np
from ..directories import only_mesh_name
import pdb
from scipy import sparse as sp


class ElementsLv0(DataManager):

    def __init__(self, M, load=False, data_name: str = 'elementsLv0'):

        data_name = data_name + '_' + only_mesh_name + '.npz'
        super().__init__(data_name, load=load)
        self.mesh = M
        if not load:
            self.run()

        self._loaded = True

    def load_elements_from_mesh(self):
        self._data['volumes'] = self.mesh.volumes.all
        self._data['faces'] = self.mesh.faces.all
        self._data['edges'] = self.mesh.edges.all
        self._data['nodes'] = self.mesh.nodes.all
        self._data['internal_faces'] = self.mesh.faces.internal
        self._data['boundary_faces'] = np.setdiff1d(self['faces'], self['internal_faces'])
        self._data['neig_faces'] = self.mesh.faces.bridge_adjacencies(self['faces'], 2, 3)
        self._data['neig_internal_faces'] = self.mesh.faces.bridge_adjacencies(self['internal_faces'], 2, 3)
        self._data['neig_boundary_faces'] = self.mesh.faces.bridge_adjacencies(self['boundary_faces'], 2, 3).flatten()
        # self._data['all_volumes'] = self.mesh.core.all_volumes
        # self._data['all_faces'] = self.mesh.core.all_faces
        # self._data['all_edges'] = self.mesh.core.all_edges
        # self._data['all_nodes'] = self.mesh.core.all_nodes

        remaped_internal_faces = np.repeat(-1, len(self._data['faces'])).astype(np.dtype(int))
        # remaped_boundary_faces = remaped_internal_faces.copy()
        remaped_internal_faces[self._data['internal_faces']] = np.arange(len(self._data['internal_faces']))
        self._data['remaped_internal_faces'] = remaped_internal_faces
        # remaped_boundary_faces[self._data['boundary_faces']] = np.arange(len(self._data['boundary_faces']))
        # self._data['remaped_boundary_faces'] = remaped_boundary_faces
        #
        self._data['volumes_face_faces'] = self.mesh.volumes.bridge_adjacencies(self._data['volumes'], 2, 2)
        # self._data['volumes_face_nodes'] = self.mesh.volumes.bridge_adjacencies(self._data['volumes'], 2, 0)
        self._data['faces_face_volumes'] = self.mesh.faces.bridge_adjacencies(self._data['faces'], 2, 3)
        self._data['volumes_face_volumes'] = self.mesh.volumes.bridge_adjacencies(self._data['volumes'], 2, 3)
        # self._data['faces_edge_edges'] = self.mesh.faces.bridge_adjacencies(self._data['faces'], 1, 1)
        # self._data['faces_node_nodes'] = self.mesh.faces.bridge_adjacencies(self._data['faces'], 0, 0)
        # self._data['edges_node_nodes'] = self.mesh.edges.bridge_adjacencies(self._data['edges'], 0, 0)
        # self._data['nodes_node_faces'] = self.mesh.nodes.bridge_adjacencies(self._data['nodes'], 0, 2)

    def run(self):
        self.load_elements_from_mesh()

    def nfaces():
        doc = "The nfaces property."

        def fget(self):
            try:
                return self._nfaces
            except AttributeError:
                nfaces = len(self._data['faces'])
                self._nfaces = nfaces
                return nfaces

        return locals()

    nfaces = property(**nfaces())

    def create_adj_matrix_volumes_to_faces(self):
        faces_volumes = self.mesh.volumes.bridge_adjacencies(self._data['volumes'], 2, 2)
        lines = []
        cols = []

        for i, faces in enumerate(faces_volumes):
            lines.append(np.repeat(i, len(faces)))
            cols.append(faces)

        lines = np.concatenate(lines)
        cols = np.concatenate(cols)
        data = np.full(len(lines), True, dtype=bool)

        self._adj_matrix_volumes_to_faces = sp.csc_matrix((data, (lines, cols)), dtype=bool,
                                                          shape=(len(self['volumes']), len(self['faces'])))

    def create_adj_matrix_faces_to_edges(self):

        edges_faces = self.mesh.faces.bridge_adjacencies(self._data['faces'], 2, 1)
        lines = []
        cols = []

        for i, edges in enumerate(edges_faces):
            lines.append(np.repeat(i, len(edges)))
            cols.append(edges)

        lines = np.concatenate(lines)
        cols = np.concatenate(cols)
        data = np.full(len(lines), True, dtype=bool)

        self._adj_matrix_faces_to_edges = sp.csc_matrix((data, (lines, cols)), dtype=bool,
                                                        shape=(len(self['faces']), len(self['edges'])))

    def create_adj_matrix_edges_to_nodes(self):

        nodes_edges = self.mesh.edges.bridge_adjacencies(self._data['edges'], 1, 0)
        lines = []
        cols = []

        for i, nodes in enumerate(nodes_edges):
            lines.append(np.repeat(i, len(nodes)))
            cols.append(nodes)

        lines = np.concatenate(lines)
        cols = np.concatenate(cols)
        data = np.full(len(lines), True, dtype=bool)

        self._adj_matrix_edges_to_nodes = sp.csc_matrix((data, (lines, cols)), dtype=bool,
                                                        shape=(len(self['edges']), len(self['nodes'])))

    def volumes_to_faces(self, volumes):
        volumes2 = self.test_instance(volumes)
        faces = []
        mat2 = self.adj_matrix_volumes_to_faces[volumes2]
        for i in range(mat2.shape[0]):
            faces.append(self['faces'][mat2[i].toarray()[0]])

        faces = np.array(faces)
        return faces

    def faces_to_volumes(self, faces):
        faces2 = self.test_instance(faces)
        volumes = []
        mat2 = self.adj_matrix_volumes_to_faces.transpose()[faces2]
        for i in range(mat2.shape[0]):
            volumes.append(self['volumes'][mat2[i].toarray()[0]])

        volumes = np.array(volumes)
        return volumes

    def faces_to_edges(self, faces):
        faces2 = self.test_instance(faces)
        edges = []
        mat2 = self.adj_matrix_faces_to_edges[faces2]
        for i in range(mat2.shape[0]):
            edges.append(self['edges'][mat2[i].toarray()[0]])

        edges = np.array(edges)
        return edges

    def edges_to_faces(self, edges):
        edges2 = self.test_instance(edges)
        faces = []
        mat2 = self.adj_matrix_faces_to_edges.transpose()[edges2]
        for i in range(mat2.shape[0]):
            faces.append(self['faces'][mat2[i].toarray()[0]])

        faces = np.array(faces)
        return faces

    def volumes_to_edges(self, volumes):
        volumes2 = self.test_instance(volumes)
        mat2 = self.adj_matrix_volumes_to_faces * self.adj_matrix_faces_to_edges[volumes2]
        nodes = []
        for i in range(mat2.shape[0]):
            nodes.append(self['nodes'][mat2[i].toarray()[0]])

        nodes = np.array(nodes)
        return nodes

    def edges_to_volumes(self, edges):
        edges2 = self.test_instance(edges)
        mat2 = self.adj_matrix_faces_to_edges.transpose() * self.adj_matrix_volumes_to_faces.transpose()[edges2]
        volumes = []
        for i in range(mat2.shape[0]):
            volumes.append(self['volumes'][mat2[i].toarray()[0]])

        volumes = np.array(volumes)
        return volumes

    def edges_to_nodes(self, edges):

        edges2 = self.test_instance(edges)
        nodes = []
        mat2 = self.adj_matrix_edges_to_nodes[edges2]
        for i in range(mat2.shape[0]):
            nodes.append(self['nodes'][mat2[i].toarray()[0]])

        nodes = np.array(nodes)
        return nodes

    def nodes_to_edges(self, nodes):

        nodes2 = self.test_instance(nodes)
        edges = []
        mat2 = self.adj_matrix_edges_to_nodes.transpose()[nodes2]
        for i in range(mat2.shape[0]):
            edges.append(self['edges'][mat2[i].toarray()[0]])

        edges = np.array(edges)
        return edges

    def faces_to_nodes(self, faces):
        faces2 = self.test_instance(faces)
        mat2 = self.adj_matrix_faces_to_edges * self.adj_matrix_edges_to_nodes[faces2]
        nodes = []
        for i in range(mat2.shape[0]):
            nodes.append(self['nodes'][mat2[i].toarray()[0]])

        nodes = np.array(nodes)
        return nodes

    def nodes_to_faces(self, nodes):
        nodes2 = self.test_instance(nodes)
        mat2 = self.adj_matrix_edges_to_nodes.transpose() * self.adj_matrix_faces_to_edges.transpose()[nodes2]
        faces = []
        for i in range(mat2.shape[0]):
            faces.append(self['faces'][mat2[i].toarray()[0]])

        faces = np.array(faces)
        return faces

    def volumes_to_nodes(self, volumes):
        volumes2 = self.test_instance(volumes)
        mat2 = self.adj_matrix_volumes_to_faces * self.adj_matrix_faces_to_edges * self.adj_matrix_edges_to_nodes[
            volumes2]
        nodes = []
        for i in range(mat2.shape[0]):
            nodes.append(self['nodes'][mat2[i].toarray()[0]])

        nodes = np.array(nodes)
        return nodes

    def nodes_to_volumes(self, nodes):
        nodes2 = self.test_instance(nodes)
        mat2 = self.adj_matrix_edges_to_nodes.transpose() * self.adj_matrix_faces_to_edges.transpose() * \
               self.adj_matrix_volumes_to_faces.transpose()[nodes2]
        volumes = []
        for i in range(mat2.shape[0]):
            volumes.append(self['volumes'][mat2[i].toarray()[0]])

        volumes = np.array(volumes)
        return volumes

    def test_instance(self, value):
        if isinstance(value, int):
            return [value]
        elif isinstance(value, list) or isinstance(value, set) or isinstance(value, tuple):
            return value
        elif isinstance(value, np.ndarray):
            return value
        else:
            raise ValueError('\ntype not suported\n')

    @property
    def adj_matrix_volumes_to_faces(self):
        try:
            return self._adj_matrix_volumes_to_faces
        except AttributeError:
            self.create_adj_matrix_volumes_to_faces()
            return self._adj_matrix_volumes_to_faces

    @property
    def adj_matrix_faces_to_edges(self):
        try:
            return self._adj_matrix_faces_to_edges
        except AttributeError:
            self.create_adj_matrix_faces_to_edges()
            return self._adj_matrix_faces_to_edges

    @property
    def adj_matrix_edges_to_nodes(self):
        try:
            return self._adj_matrix_edges_to_nodes
        except AttributeError:
            self.create_adj_matrix_edges_to_nodes()
            return self._adj_matrix_edges_to_nodes
