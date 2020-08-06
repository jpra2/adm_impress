from packs.data_class import DataManager
from packs.data_class import SparseDataManager
from packs.errors.err import NameLoadError
import numpy as np
import scipy.sparse as sp
import pdb

class AdjacenciesMatrix(SparseDataManager):
    pass

class MeshElements(DataManager):
    pass

class ElementsData:
    _tipos = ['array', 'sparse']

    def __init__(self, load=False):
        self.matrices = AdjacenciesMatrix(load=load)
        self.mesh_elements = MeshElements(load=load)

    def set_mesh_elements(self, volumes=np.array([]), faces=np.array([]), edges=np.array([]), nodes=np.array([]), boundary_faces=np.array([]), boundary_nodes=np.array([])):
        self.mesh_elements['volumes'] = volumes.copy()
        self.mesh_elements['faces'] = faces.copy()
        self.mesh_elements['edges'] = edges.copy()
        self.mesh_elements['nodes'] = nodes.copy()
        self.mesh_elements['boundary_faces'] = boundary_faces.copy()
        self.mesh_elements['boundary_nodes'] = boundary_nodes.copy()

    def create_adj_matrix_volumes_to_faces(self, volumes, faces, faces_volumes):

        lines = []
        cols = []

        for i, f2 in enumerate(faces_volumes):
            lines.append(np.repeat(i, len(f2)))
            cols.append(f2)

        lines = np.concatenate(lines)
        cols = np.concatenate(cols)
        data = np.full(len(lines), True, dtype=bool)

        self.matrices['Matrix_volumes_to_faces'] = sp.csc_matrix((data, (lines, cols)), dtype=bool, shape=(len(volumes), len(faces)))

    def create_adj_matrix_faces_to_edges(self, faces, edges, edges_faces):

        lines = []
        cols = []

        for i, e2 in enumerate(edges_faces):
            lines.append(np.repeat(i, len(e2)))
            cols.append(e2)

        lines = np.concatenate(lines)
        cols = np.concatenate(cols)
        data = np.full(len(lines), True, dtype=bool)

        self.matrices['Matrix_faces_to_edges'] = sp.csc_matrix((data, (lines, cols)), dtype=bool, shape=(len(faces), len(edges)))

    def create_adj_matrix_edges_to_nodes(self, edges, nodes, nodes_edges):

        lines = []
        cols = []

        for i, n2 in enumerate(nodes_edges):
            lines.append(np.repeat(i, len(n2)))
            cols.append(n2)

        lines = np.concatenate(lines)
        cols = np.concatenate(cols)
        data = np.full(len(lines), True, dtype=bool)

        self.matrices['Matrix_edges_to_nodes'] = sp.csc_matrix((data, (lines, cols)), dtype=bool, shape=(len(edges), len(nodes)))

    def insert(self, name, data, tipo):
        '''
            insere alguma informacao adiocional
            tipo: tipo da informacao - "array" ou "sparse"
            data: dado
            name: nome
        '''

        test = self.test_name(name)
        self.test_tipo(tipo)
        if test == 0:
            if tipo == 'array':
                self.mesh_elements[name] = data
            elif tipo == 'sparse':
                self.matrices[name] = data
        else:
            print(f'{name} already exists in {test[1]}')

    def volumes_to_faces(self, volumes):
        volumes2 = self.test_instance(volumes)
        faces = []
        mat2 = self.matrices['Matrix_volumes_to_faces'][volumes2]
        for i in range(mat2.shape[0]):

            faces.append(self.faces[mat2[i].toarray()[0]])

        faces = np.array(faces)
        return faces

    def faces_to_volumes(self, faces):
        faces2 = self.test_instance(faces)
        volumes = []
        mat2 = self.matrices['Matrix_volumes_to_faces'].transpose()[faces2]
        for i in range(mat2.shape[0]):
            volumes.append(self.volumes[mat2[i].toarray()[0]])

        volumes = np.array(volumes)
        return volumes

    def faces_to_edges(self, faces):
        faces2 = self.test_instance(faces)
        edges = []
        mat2 = self.matrices['Matrix_faces_to_edges'][faces2]
        for i in range(mat2.shape[0]):

            edges.append(self.edges[mat2[i].toarray()[0]])

        edges = np.array(edges)
        return edges

    def edges_to_faces(self, edges):
        edges2 = self.test_instance(edges)
        faces = []
        mat2 = self.matrices['Matrix_faces_to_edges'].transpose()[edges2]
        for i in range(mat2.shape[0]):

            faces.append(self.faces[mat2[i].toarray()[0]])

        faces = np.array(faces)
        return faces

    def volumes_to_edges(self, volumes, M=None):
        volumes2 = self.test_instance(volumes)
        mat2 = self.matrices['Matrix_volumes_to_faces']*self.matrices['Matrix_faces_to_edges']
        mat2 = mat2[volumes2]
        edges = []
        for i in range(mat2.shape[0]):
            t = mat2[i].toarray()[0]
            edges.append(self.edges[t])

        edges = np.array(edges)

        return edges

    def edges_to_volumes(self, edges):
        edges2 = self.test_instance(edges)
        mat2 = self.matrices['Matrix_faces_to_edges'].transpose()*self.matrices['Matrix_volumes_to_faces'].transpose()[edges2]
        volumes = []
        for i in range(mat2.shape[0]):

            volumes.append(self.volumes[mat2[i].toarray()[0]])

        volumes = np.array(volumes)
        return volumes

    def edges_to_nodes(self, edges):

        edges2 = self.test_instance(edges)
        nodes = []
        mat2 = self.matrices['Matrix_edges_to_nodes'][edges2]
        for i in range(mat2.shape[0]):

            nodes.append(self.nodes[mat2[i].toarray()[0]])

        nodes = np.array(nodes)
        return nodes

    def nodes_to_edges(self, nodes):

        nodes2 = self.test_instance(nodes)
        edges = []
        mat2 = self.matrices['Matrix_edges_to_nodes'].transpose()[nodes2]
        for i in range(mat2.shape[0]):

            edges.append(self.edges[mat2[i].toarray()[0]])

        edges = np.array(edges)
        return edges

    def faces_to_nodes(self, faces):
        faces2 = self.test_instance(faces)
        mat2 = self.matrices['Matrix_faces_to_edges']*self.matrices['Matrix_edges_to_nodes']
        mat2 = mat2[faces2]
        nodes = []
        for i in range(mat2.shape[0]):

            nodes.append(self.nodes[mat2[i].toarray()[0]])

        nodes = np.array(nodes)
        return nodes

    def nodes_to_faces(self, nodes):
        nodes2 = self.test_instance(nodes)
        mat2 = self.matrices['Matrix_edges_to_nodes'].transpose()*self.matrices['Matrix_faces_to_edges'].transpose()[nodes2]
        faces = []
        for i in range(mat2.shape[0]):

            faces.append(self['faces'][mat2[i].toarray()[0]])

        faces = np.array(faces)
        return faces

    def volumes_to_nodes(self, volumes):
        volumes2 = self.test_instance(volumes)
        mat2 = self.matrices['Matrix_volumes_to_faces']*self.matrices['Matrix_faces_to_edges']*self.matrices['Matrix_edges_to_nodes'][volumes2]
        nodes = []
        for i in range(mat2.shape[0]):

            nodes.append(self['nodes'][mat2[i].toarray()[0]])

        nodes = np.array(nodes)
        return nodes

    def nodes_to_volumes(self, nodes):
        nodes2 = self.test_instance(nodes)
        mat2 = self.matrices['Matrix_edges_to_nodes'].transpose()*self.matrices['Matrix_faces_to_edges'].transpose()*self.matrices['Matrix_volumes_to_faces'].transpose()[nodes2]
        volumes = []
        for i in range(mat2.shape[0]):

            volumes.append(self['volumes'][mat2[i].toarray()[0]])

        volumes = np.array(volumes)
        return volumes

    def test_instance(self, value):
        if isinstance(value, int):
            return [value]
        elif isinstance(value,list) or isinstance(value, set) or isinstance(value,tuple):
            return value
        elif isinstance(value,np.ndarray):
            return value
        else:
            raise ValueError('\ntype not suported\n')

    def test_name(self, name):
        if name in self.mesh_elements.keys():
            raise NameError(f'{name} already exists in mesh_elements')
        elif name in self.matrices.keys():
            raise NameError(f'{name} already exists in matrices')
        else:
            return 0

    def test_tipo(self, tipo):
        if tipo in ElementsData._tipos:
            pass
        else:
            raise NameError(f'tipo nao esta em {Elements._tipos}')

    def get(self, name):
        try:
            return self.mesh_elements[name]
        except KeyError:
            try:
                return self.matrices[name]
            except KeyError:
                raise NameError(f'{name} not in data')

    def export(self):
        self.mesh_elements.export_to_npz()
        self.matrices.export()

    @property
    def volumes(self):
        return self.mesh_elements['volumes']

    @property
    def faces(self):
        return self.mesh_elements['faces']

    @property
    def edges(self):
        return self.mesh_elements['edges']

    @property
    def nodes(self):
        return self.mesh_elements['nodes']

    @property
    def internal_faces(self):
        return np.setdiff1d(self.mesh_elements['faces'], self.mesh_elements['boundary_faces'])

    @property
    def boundary_faces(self):
        return self.mesh_elements['boundary_faces']

    @property
    def boundary_nodes(self):
        return self.mesh_elements['boudary_nodes']

    @property
    def internal_nodes(self):
        return np.setdiff1d(self.mesh_elements['nodes'], self.mesh_elements['boundary_nodes'])
