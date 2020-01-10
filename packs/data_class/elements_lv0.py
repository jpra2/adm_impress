from .data_manager import DataManager
import numpy as np

class ElementsLv0(DataManager):

    def __init__(self, M, load=False, data_name: str='elementsLv0.npz'):
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
        self._data['all_volumes'] = self.mesh.core.all_volumes
        self._data['all_faces'] = self.mesh.core.all_faces
        self._data['all_edges'] = self.mesh.core.all_edges
        self._data['all_nodes'] = self.mesh.core.all_nodes

        remaped_internal_faces = np.repeat(-1, len(self._data['faces'])).astype(np.dtype(int))
        remaped_boundary_faces = remaped_internal_faces.copy()
        remaped_internal_faces[self._data['internal_faces']] = np.arange(len(self._data['internal_faces']))
        self._data['remaped_internal_faces'] = remaped_internal_faces
        remaped_boundary_faces[self._data['boundary_faces']] = np.arange(len(self._data['boundary_faces']))
        self._data['remaped_boundary_faces'] = remaped_boundary_faces

    def run(self):
        self.load_elements_from_mesh()
