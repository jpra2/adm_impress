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

        self['volumes'] = self.mesh.volumes.all
        self['faces'] = self.mesh.faces.all
        self['edges'] = self.mesh.edges.all
        self['nodes'] = self.mesh.nodes.all
        self['internal_faces'] = self.mesh.faces.internal
        self['boundary_faces'] = np.setdiff1d(self['faces'], self['internal_faces'])
        self['neig_faces'] = self.mesh.faces.bridge_adjacencies(self['faces'], 2, 3)
        self['neig_internal_faces'] = self.mesh.faces.bridge_adjacencies(self['internal_faces'], 2, 3)
        self['neig_boundary_faces'] = self.mesh.faces.bridge_adjacencies(self['boundary_faces'], 2, 3).flatten()
        self['all_volumes'] = self.mesh.core.all_volumes
        self['all_faces'] = self.mesh.core.all_faces
        self['all_edges'] = self.mesh.core.all_edges
        self['all_nodes'] = self.mesh.core.all_nodes

    def run(self):
        self.load_elements_from_mesh()
