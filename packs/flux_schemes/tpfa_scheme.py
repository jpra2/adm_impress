from ..data_class.data_manager import DataManager
from ..directories import data_loaded
import numpy as np
import scipy.sparse as sp

class TpfaScheme(DataManager):

    def __init__(self, data_impress, elements_lv0, data_name: str='TpfaScheme.npz', load=False):
        super().__init__(data_name, load=load)
        self.data_impress = data_impress
        self.elements_lv0 = elements_lv0
        self.n_volumes = self.data_impress.len_entities['volumes']
        self.gravity = data_loaded['gravity']

    def get_transmissibility_matrix_without_boundary_conditions(self) -> None:

        vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        v0 = vols_viz_internal_faces
        internal_faces = self.elements_lv0['internal_faces']
        transmissibility_faces = self.data_impress['transmissibility']
        transmissibility_internal_faces = transmissibility_faces[internal_faces]
        t0 = transmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, self.n_volumes))

        self['Tini'] = T
