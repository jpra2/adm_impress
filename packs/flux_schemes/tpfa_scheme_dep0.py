
from ..data_class.data_manage import dataManager
from .. import directories as direc
import numpy as np
import scipy.sparse as sp

class tpfaScheme:

    def __init__(self, M: 'mesh object', data_name: 'nome do arquivo para ser gravado') -> None:

        self.mesh = M
        self.gravity = direc.data_loaded['gravity']
        self.n_nodes = len(M.data.elements_lv0[direc.entities_lv0[0]])
        self.n_faces = len(M.data.elements_lv0[direc.entities_lv0[1]])
        self.n_edges = len(M.data.elements_lv0[direc.entities_lv0[2]])
        self.n_volumes = len(M.data.elements_lv0[direc.entities_lv0[3]])
        self.data = dataManager(data_name)
        M.scheme = self

    def get_transmissibility_matrix_without_boundary_conditions(self) -> None:
        '''
        transmissibility matrix without boundary conditions
        '''

        M = self.mesh
        vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[2]]
        v0 = vols_viz_internal_faces
        internal_faces = M.data.elements_lv0[direc.entities_lv0_0[0]]
        transmissibility_faces = M.data.variables[M.data.variables_impress['transmissibility']]
        transmissibility_internal_faces = transmissibility_faces[internal_faces]
        t0 = transmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, self.n_volumes))

        self.data['Tini'] = T
