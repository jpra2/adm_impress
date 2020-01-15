
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

    def get_transmissibility_matrix_without_boundary_conditions_compositional(self) -> None:
        #z - molar fraction of the component in each phase - ex: z = [[zc1.liq, zc2.liq],[zc1.vap,zc2.vap]
        molar_fractions = np.zeros([fluid_properties.Nc,2])
        molar_fractions[0] = fluid_properties.x; molar_fractions[1] = fluid_properties.y
        mass_densities = np.array([fluid_properties.rho_L,fluid_properties.rho_V])
        mobilities = relative_permeabilities/phase_viscosities #2 column vectors
        transmissibility = molar_fractions*mass_densities*mobilities*M.data[M.data.variables_impress['pretransmissibility']]
