import scipy.sparse as sp
import numpy as np


class GlobalIMPECPressureSolver:

    @staticmethod
    def mount_transmissibility_no_bc(xkj_internal_faces, Csi_j_internal_faces, mobilities_internal_faces,
                                     Vbulk, porosity, Cf, delta_t, dVtdP, dVtdk,
                                     pretransmissibility_internal_faces, n_volumes, n_components,
                                     internal_faces_adjacencies):
        """

        @param xkj_internal_faces:
        @param Csi_j_internal_faces:
        @param mobilities_internal_faces:
        @param Vbulk: Volume of each mesh volumes
        @param porosity: porosity of volumes
        @param Cf:
        @param delta_t: time step
        @param dVtdP:
        @param dVtdk:
        @param pretransmissibility_internal_faces: static params of face transmissibility
        @param n_volumes: number of volumes
        @param n_components: number of components
        @param internal_faces_adjacencies: volumes adjacencies of internal faces
        @return: T: transmissibility matrix without boundary conditions - T
        """

        v0 = internal_faces_adjacencies

        t0_internal_faces_prod = xkj_internal_faces * Csi_j_internal_faces * mobilities_internal_faces
        t0 = t0_internal_faces_prod.sum(axis=1)
        t0 = t0 * pretransmissibility_internal_faces

        T = sp.lil_matrix((n_volumes, n_volumes))

        for i in range(n_components):
            lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
            cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
            data = np.array([-t0[i, :], -t0[i, :], +t0[i, :], +t0[i, :]]).flatten()

            # Ta = (sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes))).toarray()
            Ta = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, n_volumes))
            # T += Ta * self.dVtk[i,:, np.newaxis]
            T += Ta.multiply(dVtdk[i, :, np.newaxis])

        T *= delta_t

        T.setdiag(T.diagonal().flatten() + (Vbulk * porosity * Cf - dVtdP))

        return T
