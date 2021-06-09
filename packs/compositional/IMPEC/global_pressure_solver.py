import scipy.sparse as sp
import numpy as np


class GlobalIMPECPressureSolver:

    @staticmethod
    def mount_transmissibility_no_bc(xkj_internal_faces, Csi_j_internal_faces, mobilities_internal_faces,
                                     Vbulk, porosity, Cf, delta_t, dVtdP, dVtdk,
                                     pretransmissibility_internal_faces, n_volumes, n_components,
                                     internal_faces_adjacencies) -> sp.csc_matrix:
        """

        @param xkj_internal_faces: concentration of component in phase (n_components, n_phases, n_internal_faces)
        @param Csi_j_internal_faces: molar densities of phases (1, n_phases, n_internal_faces)
        @param mobilities_internal_faces: (1, n_phases, n_internal_faces)
        @param Vbulk: Volume of each mesh volumes (n_volumes)
        @param porosity: porosity of volumes (n_volumes)
        @param Cf: rock compressibility (cte)
        @param delta_t: time step
        @param dVtdP: volume derivatives with respect to pressure (n_volumes)
        @param dVtdk: volume derivatives with respect to component mol number (n_components, n_volumes)
        @param pretransmissibility_internal_faces: static params of face transmissibility
        @param n_volumes: number of volumes
        @param n_components: number of hydrocarbon components
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

        return T.tocsc()

    @staticmethod
    def mount_independent_term(Vbulk, porosity, Cf, dVtdP, P, n_volumes, n_components, n_phases, internal_faces_adjacencies, dVtdk, z_centroids, xkj_internal_faces, Csi_j_internal_faces, mobilities_internal_faces, pretransmissibility_internal_faces, Pcap, Vp, Vt):
        pressure_term = GlobalIMPECPressureSolver.pressure_independent_term(Vbulk, porosity, Cf, dVtdP, P)
        gravity_term = GlobalIMPECPressureSolver.gravity_independent_term(n_volumes, n_components, internal_faces_adjacencies, dVtdk, z_centroids)
        capillary_term = GlobalIMPECPressureSolver.capillary_independent_term(xkj_internal_faces, Csi_j_internal_faces, mobilities_internal_faces, pretransmissibility_internal_faces, n_components, n_phases, n_volumes, dVtdk, Pcap)
        volume_term = GlobalIMPECPressureSolver.volume_discrepancy_independent_term(Vp, Vt)


    def update_independent_terms(self, M, fprop, wells, delta_t):
        self.pressure_term = self.pressure_independent_term(fprop)
        self.capillary_term, self.gravity_term = self.capillary_and_gravity_independent_term(fprop)
        self.volume_term = self.volume_discrepancy_independent_term(fprop)
        well_term = self.well_term(fprop, wells)
        independent_terms = self.pressure_term - self.volume_term  + delta_t * \
        well_term - delta_t * (self.capillary_term + self.gravity_term)
        independent_terms[wells['ws_p']] = wells['values_p'] + ctes.g * \
        fprop.rho_j[0,0,wells['ws_p']] * (ctes.z[wells['ws_p']] - ctes.z[ctes.bhp_ind])
        return independent_terms

    @staticmethod
    def pressure_independent_term(Vbulk, porosity, Cf, dVtdP, P):
        """
        @param Vbulk: Volume of each mesh volumes (n_volumes)
        @param porosity: porosity of volumes (n_volumes)
        @param Cf: rock compressibility (cte)
        @param dVtdP: volume derivatives with respect to pressure (n_volumes)
        @param P: pressure of volumes
        @return: pressure independent term
        """
        pressure_term = Vbulk * porosity * Cf - dVtdP
        pressure_term *= P
        return pressure_term

    @staticmethod
    def gravity_independent_term(xkj_internal_faces, Csi_j_internal_faces, mobilities_internal_faces, rho_j_internal_faces, pretransmissibility_internal_faces, n_volumes, n_components, internal_faces_adjacencies, dVtdk, z_centroids, g):
        """

        @param g: gravity acc in z direction (scalar)
        @param xkj_internal_faces: concentration of component in phase (n_components, n_phases, n_internal_faces)
        @param Csi_j_internal_faces: molar densities of phases (1, n_phases, n_internal_faces)
        @param mobilities_internal_faces: (1, n_phases, n_internal_faces)
        @param rho_j_internal_faces: density of phase in internal faces (1, n_phases, n_internal_faces)
        @param pretransmissibility_internal_faces: static params of face transmissibility
        @param n_volumes: number of volumes
        @param n_components: number of components
        @param internal_faces_adjacencies: volumes adjacencies of internal faces
        @param dVtdk: volume derivatives with respect to component mol number (n_components, n_volumes)
        @param z_centroids: centroids of volumes
        @return: gravity independent term
        """

        v0 = internal_faces_adjacencies
        t0_internal_faces_prod = xkj_internal_faces * Csi_j_internal_faces * mobilities_internal_faces
        t0_j = t0_internal_faces_prod * pretransmissibility_internal_faces
        t0_k = g * np.sum(rho_j_internal_faces * t0_j, axis=1)

        grav = np.zeros([n_volumes, n_volumes])

        for i in range(n_components):
            lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
            cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
            data = np.array([t0_k[i, :], t0_k[i, :], -t0_k[i, :], -t0_k[i, :]]).flatten()
            t0_rho = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, n_volumes)).toarray()
            grav += t0_rho * dVtdk[i, :]

        gravity_term = grav @ z_centroids
        return gravity_term

    @staticmethod
    def capillary_independent_term(xkj_internal_faces, Csi_j_internal_faces, mobilities_internal_faces, pretransmissibility_internal_faces, n_components, n_phases, n_volumes, dVtdk, Pcap, internal_faces_adjacencies):
        """

        @param xkj_internal_faces: concentration of component in phase (n_components, n_phases, n_internal_faces)
        @param Csi_j_internal_faces: molar densities of phases (1, n_phases, n_internal_faces)
        @param mobilities_internal_faces: (1, n_phases, n_internal_faces)
        @param pretransmissibility_internal_faces: static params of face transmissibility
        @param n_components: number of components
        @param n_phases: number of phases
        @param n_volumes: number of volumes
        @param dVtdk: volume derivatives with respect to component mol number (n_components, n_volumes)
        @param Pcap: capipllary pressure of phases (n_phases, n_volumes)
        @param internal_faces_adjacencies: volumes adjacencies of internal faces
        @return: capillary independent term
        """

        t0_internal_faces_prod = xkj_internal_faces * Csi_j_internal_faces * mobilities_internal_faces
        t0_j = t0_internal_faces_prod * pretransmissibility_internal_faces
        v0 = internal_faces_adjacencies

        cap = np.zeros([n_volumes])
        for i in range(n_components):
            for j in range(n_phases):
                lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
                cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
                data = np.array([t0_j[i, j, :], t0_j[i, j, :], -t0_j[i, j, :], -t0_j[i, j, :]]).flatten()
                t0 = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, n_volumes)) * dVtdk[i, :]
                cap += t0 @ Pcap[j, :]

        return cap

    @staticmethod
    def volume_discrepancy_independent_term(Vp, Vt):
        """

        @param Vp: porous volume (n_volumes)
        @param Vt: fluid volume (n_volumes)
        @return: volume discrepancy term (n_volumes)
        """
        volume_discrepancy_term = Vp - Vt
        return volume_discrepancy_term


