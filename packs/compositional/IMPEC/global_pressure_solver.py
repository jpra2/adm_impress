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
    def mount_independent_term(Vbulk, porosity, Cf, dVtdP, P, n_volumes, n_components, n_phases, internal_faces_adjacencies, dVtdk, z_centroids, xkj_internal_faces, Csi_j_internal_faces, mobilities_internal_faces, pretransmissibility_internal_faces, Pcap, Vp, Vt, well_volumes_flux_prescription, values_flux_prescription, delta_t, g, well_volumes_pressure_prescription, pressure_prescription, bhp_ind, rho_j, rho_j_internal_faces):
        """

        @param Vbulk: Volume of each mesh volumes (n_volumes)
        @param porosity: porosity of volumes (n_volumes)
        @param Cf: rock compressibility (cte)
        @param dVtdP: volume derivatives with respect to pressure (n_volumes)
        @param P: pressure field
        @param n_volumes: number of volumes
        @param n_components: number of components
        @param n_phases: number of phases
        @param internal_faces_adjacencies: volumes adjacencies of internal faces
        @param dVtdk: volume derivatives with respect to component mol number (n_components, n_volumes)
        @param z_centroids:  negative z centroids of volumes = -z_mesh_centroids
        @param xkj_internal_faces: concentration of component in phase (n_components, n_phases, n_internal_faces)
        @param Csi_j_internal_faces: molar densities of phases (1, n_phases, n_internal_faces)
        @param mobilities_internal_faces: (1, n_phases, n_internal_faces)
        @param pretransmissibility_internal_faces: static params of face transmissibility
        @param Pcap: capipllary pressure of phases (n_phases, n_volumes)
        @param Vp: porous volume (n_volumes)
        @param Vt: fluid volume (n_volumes)
        @param well_volumes_flux_prescription: gids of volumes with flux prescription
        @param values_flux_prescription: values of flux prescription
        @param delta_t: time step
        @param g: gravity acc in z direction
        @param well_volumes_pressure_prescription: gids of volumes with pressure prescription
        @param pressure_prescription: pressure prescription
        @param bhp_ind: z minimum
        @param rho_j: density of phases (n_phases, n_volumes)
        @param rho_j_internal_faces: density of phases in internal_faces (1, n_phases, n_volumes)
        @return: independent terms
        """
        # import pdb; pdb.set_trace()
        pressure_term = GlobalIMPECPressureSolver.pressure_independent_term(Vbulk, porosity, Cf, dVtdP, P)
        gravity_term = GlobalIMPECPressureSolver.gravity_independent_term(xkj_internal_faces, Csi_j_internal_faces, mobilities_internal_faces, rho_j_internal_faces, pretransmissibility_internal_faces, n_volumes, n_components, internal_faces_adjacencies, dVtdk, z_centroids, g)
        capillary_term = GlobalIMPECPressureSolver.capillary_independent_term(xkj_internal_faces, Csi_j_internal_faces, mobilities_internal_faces, pretransmissibility_internal_faces, n_components, n_phases, n_volumes, dVtdk, Pcap, internal_faces_adjacencies)
        volume_term = GlobalIMPECPressureSolver.volume_discrepancy_independent_term(Vp, Vt)
        well_term = GlobalIMPECPressureSolver.well_term(well_volumes_flux_prescription, values_flux_prescription, n_components, n_volumes, dVtdk)
        independent_terms = pressure_term - volume_term + delta_t * well_term - delta_t * (capillary_term + gravity_term)
        if len(well_volumes_pressure_prescription) > 0:
            independent_terms[well_volumes_pressure_prescription] = pressure_prescription + g * rho_j[0, 0, well_volumes_pressure_prescription] * (z_centroids[well_volumes_pressure_prescription] - z_centroids[bhp_ind])
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

    @staticmethod
    def well_term(well_volumes_flux_prescription, values_flux_prescription, n_components, n_volumes, dVtdk):
        """

        @param well_volumes_flux_prescription: gids of volumes with flux prescription
        @param values_flux_prescription: values of flux prescription
        @param n_components: number of components
        @param n_volumes: number of volumes
        @param dVtdk: volume derivatives with respect to component mol number (n_components, n_volumes)
        @return: well term
        """
        _q = np.zeros([n_components, n_volumes])
        well_term = np.zeros(n_volumes)

        if len(well_volumes_flux_prescription) > 0:
            _q[:, well_volumes_flux_prescription] = values_flux_prescription
            well_term[well_volumes_flux_prescription] = np.sum(dVtdk[:, well_volumes_flux_prescription] *
                                              _q[:, well_volumes_flux_prescription], axis=0)
        return well_term

    @staticmethod
    def update_flux(Ft_internal_faces, rho_j_internal_faces, mobilities_internal_faces, Pcap, adjacencies_internal_faces, z_centroids, pretransmissibility_internal_faces, g, xkj_internal_faces, Csi_j_internal_faces, n_components, n_volumes):
        ''' Main function that calls others
        @param Ft_internal_faces: total flux internal faces
        @param rho_j_internal_faces: phase density internal faces
        @param mobilities_internal_faces: mobilities internal faces
        @param Pcap: capillary pressure
        @param adjacencies_internal_faces: volumes adjacencies of internal faces
        @param z_centroids: volumes centroids
        @param pretransmissibility_internal_faces:
        @param g: gravity acc
        @param xkj_internal_faces: composition in phase internal faces
        @param Csi_j_internal_faces: molar density internal faces
        @param n_components: number of components
        @param n_volumes: number of volumes
        @return: molar balance in volumes
        '''
        v0 = adjacencies_internal_faces
        Fj_internal_faces = GlobalIMPECPressureSolver.update_Fj_internal_faces(
            Ft_internal_faces,
            rho_j_internal_faces,
            mobilities_internal_faces,
            Pcap[:, v0],
            z_centroids[v0],
            pretransmissibility_internal_faces,
            g
        )

        Fk_internal_faces = GlobalIMPECPressureSolver.update_Fk_internal_faces(
            xkj_internal_faces,
            Csi_j_internal_faces,
            Fj_internal_faces
        )

        Fk_vols_total = GlobalIMPECPressureSolver.update_flux_volumes(
            Fk_internal_faces,
            n_components,
            v0,
            n_volumes
        )

        return Fk_vols_total

    @staticmethod
    def update_Fj_internal_faces(Ft_internal_faces, rho_j_internal_faces,
        mobilities_internal_faces, Pcap_face, z_face,
        pretransmissibility_internal_faces, g):
        ''' Function to calculate phase flux '''

        frj = mobilities_internal_faces[0,...] / \
        np.sum(mobilities_internal_faces[0,...], axis = 0)

        Fj_internal_faces = frj[np.newaxis,...] * (Ft_internal_faces +
            pretransmissibility_internal_faces * (np.sum(mobilities_internal_faces *
            (Pcap_face[:,:,1] - Pcap_face[:,:,0] - g * rho_j_internal_faces *
            (z_face[:,1] - z_face[:,0])), axis=1) - np.sum(mobilities_internal_faces,
            axis=1) * (Pcap_face[:,:,1] - Pcap_face[:,:,0] - g *
            rho_j_internal_faces * (z_face[:,1] - z_face[:,0]))))

        return Fj_internal_faces
        # M.flux_faces[M.faces.internal] = Ft_internal_faces * M.faces.normal[M.faces.internal].T

    @staticmethod
    def update_Fk_internal_faces(xkj_internal_faces, Csi_j_internal_faces, Fj_internal_faces):
        ''' Function to compute component flux '''

        Fk_internal_faces = np.sum(xkj_internal_faces * Csi_j_internal_faces *
        Fj_internal_faces, axis = 1)
        return Fk_internal_faces

    @staticmethod
    def update_flux_volumes(Fk_internal_faces, n_components, adjacencies_internal_faces, n_volumes):
        ''' Function to compute component flux balance through the control \
        volume interfaces'''
        v0 = adjacencies_internal_faces

        cx = np.arange(n_components)
        lines = np.array([np.repeat(cx,len(v0[:,0])), np.repeat(cx,len(v0[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(v0[:,0], n_components), np.tile(v0[:,1], n_components)]).flatten()
        data = np.array([-Fk_internal_faces, Fk_internal_faces]).flatten()
        # fprop.Fk_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (n_components, n_volumes)).toarray()
        Fk_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (n_components, n_volumes)).toarray()
        return Fk_vols_total