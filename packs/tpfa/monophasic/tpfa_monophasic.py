from packs.tpfa.monophasic.monophasic_fluid_properties import MonophasicFluidProperties
import scipy.sparse as sp
import numpy as np



class TpfaMonophasic:

    def mount_transmissibility_matrix(self, transmissibility_internal_faces, internal_faces, volumes_adj_internal_faces, volumes):
        n_volumes = len(volumes)

        v0 = volumes_adj_internal_faces
        t0 = transmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, n_volumes))

        return T.tolil()

    def get_linear_problem(self, volumes_dirichlet, volumes_neumman, values_dirichlet, values_neumman, gravity_source_term_volumes, transmissibility_matrix_without_boundary_conditions):

        b = gravity_source_term_volumes.copy()
        b[volumes_dirichlet] = values_dirichlet
        b[volumes_neumman] += values_neumman

        T2 = transmissibility_matrix_without_boundary_conditions.copy()
        T2 [volumes_dirichlet] = 0
        T2[volumes_dirichlet, volumes_dirichlet] = 1.0

        return T2.tocsc(), b

    def get_transmissibility_faces(self, areas_faces, internal_faces, boundary_faces, volumes_adj_internal_faces, volumes_adj_boundary_faces, mi, keq_faces, hi):

        transmissibility = np.zeros(len(areas_faces))
        transmissibility[internal_faces] = (1/mi)*keq_faces[internal_faces]/hi.sum(axis=1)
        # transmissibility[boundary_faces] = keq_faces[boundary_faces]
        # TODO: atualizar essa funcao para o caso de precisar das transmissibilidades nas faces do contorno
        transmissibility[boundary_faces] = np.zeros(len(boundary_faces))

        return transmissibility

    @staticmethod
    def get_total_flux_volumes(flux_internal_faces, volumes, volumes_adj_internal_faces):

        v0 = volumes_adj_internal_faces
        n_volumes = len(volumes)
        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_internal_faces, -flux_internal_faces]).flatten()
        flux_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
        return flux_volumes
