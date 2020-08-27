import numpy as np
import pdb
from packs.tpfa.biphasic import TpfaBiphasicCons
from packs.common_files.common_infos import CommonInfos
from scipy import sparse as sp
from packs.solvers.solvers_scipy.solver_sp import SolverSp
from packs.properties import PhisicalProperties

class ConservationTest:

    def conservation_with_gravity(self,
        volumes, primal_ids, g_source_internal_faces, pressure, global_transmissibility,
        coarse_faces, coarse_intersect_faces, areas, coarse_primal_ids_org, volumes_adj_internal_faces,
        coarse_internal_faces, map_internal_faces, abs_u_normal_faces,
        mobility_w_internal_faces, mobility_o_internal_faces, gravity_vector, keq_faces,
        rho_w, rho_o, hi, g_source_total_volumes, coarse_vertexes, g_source_total_internal_faces, g_velocity_faces):

        all_primal_ids = np.unique(primal_ids)
        cids = coarse_primal_ids_org.copy()
        flux_coarse_volumes = np.zeros(len(all_primal_ids))
        n_volumes = len(volumes)
        local_ids = np.zeros(len(volumes), dtype=int)
        velocity_internal_faces = np.zeros((len(mobility_o_internal_faces), 3))
        flux_internal_faces = np.zeros(len(mobility_o_internal_faces))
        solver = SolverSp()

        for primal_id in all_primal_ids:
            vols = volumes[primal_ids == primal_id]
            set_vols = set(vols)
            n_vols = len(vols)
            local_ids[vols] = np.arange(n_vols)
            g_source_vols = g_source_total_volumes[vols]
            T_local = CommonInfos.get_local_t(global_transmissibility, vols)
            faces_coarse_volume = coarse_faces[cids == primal_id][0]
            intersect_coarse_faces = coarse_intersect_faces[cids == primal_id][0]
            internal_faces_coarse_volume = coarse_internal_faces[cids == primal_id][0]
            vertex = coarse_vertexes[cids == primal_id][0]
            map_intersect_faces = map_internal_faces[intersect_coarse_faces]
            map_coarse_internal_faces = map_internal_faces[internal_faces_coarse_volume]
            volumes_adj_intersect = volumes_adj_internal_faces[map_intersect_faces]

            delta_p = pressure[volumes_adj_intersect]
            mob_w_int_f = mobility_w_internal_faces[map_intersect_faces]
            mob_o_int_f = mobility_o_internal_faces[map_intersect_faces]
            abs_u_normal_intersect_faces = abs_u_normal_faces[intersect_coarse_faces]
            keq_intersect_faces = keq_faces[intersect_coarse_faces]
            areas_intersect_faces = areas[intersect_coarse_faces]
            hi2 = hi[map_intersect_faces]
            ni = len(intersect_coarse_faces)
            adj_intersect_intern_vols = np.zeros(ni, dtype=int)

            for i, vadjs in enumerate(volumes_adj_intersect):
                if set([vadjs[0]]) & set_vols:
                    adj_intersect_intern_vols[i] = vadjs[0]
                else:
                    adj_intersect_intern_vols[i] = vadjs[1]

            vols2 = np.setdiff1d(vols, adj_intersect_intern_vols[i])

            resp, resp2 = self.get_upscale_flux(
                delta_p,
                abs_u_normal_intersect_faces,
                keq_intersect_faces,
                mob_w_int_f,
                mob_o_int_f,
                hi2,
                rho_w,
                rho_o,
                gravity_vector,
                areas_intersect_faces,
                ni
            )

            velocity_internal_faces[map_intersect_faces] = resp2
            flux_internal_faces[map_intersect_faces] = resp


            flux_coarse_volumes[primal_id] = resp.sum()
            #
            ################################################################
            ## neuman problem
            adj_internal_faces_coarse = volumes_adj_internal_faces[map_coarse_internal_faces]
            local_adj_internal = local_ids[adj_internal_faces_coarse]
            g_source_internal = g_source_total_internal_faces[map_coarse_internal_faces]
            g_source_local_volumes = PhisicalProperties.get_total_g_source_volumes(local_ids[vols], local_adj_internal, g_source_internal)
            local_source_term = g_source_local_volumes.copy()
            local_source_term[local_ids[adj_intersect_intern_vols]] += resp
            # local_source_term = g_source_total_volumes[vols]
            # local_source_term[local_ids[vols2]] = g_source_total_volumes[vols2]
            # local_source_term[local_ids[adj_intersect_intern_vols]] = resp
            # local_source_term += g_source_total_volumes[vols]
            local_source_term[local_ids[vertex]] = pressure[vertex]
            T_local[vertex] = 0
            T_local[vertex, vertex] = 1
            p2 = solver.direct_solver(T_local, local_source_term)

            delta_p = p2[local_adj_internal]
            abs_u_normal_internal_faces = abs_u_normal_faces[internal_faces_coarse_volume]
            keq_internal = keq_faces[internal_faces_coarse_volume]
            mob_w_int_f = mobility_w_internal_faces[map_coarse_internal_faces]
            mob_o_int_f = mobility_o_internal_faces[map_coarse_internal_faces]
            hi2 = hi[map_coarse_internal_faces]
            areas_internal = areas[internal_faces_coarse_volume]
            ni_internal = len(internal_faces_coarse_volume)

            local_flux, local_velocity  = self.get_upscale_flux(
                delta_p,
                abs_u_normal_internal_faces,
                keq_internal,
                mob_w_int_f,
                mob_o_int_f,
                hi2,
                rho_w,
                rho_o,
                gravity_vector,
                areas_internal,
                ni_internal
            )

            velocity_internal_faces[map_coarse_internal_faces] = local_velocity
            flux_internal_faces[map_coarse_internal_faces] = local_flux
            ################################################################

        return flux_internal_faces, velocity_internal_faces

    def get_upscale_flux(self, delta_p, abs_u_normal_intersect_faces, keq_intersect_faces,
            mob_w_int_f, mob_o_int_f, hi2, rho_w, rho_o, gravity_vector, areas_intersect_faces, ni):

        delta_p2 = delta_p[:, 1] - delta_p[:, 0]
        delta_p2 = delta_p2.reshape(ni, 1)
        pressure_direction = -delta_p2*abs_u_normal_intersect_faces
        keq = keq_intersect_faces.reshape(ni, 1)
        mob_w_int_f = mob_w_int_f.reshape(ni, 1)
        mob_o_int_f = mob_o_int_f.reshape(ni, 1)
        mob_t = mob_w_int_f + mob_o_int_f
        hi2 = hi2.sum(axis=1).reshape(ni, 1)
        k1 = pressure_direction*(mob_t)*keq
        k2 = (rho_w*mob_w_int_f + rho_o*mob_o_int_f)*gravity_vector*keq*hi2
        # resp = ((k1 + k2)*abs_u_normal_intersect_faces).sum(axis=1)
        resp2 = (k1 + k2)*abs_u_normal_intersect_faces
        resp = resp2.sum(axis=1)*areas_intersect_faces
        return resp, resp2

    def get_local_velocity(self, pressure, volumes, volumes_adj_faces, keq_faces, g_total_velocity_faces):

        pass
