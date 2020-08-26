import numpy as np
import pdb
from packs.tpfa.biphasic import TpfaBiphasicCons

class ConservationTest:

    def conservation_with_gravity(self,
        volumes, primal_ids, g_source_internal_faces, pressure, global_transmissibility,
        coarse_faces, coarse_intersect_faces, areas, coarse_primal_ids_org, volumes_adj_internal_faces,
        coarse_internal_faces, map_internal_faces, abs_u_normal_faces,
        mobility_w_internal_faces, mobility_o_internal_faces, gravity_vector, keq_faces,
        rho_w, rho_o, hi):


        all_primal_ids = np.unique(primal_ids)
        cids = coarse_primal_ids_org.copy()
        flux_coarse_volumes = np.zeros(len(all_primal_ids))

        for primal_id in all_primal_ids:
            vols = volumes[primal_ids == primal_id]
            faces_coarse_volume = coarse_faces[cids == primal_id][0]
            intersect_coarse_faces = coarse_intersect_faces[cids == primal_id][0]
            internal_faces_coarse_volume = coarse_internal_faces[cids == primal_id][0]
            volumes_adj_intersect = volumes_adj_internal_faces[map_internal_faces[intersect_coarse_faces]]
            delta_p = pressure[volumes_adj_intersect]
            ni = len(intersect_coarse_faces)
            abs_u_normal_intersect_faces = abs_u_normal_faces[intersect_coarse_faces]

            delta_p = delta_p[:, 1] - delta_p[:, 0]
            delta_p = delta_p.reshape(ni, 1)
            pressure_direction = -delta_p*abs_u_normal_intersect_faces
            mob_w_int_f = mobility_w_internal_faces[map_internal_faces[intersect_coarse_faces]]
            mob_o_int_f = mobility_o_internal_faces[map_internal_faces[intersect_coarse_faces]]
            keq = keq_faces[intersect_coarse_faces]
            mob_w_int_f = mob_w_int_f.reshape(ni, 1)
            mob_o_int_f = mob_o_int_f.reshape(ni, 1)
            keq = keq.reshape(ni, 1)
            mob_t = mob_w_int_f + mob_o_int_f
            hi2 = hi[map_internal_faces[intersect_coarse_faces]]
            hi2 = hi2.sum(axis=1).reshape(ni, 1)

            k1 = pressure_direction*(mob_t)*keq
            k2 = (rho_w*mob_w_int_f + rho_o*mob_o_int_f)*gravity_vector*keq*hi2
            resp = ((k1 + k2)*abs_u_normal_intersect_faces).sum(axis=1)
            resp = resp*areas[intersect_coarse_faces]
            flux_coarse_volumes[primal_id] = resp.sum()

        pdb.set_trace()
