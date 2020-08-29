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
        rho_w, rho_o, hi, g_source_total_volumes, coarse_vertexes, g_source_total_internal_faces, g_velocity_internal_faces,
        wells_p=[], vals_p=[], wells_q=[], vals_q=[]):

        all_primal_ids = np.unique(primal_ids)
        cids = coarse_primal_ids_org.copy()
        flux_coarse_volumes = np.zeros(len(all_primal_ids))
        test_vector = flux_coarse_volumes.copy()
        n_volumes = len(volumes)
        local_ids = np.zeros(len(volumes), dtype=int)
        test_conserv = local_ids.astype(float)
        local_pressure = test_conserv.copy()
        n_internal_faces = len(mobility_o_internal_faces)
        velocity_internal_faces = np.zeros((n_internal_faces, 3))
        flux_internal_faces = np.zeros(n_internal_faces)
        solver = SolverSp()
        transmissibility_internal_faces = np.absolute(global_transmissibility[volumes_adj_internal_faces[:,0], volumes_adj_internal_faces[:,1]].toarray().flatten())
        set_p = set(wells_p)
        set_q = set(wells_q)
        if len(set_p) > 0:
            presc_p = np.zeros(n_volumes)
            presc_p[wells_p] = vals_p

        if len(wells_q) > 0:
            presc_q = np.zeros(n_volumes)
            presc_q[wells_q] = vals_q

        contador = 0

        for primal_id in all_primal_ids:

            flux2 = np.zeros(n_internal_faces)
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
            g_velocity_intersect_faces = g_velocity_internal_faces[map_intersect_faces]
            ni = len(intersect_coarse_faces)
            adj_intersect_intern_vols = np.zeros(ni, dtype=int)
            transmissibility_intersect = transmissibility_internal_faces[map_intersect_faces]
            # pdb.set_trace()

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
                ni,
                g_velocity_intersect_faces,
                transmissibility_intersect
            )

            velocity_internal_faces[map_intersect_faces] = resp2
            flux_internal_faces[map_intersect_faces] = resp
            flux2[map_intersect_faces] = resp

            flux_coarse_volumes[primal_id] = resp.sum()
            #
            # ################################################################
            # neuman problemg_velocity_intersect_faces
            adj_internal_faces_coarse = volumes_adj_internal_faces[map_coarse_internal_faces]
            local_adj_internal = local_ids[adj_internal_faces_coarse]
            presc_flux_vols = PhisicalProperties.get_total_g_source_volumes(volumes, volumes_adj_internal_faces, flux2)[vols]
            g_source_internal = g_source_total_internal_faces[map_coarse_internal_faces]
            g_source_local_volumes = PhisicalProperties.get_total_g_source_volumes(local_ids[vols], local_adj_internal, g_source_internal)
            local_source_term = g_source_local_volumes + presc_flux_vols

            if set_p & set(vols):
                volsp = np.intersect1d(vols, wells_p)
                T_local[local_ids[volsp]] = 0
                T_local[local_ids[volsp], local_ids[volsp]] = 1
                local_source_term[local_ids[volsp]] = presc_p[volsp]

            if set_q & set(wells_q):
                volsq = np.intersect1d(vols, wells_q)
                volsq = np.setdiff1d(volsq, [vertex])
                local_source_term[local_ids[volsq]] += presc_q[volsq]

            # local_source_term = g_source_local_volumes.copy()
            # local_source_term[local_ids[adj_intersect_intern_vols]] += resp
            # local_source_term[local_ids[adj_intersect_intern_vols]] -= resp
            # local_source_term = g_source_total_volumes[vols]
            # local_source_term[local_ids[vols2]] = g_source_total_volumes[vols2]
            # local_source_term[local_ids[adj_intersect_intern_vols]] = resp
            # local_source_term += g_source_total_volumes[vols]
            local_source_term[local_ids[vertex]] = pressure[vertex]
            T_local[local_ids[vertex]] = 0
            T_local[local_ids[vertex], local_ids[vertex]] = 1
            # if contador == 0:
            #     # pvs = pressure[vols]
            #     # ids = local_ids[vols]
            #     # ids = ids[pvs == 1]
            #     # T_local[ids] = 0
            #     # T_local[ids, ids] = 1
            #     # local_source_term[ids] = pvs[ids]
            #     local_source_term[np.absolute(local_source_term) < 1e-4] = 0
            p2 = solver.direct_solver(T_local, local_source_term)
            local_pressure[vols] = local_source_term
            pvols = pressure[vols]
            test = np.allclose(p2, pvols)
            print(test)
            if not test:
                test_vector[primal_id] = 1.0
            ################################################################

            delta_p2 = p2[local_adj_internal]
            abs_u_normal_internal_faces = abs_u_normal_faces[internal_faces_coarse_volume]
            keq_internal = keq_faces[internal_faces_coarse_volume]
            mob_w_int_f = mobility_w_internal_faces[map_coarse_internal_faces]
            mob_o_int_f = mobility_o_internal_faces[map_coarse_internal_faces]
            hi2 = hi[map_coarse_internal_faces]
            areas_internal = areas[internal_faces_coarse_volume]
            ni_internal = len(internal_faces_coarse_volume)
            g_velocity_internal = g_velocity_internal_faces[map_coarse_internal_faces]
            transmissibility_internal = transmissibility_internal_faces[map_coarse_internal_faces]

            local_flux, local_velocity  = self.get_upscale_flux(
                delta_p2,
                abs_u_normal_internal_faces,
                keq_internal,
                mob_w_int_f,
                mob_o_int_f,
                hi2,
                rho_w,
                rho_o,
                gravity_vector,
                areas_internal,
                ni_internal,
                g_velocity_internal,
                transmissibility_internal
            )

            velocity_internal_faces[map_coarse_internal_faces] = local_velocity
            flux_internal_faces[map_coarse_internal_faces] = local_flux
            ###############################################################

            contador += 1

        # pdb.set_trace()

        # return flux_internal_faces, velocity_internal_faces
        # pdb.set_trace()
        return flux_coarse_volumes, test_vector, local_pressure

    def get_upscale_flux(self, delta_p, abs_u_normal_intersect_faces, keq_intersect_faces,
            mob_w_int_f, mob_o_int_f, hi2, rho_w, rho_o, gravity_vector, areas_intersect_faces, ni, g_velocity_intersect_faces, transmissibility_faces):

        delta_p2 = delta_p[:, 1] - delta_p[:, 0]
        delta_p2 = delta_p2.reshape(ni, 1)
        pressure_direction = -delta_p2*abs_u_normal_intersect_faces
        # keq = keq_intersect_faces.reshape(ni, 1)
        # mob_w_int_f = mob_w_int_f.reshape(ni, 1)
        # mob_o_int_f = mob_o_int_f.reshape(ni, 1)
        # mob_t = mob_w_int_f + mob_o_int_f
        # hi2 = hi2.sum(axis=1).reshape(ni, 1)
        # k1 = pressure_direction*(mob_t)*keq
        transm = transmissibility_faces.reshape(ni, 1)
        k1 = pressure_direction*transm
        # k2 = (rho_w*mob_w_int_f + rho_o*mob_o_int_f)*gravity_vector*keq*hi2
        k2 = g_velocity_intersect_faces
        # resp = ((k1 + k2)*abs_u_normal_intersect_faces).sum(axis=1)
        # resp2 = (k1 + k2)*abs_u_normal_intersect_faces
        resp2 = (k1 + k2)
        resp3 = resp2*abs_u_normal_intersect_faces
        resp = resp3.sum(axis=1)*areas_intersect_faces
        return resp, resp2

    def get_local_velocity(self, pressure, volumes, volumes_adj_faces, keq_faces, g_total_velocity_faces):

        pass
