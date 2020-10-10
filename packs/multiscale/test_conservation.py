import numpy as np
import pdb
from packs.tpfa.biphasic import TpfaBiphasicCons
from packs.common_files.common_infos import CommonInfos
from scipy import sparse as sp
from packs.solvers.solvers_scipy.solver_sp import SolverSp
from packs.properties import PhisicalProperties
from numba import jit


class ConservationTest:

    def conservation_with_gravity(self,
        volumes, primal_ids, pressure, global_transmissibility,
        coarse_faces, coarse_intersect_faces, areas, coarse_primal_ids_org, volumes_adj_internal_faces,
        coarse_internal_faces, map_internal_faces, abs_u_normal_faces,
        gravity_vector, keq_faces,
        rho_w, rho_o, hi, g_source_total_volumes, coarse_vertexes, g_source_total_internal_faces, g_velocity_internal_faces, faces_volumes, boundary_faces,
        wells_p=[], vals_p=[], wells_q=[], vals_q=[], loop=0):

        min_value = 1e-8

        all_primal_ids = np.unique(primal_ids)
        cids = coarse_primal_ids_org.copy()
        flux_coarse_volumes = np.zeros(len(all_primal_ids))
        # test_vector = flux_coarse_volumes.copy()
        n_volumes = len(volumes)
        local_ids = np.zeros(len(volumes), dtype=int)
        local_pressure = local_ids.astype(float)
        n_internal_faces = len(volumes_adj_internal_faces)
        velocity_internal_faces = np.zeros((n_internal_faces, 3))
        flux_internal_faces = np.zeros(n_internal_faces)
        solver = SolverSp()
        transmissibility_internal_faces = np.absolute(global_transmissibility[volumes_adj_internal_faces[:,0], volumes_adj_internal_faces[:,1]].toarray().flatten())
        set_p = set(wells_p)
        set_q = set(wells_q)
        set_all_wells = set_p | set_q
        if len(set_p) > 0:
            presc_p = np.zeros(n_volumes)
            presc_p[wells_p] = vals_p

        if len(wells_q) > 0:
            presc_q = np.zeros(n_volumes)
            presc_q[wells_q] = vals_q

        flux2 = np.zeros(n_internal_faces)
        flux3 = flux2.copy()
        pare = True

        for primal_id in all_primal_ids:

            flux2[:] = 0.0
            coarse_flux_2 = 0.0
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
            abs_u_normal_intersect_faces = abs_u_normal_faces[intersect_coarse_faces]
            areas_intersect_faces = areas[intersect_coarse_faces]
            g_velocity_intersect_faces = g_velocity_internal_faces[map_intersect_faces]
            ni = len(intersect_coarse_faces)
            adj_intersect_intern_vols = np.zeros(ni, dtype=int)
            transmissibility_intersect = transmissibility_internal_faces[map_intersect_faces]

            for i, vadjs in enumerate(volumes_adj_intersect):
                if set([vadjs[0]]) & set_vols:
                    adj_intersect_intern_vols[i] = vadjs[0]
                else:
                    adj_intersect_intern_vols[i] = vadjs[1]

            # vols2 = np.setdiff1d(vols, adj_intersect_intern_vols[i])

            resp, resp2 = self.get_upscale_flux(
                delta_p,
                abs_u_normal_intersect_faces,
                areas_intersect_faces,
                ni,
                g_velocity_intersect_faces,
                transmissibility_intersect
            )

            velocity_internal_faces[map_intersect_faces] = resp2
            flux_internal_faces[map_intersect_faces] = resp
            flux2[map_intersect_faces] = resp

            if set_all_wells & set(vols):
                flux3[:] = 0.0
                volsp3 = np.array(list(set_all_wells & set(vols)))
                faces_volsp3 = faces_volumes[volsp3]
                all_intfaces_volsp3 = np.setdiff1d(np.unique(faces_volsp3.flatten()), boundary_faces)
                volumes_adjs_faces_volsp3 = volumes_adj_internal_faces[map_internal_faces[all_intfaces_volsp3]]

                delta_p3 = pressure[volumes_adjs_faces_volsp3]
                abs_u_normal_faces_volsp3 = abs_u_normal_faces[all_intfaces_volsp3]
                areas_faces_volsp3 = areas[all_intfaces_volsp3]
                g_velocity_faces_volsp3 = g_velocity_internal_faces[map_internal_faces[all_intfaces_volsp3]]
                ni3 = len(all_intfaces_volsp3)
                transmissibility_faces_volsp3 = transmissibility_internal_faces[map_internal_faces[all_intfaces_volsp3]]

                flux_faces_volsp3, vel_faces_volsp3 = self.get_upscale_flux(
                    delta_p3,
                    abs_u_normal_faces_volsp3,
                    areas_faces_volsp3,
                    ni3,
                    g_velocity_faces_volsp3,
                    transmissibility_faces_volsp3
                )

                flux3[map_internal_faces[all_intfaces_volsp3]] = flux_faces_volsp3

                presc_flux_wells = PhisicalProperties.get_total_g_source_volumes(volumes, volumes_adj_internal_faces, flux3)[volsp3]
                coarse_flux_2 = presc_flux_wells.sum()

                volsp3[:] = -1


            #
            # ################################################################
            # neuman problemg_velocity_intersect_faces
            presc_flux_vols = PhisicalProperties.get_total_g_source_volumes(volumes, volumes_adj_internal_faces, flux2)[vols]
            c_flux = presc_flux_vols.sum() - coarse_flux_2
            # if abs(c_flux) > min_value:
            # if loop>6 and pare==True:
            #     pdb.set_trace()

            # pdb.set_trace()

            flux_coarse_volumes[primal_id] = c_flux


            adj_internal_faces_coarse = volumes_adj_internal_faces[map_coarse_internal_faces]
            local_adj_internal = local_ids[adj_internal_faces_coarse]
            g_source_internal = g_source_total_internal_faces[map_coarse_internal_faces]
            g_source_local_volumes = PhisicalProperties.get_total_g_source_volumes(local_ids[vols], local_adj_internal, g_source_internal)
            # import pdb; pdb.set_trace()
            local_source_term = presc_flux_vols + g_source_local_volumes

            if set_p & set(vols):
                volsp = np.intersect1d(vols, wells_p)
                T_local[local_ids[volsp]] = 0
                T_local[local_ids[volsp], local_ids[volsp]] = 1
                local_source_term[local_ids[volsp]] = presc_p[volsp]

            if set_q & set(vols):
                volsq = np.intersect1d(vols, wells_q)
                volsq = np.setdiff1d(volsq, [vertex])
                local_source_term[local_ids[volsq]] += presc_q[volsq]

            local_source_term[local_ids[vertex]] = pressure[vertex]
            T_local[local_ids[vertex]] = 0
            T_local[local_ids[vertex], local_ids[vertex]] = 1
            p2 = solver.direct_solver(T_local, local_source_term)
            local_pressure[vols] = p2
            ################################################################
            delta_p2 = p2[local_adj_internal]
            abs_u_normal_internal_faces = abs_u_normal_faces[internal_faces_coarse_volume]
            areas_internal = areas[internal_faces_coarse_volume]
            ni_internal = len(internal_faces_coarse_volume)
            g_velocity_internal = g_velocity_internal_faces[map_coarse_internal_faces]
            transmissibility_internal = transmissibility_internal_faces[map_coarse_internal_faces]

            local_flux, local_velocity  = self.get_upscale_flux(
                delta_p2,
                abs_u_normal_internal_faces,
                areas_internal,
                ni_internal,
                g_velocity_internal,
                transmissibility_internal
            )

            velocity_internal_faces[map_coarse_internal_faces] = local_velocity
            flux_internal_faces[map_coarse_internal_faces] = local_flux
            ###############################################################

        return flux_coarse_volumes, flux_internal_faces, velocity_internal_faces, local_pressure
        # return flux_coarse_volumes, test_vector, local_pressure

    def get_upscale_flux(self, delta_p, abs_u_normal_intersect_faces,
            areas_intersect_faces, ni, g_velocity_intersect_faces, transmissibility_faces):

        delta_p2 = delta_p[:, 1] - delta_p[:, 0]
        delta_p2 = delta_p2.reshape(ni, 1)
        pressure_direction = -delta_p2*abs_u_normal_intersect_faces
        transm = transmissibility_faces.reshape(ni, 1)
        v1 = pressure_direction*transm
        v2 = g_velocity_intersect_faces
        velocity = (v1 + v2)
        flux = velocity*abs_u_normal_intersect_faces
        flux = flux.sum(axis=1)*areas_intersect_faces
        return flux, velocity
