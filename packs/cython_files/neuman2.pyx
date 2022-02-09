import numpy as np
from packs.solvers.solvers_scipy.solver_sp import SolverSp
from packs.properties import PhisicalProperties
from packs.common_files.common_infos import CommonInfos


def get_upscale_flux(
    delta_p,
    abs_u_normal_intersect_faces,
    areas_intersect_faces,
    ni,
    g_velocity_intersect_faces,
    transmissibility_faces):

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


def neuman_with_gravity(
    volumes,
    primal_ids,
    pressure,
    global_transmissibility,
    coarse_faces,
    coarse_intersect_faces,
    areas,
    coarse_primal_ids_org,
    volumes_adj_internal_faces,
    coarse_internal_faces,
    map_internal_faces,
    abs_u_normal_faces,
    g_source_total_volumes,
    coarse_vertexes,
    g_source_total_internal_faces,
    g_velocity_internal_faces,
    velocity_internal_faces,
    wells_p=[],
    vals_p=[],
    wells_q=[],
    vals_q=[]):


    all_primal_ids = np.unique(primal_ids)
    cids = coarse_primal_ids_org.copy()
    flux_coarse_volumes = np.zeros(len(all_primal_ids), dtype=np.double)

    test_vector = flux_coarse_volumes.copy()

    n_volumes = len(volumes)
    local_ids = np.zeros(len(volumes), dtype=np.long)

    local_pressure = np.zeros(n_volumes)

    n_internal_faces = len(volumes_adj_internal_faces)
    flux_internal_faces = np.zeros(n_internal_faces)
    # velocity_internal_faces = np.zeros((n_internal_faces, 3))
    # cdef double [:, :, :] velocity_internal_faces_view = velocity_internal_faces
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

    flux2 = np.zeros(n_internal_faces)

    for primal_id in all_primal_ids:

        flux2[:from packs.cython_files.neuman2 import neuman_with_gravity

t0 = time.time()
fluxc2, test_vector2, local_pressure2 = neuman_with_gravity(
    elements.volumes.astype(int),
    data_impress['GID_1'].astype(int),
    padm,
    T,
    ml_data['coarse_faces_level_1'],
    ml_data['coarse_intersect_faces_level_1'],
    geom['areas'],
    ml_data['coarse_primal_id_level_1'].astype(np.int),
    elements.get('volumes_adj_internal_faces').astype(np.int),
    ml_data['coarse_internal_faces_level_1'],
    elements.get('map_internal_faces').astype(np.int),
    geom['abs_u_normal_faces'],
    g_source_total_volumes,
    ml_data['vertex_level_1'].astype(np.int),
    g_source_total_internal_faces,
    g_total_velocity_internal_faces,
    np.zeros((len(elements.internal_faces), 3)),
    wells2['ws_p'].astype(np.int),
    wells2['values_p'],
    wells2['ws_q'].astype(np.int),
    wells2['values_q']
)
t1 = time.time()
print(t1 - t0)] = 0.0
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

        resp, resp2 = get_upscale_flux(
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

        #
        # ################################################################
        # neuman problemg_velocity_intersect_faces
        presc_flux_vols = PhisicalProperties.get_total_g_source_volumes(volumes, volumes_adj_internal_faces, flux2)[vols]
        flux_coarse_volumes[primal_id] = presc_flux_vols.sum()

        adj_internal_faces_coarse = volumes_adj_internal_faces[map_coarse_internal_faces]
        local_adj_internal = local_ids[adj_internal_faces_coarse]
        g_source_internal = g_source_total_internal_faces[map_coarse_internal_faces]
        g_source_local_volumes = PhisicalProperties.get_total_g_source_volumes(local_ids[vols], local_adj_internal, g_source_internal)
        local_source_term = g_source_local_volumes + presc_flux_vols

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
        # local_pressure[vols] = p2
        local_pressure[vols] = p2
        ################################################################
        delta_p2 = p2[local_adj_internal]
        abs_u_normal_internal_faces = abs_u_normal_faces[internal_faces_coarse_volume]
        areas_internal = areas[internal_faces_coarse_volume]
        ni_internal = len(internal_faces_coarse_volume)
        g_velocity_internal = g_velocity_internal_faces[map_coarse_internal_faces]
        transmissibility_internal = transmissibility_internal_faces[map_coarse_internal_faces]

        local_flux, local_velocity  = get_upscale_flux(
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

    # return flux_internal_faces, velocity_internal_faces
    return flux_coarse_volumes, test_vector, local_pressure

    # return 0, 0, 0
