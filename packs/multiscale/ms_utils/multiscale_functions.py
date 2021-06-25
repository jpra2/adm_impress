from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh
import numpy as np
import scipy.sparse as sp
from packs.compositional.IMPEC.global_pressure_solver import GlobalIMPECPressureSolver as Gips


def insert_prolongation_operator_impress(
        m_object: FineScaleMesh,
        prolongation,
        level,
        fine_scale_gid,
        previous_gid_level,
        coarse_ids=None
):
    range_coarse_volumes = range(prolongation.shape[1])
    n_fine_vols = prolongation.shape[0]
    if coarse_ids:
        range_coarse_volumes = coarse_ids

    dict_tag = {
        'var_type': 'volumes',
        'data_size': 1,
        'data_format': 'float',
        'data_density': 'sparse',
        'entity_index': None,
        'create': True
    }

    for i in range_coarse_volumes:
        tag_name = 'PROLONGATION_LEVEL_' + str(level) + '_COARSE_ID_' + str(i)
        dict_tag.update({'name_tag': tag_name})
        var = m_object.create_variable(**dict_tag)
        data = prolongation[:, i].toarray()
        insert_data_from_level(
            var,
            data,
            n_fine_vols,
            fine_scale_gid,
            previous_gid_level,
            level
        )


def insert_restriction_operator_impress(
        m_object: FineScaleMesh,
        restriction,
        level,
        fine_scale_gid,
        previous_gid_level,
        coarse_ids=None
):
    range_coarse_volumes = range(restriction.shape[0])
    n_fine_vols = restriction.shape[1]
    if coarse_ids:
        range_coarse_volumes = coarse_ids

    dict_tag = {
        'var_type': 'volumes',
        'data_size': 1,
        'data_format': 'int',
        'data_density': 'sparse',
        'entity_index': None,
        'create': True
    }

    for i in range_coarse_volumes:
        tag_name = 'RESTRICTION_LEVEL_' + str(level) + '_COARSE_ID_' + str(i)
        dict_tag.update({'name_tag': tag_name})
        var = m_object.create_variable(**dict_tag)
        data = restriction[i, :].toarray().flatten().reshape([n_fine_vols, 1]).astype(int)
        insert_data_from_level(
            var,
            data,
            n_fine_vols,
            fine_scale_gid,
            previous_gid_level,
            level
        )


def insert_data_from_level(moab_variable, data, n_fine_vols, fine_scale_gid, previous_gid_level, level):
    if level > 1:
        for gid_fine in range(n_fine_vols):
            if data[gid_fine] == 0:
                continue
            moab_variable[fine_scale_gid[previous_gid_level == gid_fine]] = data[gid_fine]
    else:
        moab_variable[:] = data


def multilevel_pressure_solver(
        fine_scale_transmissibility,
        fine_scale_source_term,
        prolongation_list: list,
        restriction_list: list
):
    from packs.solvers.solvers_scipy.solver_sp import SolverSp
    solver = SolverSp()

    assert isinstance(prolongation_list, list)
    assert isinstance(restriction_list, list)
    n = len(prolongation_list)
    assert n == len(restriction_list)

    t2 = fine_scale_transmissibility.copy()
    q2 = fine_scale_source_term.copy()

    for OP, OR in zip(prolongation_list, restriction_list):
        assert OP.shape[0] == OR.shape[1]
        assert OP.shape[1] == OR.shape[0]
        assert t2.shape[0] == OP.shape[0]
        t2 = OR * (t2 * OP)
        q2 = OR * q2

    solution = solver.direct_solver(t2, q2)
    prolongation_list.reverse()

    for OP in prolongation_list:
        solution = OP * solution

    prolongation_list.reverse()

    return solution


def print_mesh_volumes_data(m_object: FineScaleMesh, file_name):
    m1 = m_object.core.mb.create_meshset()
    m_object.core.mb.add_entities(m1, m_object.core.all_volumes)
    m_object.core.mb.write_file(file_name, [m1])


def update_local_problem(neumann_subds_list, fine_scale_transmissibility_no_bc, diagonal_term, adm_pressure,
                         Ft_internal_faces, map_internal_faces, gids_volumes, adjacencies_internal_faces, all_coarse_intersect_faces, **kwargs):
    """
        @param neumann_subds_list: local primal subdomains structure
        @param fine_scale_transmissibility_no_bc: fine scale transmissibility without boundary conditions
        @param diagonal_term: diagonal term of fine scale transmissibility
        @param adm_pressure: adm pressure in fine scale
        @param Ft_internal_faces: total flux internal faces
        @param map_internal_faces: map for internal faces
        @param gids_volumes: all volumes gids
        @param adjacencies_internal_faces: volumes adjacencies of internal faces
        @param all_coarse_intersect_faces: all faces in intersection between coarse volumes
    """

    v0 = adjacencies_internal_faces
    elements_lv0 = kwargs.get('elements_lv0')
    m_object = kwargs.get('m_object')

    Tglobal = fine_scale_transmissibility_no_bc
    ft_internal_faces_for_prescription = np.zeros(Ft_internal_faces.shape)
    ft_internal_faces_for_prescription[:, map_internal_faces[all_coarse_intersect_faces]] = Ft_internal_faces[:, map_internal_faces[all_coarse_intersect_faces]]
    global_molar_flux_prescription = Gips.update_flux(
        ft_internal_faces_for_prescription,
        kwargs['rho_j_internal_faces'],
        kwargs['mobilities_internal_faces'],
        kwargs['Pcap'],
        v0,
        kwargs['z_centroids'],
        kwargs['pretransmissibility_internal_faces'],
        kwargs['g'],
        kwargs['xkj_internal_faces'],
        kwargs['Csi_j_internal_faces'],
        kwargs['n_components'],
        kwargs['n_volumes']
    )
    m_object.data['flux_volumes'][:] = global_molar_flux_prescription[0]

    for subd in neumann_subds_list:
        subd.Tlocal_no_bc = update_local_transmissibility(Tglobal, subd.volumes, diagonal_term[subd.volumes])
        subd.adm_pressure = adm_pressure[subd.volumes]
        local_flux_prescription = global_molar_flux_prescription[:, subd.volumes]

        subd.local_rhs = Gips.mount_independent_term(
            kwargs['Vbulk'][subd.volumes],
            kwargs['porosity'][subd.volumes],
            kwargs['Cf'],
            kwargs['dVtdP'][subd.volumes],
            kwargs['P'][subd.volumes],
            len(subd.volumes),
            kwargs['n_components'],
            kwargs['n_phases'],
            subd.map_volumes[subd.adj_intern_local_faces],
            kwargs['dVtdk'][:, subd.volumes],
            kwargs['z_centroids'][subd.volumes],
            kwargs['xkj_internal_faces'][:, :, elements_lv0['remaped_internal_faces'][subd.intern_local_faces]],
            kwargs['Csi_j_internal_faces'][:, :, elements_lv0['remaped_internal_faces'][subd.intern_local_faces]],
            kwargs['mobilities_internal_faces'][: , :, elements_lv0['remaped_internal_faces'][subd.intern_local_faces]],
            kwargs['pretransmissibility_internal_faces'][elements_lv0['remaped_internal_faces'][subd.intern_local_faces]],
            kwargs['Pcap'][:, subd.volumes],
            kwargs['Vp'][subd.volumes],
            kwargs['Vt'][subd.volumes],
            subd.map_volumes[subd.ind_neum],
            local_flux_prescription[:, subd.map_volumes[subd.ind_neum]],
            kwargs['delta_t'],
            kwargs['g'],
            subd.map_volumes[subd.ind_diric],
            subd.adm_pressure[subd.map_volumes[subd.ind_diric]],
            subd.map_volumes[subd.ind_diric],
            kwargs['rho_j'][:, :, subd.volumes],
            kwargs['rho_j_internal_faces'][:, :, elements_lv0['remaped_internal_faces'][subd.intern_local_faces]]
        )
        # local_rhs += subd.flux_prescription
        # local_rhs[subd.map_volumes[subd.ind_diric]] = adm_pressure[subd.ind_diric]
        # subd.local_rhs = local_rhs
        subd.Tlocal = set_matrix_pressure_prescription(
            subd.Tlocal_no_bc,
            subd.map_volumes[subd.ind_diric]
        )


def update_local_transmissibility(Tglobal, gids, diagonal_term) -> sp.csc_matrix:
    """
        @param Tglobal: global fine scale transmissibility
        @param gids: global ids of local volumes
        @param diagonal_term: diagonal term of local volumes
        @return: Tlocal local transmisssibility matrix
    """
    n = len(gids)
    lines = np.repeat(gids, n).reshape(n, n)
    cols = np.tile(gids, n).reshape(n, n)

    Tlocal = Tglobal[lines, cols]
    Tlocal.setdiag(np.zeros(n))
    diagonal = np.array(-Tlocal.sum(axis=1)).flatten() + diagonal_term
    Tlocal.setdiag(diagonal)
    return Tlocal.tocsc()


def update_flux_prescription(n_gids, ft_internal_faces, intersect_faces, adjacencies_intersect_faces,
                             map_internal_faces, volumes_in_primal):
    """

    @param n_gids: number of global gids
    @param ft_internal_faces: total flux internal faces
    @param intersect_faces: intersect faces of primal
    @param adjacencies_intersect_faces: volumes adjacencies of intersect faces
    @param map_internal_faces: remap internal faces
    @param volumes_in_primal: fine volumes in primal coarse volume
    @return: local flux prescription
    """
    # k = 4.42707412e3
    k = 1
    v0 = adjacencies_intersect_faces
    flux_internal_faces = ft_internal_faces[0][map_internal_faces[intersect_faces]]

    flux2 = np.bincount(
        np.concatenate([v0[:, 0], v0[:, 1]]),
        weights=np.concatenate([-flux_internal_faces, flux_internal_faces])
    )
    return flux2[volumes_in_primal]*k

def set_matrix_pressure_prescription(transmissibility_matrix_no_bc: sp.csc_matrix, gids_pressure_prescription):
    t2 = transmissibility_matrix_no_bc.copy()
    t2[gids_pressure_prescription] = 0
    t2[gids_pressure_prescription, gids_pressure_prescription] = 1
    return t2
