from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh

def insert_prolongation_operator_impress(
        m_object: FineScaleMesh,
        prolongation,
        level,
        gid_level,
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
        var[:] = prolongation[:,i].toarray()
        # TODO
        # falta implementar para n levels
        ## esta so da malha fina para o level 1
        ##########

def insert_restriction_operator_impress(
        m_object: FineScaleMesh,
        restriction,
        level,
        gid_level,
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
        var[:] = restriction[i,:].toarray().flatten().reshape([n_fine_vols, 1]).astype(int)
        # TODO
        # falta implementar para n levels
        ## esta so da malha fina para o level 1
        ##########




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

    return solution


def print_mesh_volumes_data(m_object: FineScaleMesh, file_name):
    m1 = m_object.core.mb.create_meshset()
    m_object.core.mb.add_entities(m1, m_object.core.all_volumes)
    m_object.core.mb.write_file(file_name, [m1])