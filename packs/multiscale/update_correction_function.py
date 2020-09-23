from packs.multiscale.correction_function.CF import get_correction_function


def get_cfs(g_source_total_volumes, volumes_without_grav_level_0, wells2, data_impress, local_lu_matrices, As):

    t0=time.time()
    b2 = g_source_total_volumes.copy()
    b2[volumes_without_grav_level_0] = 0
    b2[wells2['ws_p']] = 0.0


    # b2 = b.copy()
    # b2[wells2['ws_p']] = 0.0
    # b2[wells[]]
    # b2[data_impress['LEVEL']==0] = 0.0

    # b2[data_impress['DUAL_1'] == 2] = 0.0

    # b2=np.ones_like(b2)10000
    # b2[wells['ws_q']] = -b[wells['ws_q']]
    # b2[wells['ws_p']] = 0
    # cfs = get_correction_function(local_lu_matrices.local_lu_and_global_ids, As, np.ones_like(b2))
    cfs = get_correction_function(local_lu_matrices.local_lu_and_global_ids, As, b2)
    # cfs[data_impress['LEVEL']==0] = 0.0
    # cfs[data_impress['LEVEL'] == 0] = 0.0
    # cfs[:] = 0.0
    # cfs = get_correction_function(local_lu_matrices.local_lu_and_global_ids, As, np.zeros_like(b2))
    # cfs[wells['ws_p']] = 0
    # cfs[:] = 0
    data_impress['gama'][:] = cfs
    print('cf {}'.format(time.time()-t0))
    return cfs
