from packs.multiscale.correction_function.CF import get_correction_function
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.multiscale.operators.prolongation.AMS.Paralell.group_dual_volumes import group_dual_volumes_and_get_OP
import scipy.sparse as sp
import os

import time



def multiscale_cf_calc(g_source_total_volumes, volumes_without_grav, local_lu_and_global_ids, data_impress, As, cf_keyword=None, update_cf=False):
    
    t0 = time.time()
    b2 = g_source_total_volumes.copy()
    b2[volumes_without_grav] = 0
    cfs = get_correction_function(local_lu_and_global_ids, As, b2)
    
    if cf_keyword:
        data_impress[cf_keyword][:] = cfs
    
    print('cf {}'.format(time.time()-t0))
    
    return cfs
    

def multiscale_prolongation_operator(multilevel_operators: MultilevelOperators, load_operators: bool, multilevel_data, T, T_with_boundary, M, data_impress, neta_lim, ams_tpfa, local_lu_matrices, separated_dual_structures, transmissibility, As=None, update_op=False):
    if not update_op:
        pass
    else:
        multilevel_operators.run_paralel(T, M.multilevel_data['dual_structure_level_1'], 0, False)
        # OP=group_dual_volumes_and_get_OP(multilevel_operators, T_with_boundary, M, data_impress, T, neta_lim=neta_lim)
        As = ams_tpfa.get_as_off_diagonal(ams_tpfa.get_Twire(T))
        local_lu_matrices.update_lu_objects(separated_dual_structures, transmissibility)
        multilevel_operators['prolongation_lcd_level_1']=sp.find(multilevel_operators['prolongation_level_1'])
    
    return As


def print_mesh_volumes_results(M, data_impress, meshset_volumes, result_name):
    data_impress.update_variables_to_mesh()
    # M.core.mb.write_file('results/testt_00'+'.vtk', [meshset_volumes])
    M.core.mb.write_file(os.path.join('results', result_name), [meshset_volumes])
    import pdb; pdb.set_trace()


def cf_remove_volumes_by_level(cf, volumes):
    cf[volumes] = 0
    return cf