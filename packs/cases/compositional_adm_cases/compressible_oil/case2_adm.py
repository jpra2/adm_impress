import copy

from scipy.sparse import data
from packs.data_class import sparse_operators
from packs.directories import data_loaded
from run_compositional_adm import RunSimulationAdm
import time
from packs.multiscale.preprocess.dual_primal.create_dual_and_primal_mesh import MultilevelData
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.compositional.compositional_params import Params
from packs.adm.non_uniform.adm_method_non_nested import AdmNonNested
from packs.multiscale.preprocess.prep_neumann import NeumannSubdomains
from packs.multiscale.preprocess.dual_domains import DualSubdomain, create_dual_subdomains
from packs.utils import constants as ctes
import numpy as np
import scipy.sparse as sp
from packs.data_class.compositional_data import CompositionalData
from packs.data_class.compositional_cumulative_datamanager import CumulativeCompositionalDataManager
from packs.cases.compositional_adm_cases.compressible_oil import functions_update
from packs.multiscale.neuman_local_problems.master_local_solver import MasterLocalSolver
from packs.multiscale.operators.prolongation.AMS.paralell2.paralel_ams_new_2 import MasterLocalOperator
from packs.data_class.sparse_operators import SparseOperators

""" ---------------- LOAD STOP CRITERIA AND MESH DATA ---------------------- """


name_current = 'current_compositional_results_'
name_all = data_loaded['name_save_file'] + '_'
mesh = 'mesh/' + data_loaded['mesh_name']
n_levels = data_loaded['n_levels']
get_correction_term = data_loaded['get_correction_term']
load_operators = data_loaded['load_operators']
load_multilevel_data = data_loaded['load_multilevel_data']

# description = 'case1_finescale_'
# description = 'case2_adm_'
# description = 'case4_adm_3k'
description = 'case5_adm_3k'
compositional_data = CompositionalData(description=description)
manage_operators = SparseOperators(description=description)
cumulative_compositional_datamanager = CumulativeCompositionalDataManager(description=description)
cumulative_compositional_datamanager.create()
# datas_comp = cumulative_compositional_datamanager.load_all_datas()
# cumulative_compositional_datamanager.delete_all_datas()
# compositional_data.delete()
# import pdb; pdb.set_trace()
loop_array = np.zeros(
    1,
    dtype=[
        ('loop', int),
        ('t', float),
        ('vpi', float),
        ('simulation_time', float),
        ('oil_production', float),
        ('gas_production', float),
        ('n_volumes_update_base_functions', int),
        ('total_volumes_updated', int)
    ]
)
params = Params()

if data_loaded['use_vpi']:
    stop_criteria = max(data_loaded['compositional_data']['vpis_para_gravar_vtk'])
else: stop_criteria = data_loaded['compositional_data']['maximum_time']

loop_max = 1000
run_criteria = 0# ###############

# ###############- RUN CODE --------------------------------- """

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']

t = time.time()
sim = RunSimulationAdm(name_current, name_all)
M, data_impress, wells, fprop, load, elements_lv0 = sim.initialize(load, convert, mesh)
# import pdb; pdb.set_trace()
# load_multilevel_data = False


ml_data = MultilevelData(data_impress, M, load=load_multilevel_data, n_levels=n_levels)

# import pdb; pdb.set_trace()
# data_impress['DUAL_1'][:] = dual_ids
# data_impress['GID_1'][:] = primal_ids
# data_impress.update_variables_to_mesh()
# m = M.core.mb.create_meshset()
# M.core.mb.add_entities(m, M.core.all_volumes)
# M.core.mb.write_file('results/test_primals.vtk', [m])
# import pdb; pdb.set_trace()




ml_data.run()
data_impress.update_variables_to_mesh()

mlo = MultilevelOperators(n_levels, data_impress, elements_lv0, ml_data, load=load_operators, get_correction_term=get_correction_term)
neumann_subds = NeumannSubdomains(elements_lv0, ml_data, data_impress, wells)
adm = AdmNonNested(wells['all_wells'], n_levels, M, data_impress, elements_lv0)
# ml_data.load_tags()
# import pdb; pdb.set_trace()
dual_subdomains = create_dual_subdomains(ml_data['dual_structure_level_1'], ml_data['fine_dual_id_level_1'], ml_data['fine_primal_id_level_1'])
global_vector_update = np.full(ctes.n_volumes, True, dtype=bool)
ncoarse_ids = len(np.unique(data_impress['GID_1']))
OP_AMS = sp.lil_matrix((ctes.n_volumes, ncoarse_ids)).tocsc()

master_neumann = MasterLocalSolver(neumann_subds.neumann_subds, ctes.n_volumes)
master_local_operator = MasterLocalOperator(dual_subdomains, ctes.n_volumes)


params['area'] = data_impress['area']
params['pretransmissibility'] = data_impress['pretransmissibility']

local_problem_params = {
    'Vbulk': ctes.Vbulk,
    'porosity': ctes.porosity,
    'Cf': ctes.Cf,
    'dVtdP': None,
    'P': None,
    'n_volumes': ctes.n_volumes,
    'n_components': ctes.n_components,
    'n_phases': 3,
    'internal_faces_adjacencies': ctes.v0,
    'dVtdk': None,
    'z_centroids': ctes.z,
    'xkj_internal_faces': None,
    'Csi_j_internal_faces': None,
    'mobilities_internal_faces': None,
    'pretransmissibility_internal_faces': ctes.pretransmissibility_internal_faces,
    'Pcap': None,
    'Vp': None,
    'Vt': None,
    'well_volumes_flux_prescription': wells['ws_q'],
    'values_flux_prescription': wells['values_q'],
    'delta_t': None,
    'g': ctes.g,
    'well_volumes_pressure_prescription': wells['ws_p'],
    'pressure_prescription': wells['values_p'],
    'bhp_ind': ctes.bhp_ind,
    'rho_j': None,
    'rho_j_internal_faces': None,
    'm_object': M,
    'global_vector_update': global_vector_update,
    'OP_AMS': OP_AMS,
    'dual_subdomains': dual_subdomains,
    'master_neumann': master_neumann,
    'master_local_operator': master_local_operator
}

latest_mobility = np.zeros(fprop.mobilities.shape)
global_vector_update[:] = False
total_volumes_updated = copy.deepcopy(global_vector_update)
data_impress['LEVEL'][:] = 1

n_loops_for_acumulate = 1
n_loops_for_export = 500

assert (n_loops_for_export % n_loops_for_acumulate) == 0



while run_criteria < stop_criteria:# and loop < loop_max:
    params['pressure'] = fprop.P
    params['mobilities'] = fprop.mobilities
    params['composition'] = fprop.Csi_j
    params['mol_number'] = fprop.Nk
    params['Sg'] = fprop.Sg
    params['Sw'] = fprop.Sw
    params['So'] = fprop.So
    params['z'] = fprop.z
    params['porous_volume'] = fprop.Vp
    params['total_volume'] = fprop.Vt
    
    # global_vector_update[:] = True # update the prolongation operator in all dual volumes
    for phase in range(ctes.n_phases):
        functions_update.update_global_vector_for_latest_variable(
            global_vector_update,
            latest_mobility[:, phase, :],
            fprop.mobilities[:, phase, :],
            0.1
        )
        
    for comp in range(fprop.z.shape[0]):
        functions_update.update_global_vector_for_volumes_adjacencies_variable(
            global_vector_update, 
            elements_lv0['neig_internal_faces'], 
            fprop.z[comp, :], 
            0.1
        )
        
    for dual in dual_subdomains:
        dual: DualSubdomain
        if dual.test_update():
            latest_mobility[:, :, dual.gids] = fprop.mobilities[:, :, dual.gids]
            total_volumes_updated[dual.gids] = True
    
    functions_update.set_level0_delta_sat(
        data_impress['LEVEL'],
        fprop.Sg,
        elements_lv0['neig_internal_faces'],
        0.1
    )
    
    functions_update.set_level0_delta_sat(
        data_impress['LEVEL'],
        fprop.So,
        elements_lv0['neig_internal_faces'],
        0.1
    )
    
    functions_update.set_level0_delta_sat(
        data_impress['LEVEL'],
        fprop.Sw,
        elements_lv0['neig_internal_faces'],
        0.1
    )
        

    t0 = time.time()
    sim.run(M, wells, fprop, load, 
            multilevel_data=ml_data, 
            multilevel_operators=mlo,
            params=params,
            adm_method=adm,
            neumann_subds=neumann_subds,
            data_impress=data_impress,
            elements_lv0=elements_lv0,
            **local_problem_params)
    simulation_time = time.time() - t0

    if data_loaded['use_vpi']:
        'If using time-step unit as vpi'
        run_criteria = sim.vpi
    else:
        'If using time-step unit as second'
        run_criteria = sim.t
        if sim.time_save[-1] == 0.0:
            'if time_save = [0.0] only, it means that all time-steps are going \
            to be saved'
            t_next = sim.t + sim.delta_t
        else:
            'This is made so time to save is equal to the simulation time - this \
            only works (for now) if time to save is in seconds and not vpi'
            t_next = sim.time_save[sim.time_save > sim.t]
            if len(t_next)>1: t_next = t_next[0]

        'If the current simulation time plus the computed time-step is bigger \
        than the final simulation time, correct the time-step so the current \
        simulation time plus delta_t is equal to the final time'
        if sim.t + sim.delta_t > t_next:
            sim.delta_t = t_next - sim.t

    params['mobilities_internal_faces'] = fprop.mobilities_internal_faces
    params['composition_internal_faces'] = fprop.Csi_j_internal_faces
    params['xkj_internal_faces'] = fprop.xkj_internal_faces
    params['rho_phase_internal_faces'] = fprop.rho_j_internal_faces
    
    
    loop = sim.loop
    print(sim.t)
    
    if (loop) % n_loops_for_acumulate == 0:
             
        loop_array['loop'][0] = loop
        loop_array['t'][0] = sim.t
        loop_array['vpi'][0] = sim.vpi
        loop_array['simulation_time'][0] = simulation_time
        loop_array['oil_production'][0] = sim.oil_production
        loop_array['gas_production'][0] = sim.gas_production
        loop_array['n_volumes_update_base_functions'][0] = global_vector_update.sum()
        loop_array['total_volumes_updated'][0] = total_volumes_updated.sum()
        compositional_data.update({
            'pressure': fprop.P,
            'Sg': fprop.Sg,
            'Sw': fprop.Sw,
            'So': fprop.So,
            'global_composition': fprop.z,
            'mols': fprop.Nk,
            'xkj': fprop.xkj,
            'Vp': fprop.Vp,
            'latest_mobility': latest_mobility,
            'loop_array': loop_array
        })
        cumulative_compositional_datamanager.insert_data(compositional_data._data)
        
        manage_operators.update(
            {
                'prolongation_level_1': OP_AMS
            }
        )
    
    if (loop) % n_loops_for_export == 0:
        compositional_data.export_to_npz()
        cumulative_compositional_datamanager.export()
        manage_operators.export()
         
    global_vector_update[:] = False
    total_volumes_updated[:] = False
    data_impress['LEVEL'][:] = 1
        
    
    if loop % 2500 == 0:
        print('sleeping...')
        time.sleep(5)
        

tf = time.time()
print('Total computational time: ', tf-t) #total simulation time
import pdb; pdb.set_trace()
sim.save_infos(data_impress, M) #Save data to file
