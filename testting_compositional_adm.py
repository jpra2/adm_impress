from packs.directories import data_loaded
from run_compositional_adm import RunSimulationAdm
import time
from packs.multiscale.preprocess.dual_primal.create_dual_and_primal_mesh import MultilevelData
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.compositional.compositional_params import Params
from packs.adm.non_uniform.adm_method_non_nested import AdmNonNested
from packs.multiscale.preprocess.prep_neumann import NeumannSubdomains
from packs.multiscale.preprocess.dual_domains import create_dual_subdomains
from packs.utils import constants as ctes
import numpy as np
import scipy.sparse as sp

""" ---------------- LOAD STOP CRITERIA AND MESH DATA ---------------------- """

name_current = 'current_compositional_results_'
name_all = data_loaded['name_save_file'] + '_'
mesh = 'mesh/' + data_loaded['mesh_name']
n_levels = data_loaded['n_levels']
get_correction_term = data_loaded['get_correction_term']
load_operators = data_loaded['load_operators']
load_multilevel_data = data_loaded['load_multilevel_data']
params = Params()

if data_loaded['use_vpi']:
    stop_criteria = max(data_loaded['compositional_data']['vpis_para_gravar_vtk'])
else: stop_criteria = data_loaded['compositional_data']['maximum_time']

loop_max = 1000
run_criteria = 0# ###############
# from packs.multiscale.preprocess.dual_primal import generate_primal
# from packs.multiscale.preprocess.dual_primal import create_coarse_volumes_structure

# volumes_dimension = generate_primal.get_dimension_volumes(M.volumes.bridge_adjacencies(M.volumes.all, 3, 0), M.nodes.center[:])
# multilevel_info = {
#     # 'volumes': elements_lv0['volumes'],
#     'centroids_volumes': M.volumes.center[:],
#     'volumes_dimension': volumes_dimension,
#     # 'nodes': M.nodes.all,
#     'centroids_nodes': M.nodes.center[:],
#     'cr': np.array(data_loaded['Crs']['Cr1']),
#     # 'volumes_nodes': M.volumes.bridge_adjacencies(M.volumes.all, 3, 0),
#     'adjacencies_internal_faces': elements_lv0['neig_internal_faces'],
#     # 'internal_faces': elements_lv0['internal_faces'],
#     # 'boundary_faces': elements_lv0['boundary_faces'],
#     # 'volumes_faces': M.volumes.bridge_adjacencies(M.volumes.all, 3, 2)
# }
# primal_ids, dual_ids = generate_primal.create_dual_and_primal(**multilevel_info)

# coarse_info = {
#     # 'fine_volumes': elements_lv0['volumes'],
#     'centroids_fine_volumes': M.volumes.center[:],
#     'primal_id': primal_ids, 
#     'adjacencies_internal_faces': elements_lv0['neig_internal_faces'],
#     # 'fine_nodes': M.nodes.all, 
#     'fine_volumes_nodes': M.volumes.bridge_adjacencies(M.volumes.all, 3, 0), 
#     'centroids_fine_nodes': M.nodes.center[:],
#     'fine_internal_faces': elements_lv0['internal_faces'],
#     'fine_boundary_faces': elements_lv0['boundary_faces'],
#     # 'fine_volumes_faces': M.volumes.bridge_adjacencies(M.volumes.all, 3, 2),
#     # 'fine_centroids_faces': M.faces.center[:],
#     'adjacencies_boundary_faces': M.faces.bridge_adjacencies(M.faces.boundary[:], 2, 3)
# }
# create_coarse_volumes_structure.create_coarse_structure(**coarse_info)
# import pdb; pdb.set_trace()
# ###############- RUN CODE --------------------------------- """

load = data_loaded['load_data']
convert = data_loaded['convert_english_to_SI']

t = time.time()
sim = RunSimulationAdm(name_current, name_all)
M, data_impress, wells, fprop, load, elements_lv0 = sim.initialize(load, convert, mesh)
# import pdb; pdb.set_trace()
# load_multilevel_data = False


# ###############
# from packs.multiscale.preprocess.dual_primal import generate_primal
# from packs.multiscale.preprocess.dual_primal import create_coarse_volumes_structure

# volumes_dimension = generate_primal.get_dimension_volumes(M.volumes.bridge_adjacencies(M.volumes.all, 3, 0), M.nodes.center[:])
# multilevel_info = {
#     # 'volumes': elements_lv0['volumes'],
#     'centroids_volumes': M.volumes.center[:],
#     'volumes_dimension': volumes_dimension,
#     # 'nodes': M.nodes.all,
#     'centroids_nodes': M.nodes.center[:],
#     'cr': np.array(data_loaded['Crs']['Cr1']),
#     # 'volumes_nodes': M.volumes.bridge_adjacencies(M.volumes.all, 3, 0),
#     'adjacencies_internal_faces': elements_lv0['neig_internal_faces'],
#     # 'internal_faces': elements_lv0['internal_faces'],
#     # 'boundary_faces': elements_lv0['boundary_faces'],
#     # 'volumes_faces': M.volumes.bridge_adjacencies(M.volumes.all, 3, 2)
# }
# primal_ids, dual_ids = generate_primal.create_dual_and_primal(**multilevel_info)

# coarse_info = {
#     # 'fine_volumes': elements_lv0['volumes'],
#     'centroids_fine_volumes': M.volumes.center[:],
#     'primal_id': primal_ids, 
#     'adjacencies_internal_faces': elements_lv0['neig_internal_faces'],
#     # 'fine_nodes': M.nodes.all, 
#     'fine_volumes_nodes': M.volumes.bridge_adjacencies(M.volumes.all, 3, 0), 
#     'centroids_fine_nodes': M.nodes.center[:],
#     'fine_internal_faces': elements_lv0['internal_faces'],
#     'fine_boundary_faces': elements_lv0['boundary_faces'],
#     # 'fine_volumes_faces': M.volumes.bridge_adjacencies(M.volumes.all, 3, 2),
#     # 'fine_centroids_faces': M.faces.center[:],
#     'adjacencies_boundary_faces': M.faces.bridge_adjacencies(M.faces.boundary[:], 2, 3)
# }
# create_coarse_volumes_structure.create_coarse_structure(**coarse_info)
# import pdb; pdb.set_trace()
# ###############

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

mlo = MultilevelOperators(n_levels, data_impress, elements_lv0, ml_data, load=load_operators, get_correction_term=get_correction_term)
neumann_subds = NeumannSubdomains(elements_lv0, ml_data, data_impress, wells)
adm = AdmNonNested(wells['all_wells'], n_levels, M, data_impress, elements_lv0)
# ml_data.load_tags()
# import pdb; pdb.set_trace()
dual_subdomains = create_dual_subdomains(ml_data['dual_structure_level_1'], ml_data['fine_dual_id_level_1'], ml_data['fine_primal_id_level_1'])
global_vector_update = np.full(ctes.n_volumes, True, dtype=bool)
ncoarse_ids = len(np.unique(data_impress['GID_1']))
OP_AMS = sp.lil_matrix((ctes.n_volumes, ncoarse_ids)).tocsc()

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
    'dual_subdomains': dual_subdomains
}

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
    
    global_vector_update[:] = True # update the prolongation operator


    sim.run(M, wells, fprop, load, 
            multilevel_data=ml_data, 
            multilevel_operators=mlo,
            params=params,
            adm_method=adm,
            neumann_subds=neumann_subds,
            data_impress=data_impress,
            elements_lv0=elements_lv0,
            **local_problem_params)

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

tf = time.time()
print('Total computational time: ', tf-t) #total simulation time
import pdb; pdb.set_trace()
sim.save_infos(data_impress, M) #Save data to file
