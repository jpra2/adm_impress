from packs.tpfa.biphasic import TpfaBiphasicCons
from packs.tpfa.monophasic import TpfaMonophasic
from packs.data_class import GeometricData, ElementsData, DataManager, SparseDataManager, RockData, BiphasicData, SimulationData, WellsData, AccumulativeBiphasicData, CurrentBiphasicData
from packs.preprocess import TpfaPreprocess
from packs.running.initial_mesh_properties import initial_mesh
from packs.properties import PhisicalProperties
from packs.solvers.solvers_scipy.solver_sp import SolverSp
from copy import deepcopy
import numpy as np
import os
import pdb
import time
from packs.multiscale.test_conservation import ConservationTest
from packs.well_model.create_well import AllWells, Well, get_linear_problem

import scipy.io as sio


test_coarse_flux = False
if test_coarse_flux == True:
    conservation_test = ConservationTest()

M, _, data_impress, wells = initial_mesh()

meshset_volumes = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)
meshset_faces = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset_faces, M.core.all_faces)
dty = [('prod_o', float), ('prod_w', float), ('wor', float),
       ('delta_t', float), ('dvpi', float), ('loop', int), ('global_identifier', int)]

phisical_properties = PhisicalProperties()

def print_test_volumes(file_name):
    data_impress.update_variables_to_mesh()
    M.core.mb.write_file(file_name, [meshset_volumes])

def print_test_faces(file_name):
    data_impress.update_variables_to_mesh()
    M.core.mb.write_file(file_name, [meshset_faces])

###########################################
## preprocessamento
def preprocessar():

    wells2 = WellsData()
    wells2._data = deepcopy(wells._data)

    perms = M.permeability[:]
    perms = perms.reshape(len(perms), 3, 3)

    elements = ElementsData()
    elements.set_mesh_elements(
        M.volumes.all,
        M.faces.all,
        M.edges.all,
        M.nodes.all,
        M.faces.boundary,
        M.nodes.boundary
    )
    elements.create_adj_matrix_volumes_to_faces(elements.volumes, elements.faces, M.volumes.bridge_adjacencies(elements.volumes, 3, 2))
    elements.create_adj_matrix_faces_to_edges(elements.faces, elements.edges, M.faces.bridge_adjacencies(elements.faces, 2, 1))
    elements.create_adj_matrix_edges_to_nodes(elements.edges, elements.nodes, M.edges.bridge_adjacencies(elements.edges, 1, 0))
    elements.insert('volumes_adj_internal_faces', elements.faces_to_volumes(elements.internal_faces), 'array')
    elements.insert('volumes_adj_boundary_faces', elements.faces_to_volumes(elements.boundary_faces).flatten(), 'array')
    elements.insert('volumes_adj_volumes_by_faces', M.volumes.bridge_adjacencies(elements.volumes, 2, 3), 'array')
    map_internal_faces = np.repeat(-1, len(elements.faces))
    map_internal_faces[elements.internal_faces] = np.arange(len(elements.internal_faces))
    elements.insert('map_internal_faces', map_internal_faces, 'array')
    vb = elements.get('volumes_adj_boundary_faces')
    neighborship = np.zeros((len(elements.internal_faces)+len(elements.boundary_faces), 2), dtype=np.int)
    neighborship[elements.internal_faces] = elements.get('volumes_adj_internal_faces')
    neighborship[elements.boundary_faces] = np.array([vb, vb]).T
    elements.insert('neighborship_faces', neighborship, 'array')
    del neighborship, vb

    n_vols = len(elements.volumes)

    preprocess = TpfaPreprocess()
    geom = GeometricData()
    geom['abs_u_normal_faces'] = M.faces.normal[:]
    geom['block_dimension'] = np.repeat(1.0, n_vols*3).reshape(n_vols, 3)
    geom['volume'] = np.prod(geom['block_dimension'], axis=1)
    geom['centroid_volumes'] = M.volumes.center[:]
    geom['centroid_faces'] = M.faces.center[:]
    geom['centroid_edges'] = M.edges.center[:]
    geom['centroid_nodes'] = M.nodes.center[:]
    nodes_of_faces = elements.faces_to_nodes(elements.faces)
    points_faces = geom['centroid_nodes'][nodes_of_faces]
    geom['areas'] = preprocess.get_areas_faces(elements.faces, points_faces)
    geom['u_direction_internal_faces'] = preprocess.get_u_normal_internal_faces(elements.get('volumes_adj_internal_faces'), geom['abs_u_normal_faces'][elements.internal_faces], geom['centroid_volumes'])
    geom['hi'] = preprocess.get_h_internal_faces(geom['block_dimension'], elements.get('volumes_adj_internal_faces'), geom['abs_u_normal_faces'][elements.internal_faces])
    geom['h'] = np.zeros((len(elements.internal_faces) + len(elements.boundary_faces), 2))
    geom['h'][elements.internal_faces] = geom['hi']
    hb = preprocess.get_h_boundary_faces(geom['block_dimension'], elements.get('volumes_adj_boundary_faces'), np.absolute(geom['abs_u_normal_faces'][elements.boundary_faces]))
    geom['h'][elements.boundary_faces] = np.array([hb, hb]).T
    del hb

    rock_data = RockData()
    rock_data['permeability'] = perms.copy()
    rock_data['permeability_volumes_internal_faces_direction'] = preprocess.get_k_volumes_internal_faces_direction(elements.get('volumes_adj_internal_faces'), perms, geom['abs_u_normal_faces'][elements.internal_faces])
    rock_data['porosity'] = M.poro[:].flatten()
    rock_data['keq_faces'] = preprocess.get_equivalent_permeability_faces_from_diagonal_permeability(
        elements.volumes,
        elements.faces,
        elements.get('volumes_adj_internal_faces'),
        elements.get('volumes_adj_boundary_faces'),
        geom['abs_u_normal_faces'],
        rock_data['permeability'],
        elements.internal_faces,
        geom['block_dimension']
    )

    rock_data['k_faces'] = preprocess.get_permeability_faces_from_diagonal_permeability(
        elements.volumes,
        elements.faces,
        elements.get('volumes_adj_internal_faces'),
        elements.get('volumes_adj_boundary_faces'),
        geom['abs_u_normal_faces'],
        rock_data['permeability'],
        elements.internal_faces,
        geom['block_dimension']
    )

    biphasic_data = BiphasicData()
    biphasic_data['saturation'] = M.saturation[:].flatten()
    biphasic_data['one_sided_mobility_faces_w'] = np.zeros(geom['h'].shape)
    biphasic_data['one_sided_mobility_faces_o'] = np.zeros(geom['h'].shape)
    biphasic_data['one_sided_total_mobility_faces'] = np.zeros(geom['h'].shape)

    simulation_data = SimulationData()
    # simulation_data['nkga_internal_faces'] = phisical_properties.get_nkga(rock_data['keq_faces'][elements.internal_faces], geom['u_direction_internal_faces'], geom['areas'][elements.internal_faces])
    simulation_data['nada'] = np.array([-1])


    current_data = CurrentBiphasicData()
    current_data['current'] = np.zeros(1, dtype=dty)

    accumulate = AccumulativeBiphasicData()
    accumulate.create()
    accumulate.insert_data(current_data['current'])

    DataManager.export_all_datas_to_npz()
    SparseDataManager.export_all_datas()

    return wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate
#########################################################

########################################
## carregamento
def carregar():

    wells2 = WellsData(load=True)
    elements = ElementsData(load=True)
    geom = GeometricData(load=True)
    rock_data = RockData(load=True)
    biphasic_data = BiphasicData(load=True)
    simulation_data = SimulationData(load=True)
    current_data = CurrentBiphasicData(load=True)
    accumulate = AccumulativeBiphasicData()
    accumulate.create(global_identifier=current_data['current']['global_identifier'][0])

    return wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate
########################################

wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate = preprocessar()
# wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate = carregar()
# pdb.set_trace()





###############################################
monophasic = TpfaMonophasic()
biphasic = TpfaBiphasicCons()
###############################################

#############################
### well_model
def init_well_model(w_ids, w_centroids, w_block_dimensions, types_prescription,
                    w_values, w_directions, w_types, w_permeabilities, w_mobilities, z_wells):

    list_of_wells = []
    rho_w = biphasic.properties.rho_w
    rho_o = biphasic.properties.rho_o
    for i, w_id in enumerate(w_ids):
        centroids = w_centroids[i]
        block_dimension = w_block_dimensions[i]
        type_prescription = types_prescription[i]
        value = w_values[i]
        direction = w_directions[i]
        w_type = w_types[i]
        w_permeability = w_permeabilities[i]
        mobilities = w_mobilities[i]
        z_well = z_wells[i]

        well = Well()
        well.set_well(w_id, centroids, block_dimension, type_prescription=type_prescription, value=value,
                      direction=direction, well_type=w_type, z_well=z_well, n_phases=2,
                      well_permeability=w_permeability, rho_phases=[rho_w, rho_o])
        # well.update_req(well.calculate_req(well_permeability=w_permeability))
        # well.wi = well.calculate_WI(well.block_dimension, well.req, well.well_radius, w_permeability, well.direction)
        set_phase_infos(well)
        well.mobilities = mobilities

        list_of_wells.append(well)

    return list_of_wells

def load_well_model():
    list_wells = AllWells.load_wells_from_database()

    return list_wells

def set_phase_infos(well):
    # well.update_n_phases(2)
    well.phases = ['water', 'oil']
    # well.rho = [1000, 100]


def setting_well_model():

    saturation = biphasic_data['saturation'].copy()
    saturation[wells['ws_q_sep'][0]] = 1.0
    # saturation[wells['ws_p_sep'][0]] = 1.0
    krw, kro = biphasic.get_krw_and_kro(saturation)
    mob_w, mob_o = biphasic.get_mobilities_w_o(krw, kro)

    list_w_ids = [wells['ws_q_sep'][0], wells['ws_p_sep'][0]]
    # list_w_ids = [wells['ws_p_sep'][0], wells['ws_p_sep'][1]]
    list_w_centroids = [geom['centroid_volumes'][w_id] for w_id in list_w_ids]
    list_w_block_dimensions = [geom['block_dimension'][w_id] for w_id in list_w_ids]
    list_type_prescription = ['flow_rate', 'pressure']
    # list_type_prescription = ['pressure', 'pressure']
    # list_w_values = [10000000.0, 100000.0]
    list_w_values = [10000.0, 100000.0] # para vazao prescrita
    list_w_directions = ['z', 'z']
    list_w_types = ['injector', 'producer']
    list_w_permeability = [rock_data['permeability'][w_id] for w_id in list_w_ids]
    list_w_mobilities = [np.array([mob_w[w_id], mob_o[w_id]]).T for w_id in list_w_ids]
    centroid_nodes = geom['centroid_nodes']
    zmax = centroid_nodes.max(axis=0)[2]
    # list_z_wells = [geom['centroid_volumes'][w_id][:,2].min() for w_id in list_w_ids]
    # list_z_wells = [geom['centroid_volumes'][w_id][:,2].max() for w_id in list_w_ids]
    list_z_wells = [zmax for w_id in list_w_ids]

    all_wells = init_well_model(list_w_ids, list_w_centroids, list_w_block_dimensions,
                           list_type_prescription, list_w_values, list_w_directions,
                           list_w_types, list_w_permeability, list_w_mobilities, list_z_wells)

    AllWells.create_database()

    return all_wells


dT = 1e-9
t = 0
Tmax = 1e-5

well_models = setting_well_model()

# well_models = load_well_model()
######################################################

solver = SolverSp()
###########################################

loop_max = 100000
loop = current_data['current']['loop'][0]

centroids_save = geom['centroid_volumes']
all_pressures = np.zeros(len(centroids_save))
all_saturations = all_pressures.copy()
all_saturations[:] = biphasic_data['saturation'][:]


# pressure_mat = sio.loadmat('data/ex4/pressure.mat')['pr'].flatten()
# centroids_mat = sio.loadmat('data/ex1/centroids.mat')['rr']
# pressure_comp = np.zeros(len(pressure_mat))

while loop <= loop_max:

    biphasic_data['krw'], biphasic_data['kro'] = biphasic.get_krw_and_kro(biphasic_data['saturation'])
    mob_w, mob_o = biphasic.get_mobilities_w_o(biphasic_data['krw'], biphasic_data['kro'])
    for well in well_models:
        well.update_mobilities(np.array([mob_w[well.volumes_ids], mob_o[well.volumes_ids]]).T)

    biphasic_data['one_sided_mobility_faces_w'][elements.internal_faces] = mob_w[elements.get('neighborship_faces')[elements.internal_faces]]
    biphasic_data['one_sided_mobility_faces_w'][elements.boundary_faces] = mob_w[elements.get('neighborship_faces')[elements.boundary_faces]]
    biphasic_data['one_sided_mobility_faces_o'][elements.internal_faces] = mob_o[elements.get('neighborship_faces')[elements.internal_faces]]
    biphasic_data['one_sided_mobility_faces_o'][elements.boundary_faces] = mob_o[elements.get('neighborship_faces')[elements.boundary_faces]]
    biphasic_data['one_sided_total_mobility_faces'] = biphasic_data['one_sided_mobility_faces_w'] + biphasic_data['one_sided_mobility_faces_o']


    if loop == 0:

        biphasic_data['upwind_w'], biphasic_data['upwind_o'] = biphasic.set_initial_upwind_internal_faces(
            biphasic_data['saturation'],
            elements.get('volumes_adj_internal_faces'),
            wells['ws_inj'],
            elements.internal_faces,
            elements.volumes_to_faces(elements.volumes),
            elements.boundary_faces,
            elements.get('map_internal_faces')
        )

    else:

        biphasic_data['upwind_w'], biphasic_data['upwind_o'] = biphasic.update_upwind_phases(
            geom['centroid_volumes'],
            elements.internal_faces,
            elements.get('volumes_adj_internal_faces'),
            biphasic_data['saturation'],
            flux_w_internal_faces,
            flux_o_internal_faces,
            total_flux_internal_faces,
            gravity = phisical_properties._gravity
        )

    upwind_w_vec, upwind_o_vec = biphasic.visualize_upwind_vec(
        elements.internal_faces,
        geom['abs_u_normal_faces'][elements.internal_faces],
        geom['centroid_volumes'],
        elements.get('volumes_adj_internal_faces'),
        biphasic_data['upwind_w'],
        biphasic_data['upwind_o']
    )

    data_impress['upwind_w_faces_vec'][elements.internal_faces] = upwind_w_vec
    data_impress['upwind_o_faces_vec'][elements.internal_faces] = upwind_o_vec

    ###
    biphasic_data['mob_w_internal_faces'] = mob_w[elements.get('volumes_adj_internal_faces')[biphasic_data['upwind_w']]]
    biphasic_data['mob_o_internal_faces'] = mob_o[elements.get('volumes_adj_internal_faces')[biphasic_data['upwind_o']]]

    # biphasic_data['transmissibility_faces'] = biphasic.get_transmissibility_faces(
    #     geom['areas'],
    #     elements.internal_faces,
    #     elements.boundary_faces,
    #     elements.get('volumes_adj_internal_faces'),
    #     elements.get('volumes_adj_boundary_faces'),
    #     biphasic_data['mob_w_internal_faces'],
    #     biphasic_data['mob_o_internal_faces'],
    #     mob_w,
    #     mob_o,
    #     rock_data['keq_faces']
    # )

    biphasic_data['transmissibility_faces'] = biphasic.get_transmissibility_faces(
        geom['areas'],
        elements.internal_faces,
        elements.boundary_faces,
        elements.get('volumes_adj_internal_faces'),
        elements.get('volumes_adj_boundary_faces'),
        mob_w,
        mob_o,
        rock_data['keq_faces'],
        rock_data['k_faces'],
        biphasic_data['one_sided_total_mobility_faces'],
        geom['h']
    )

    biphasic_data['g_velocity_w_internal_faces'], biphasic_data['g_velocity_o_internal_faces'] = biphasic.get_g_velocity_w_o_internal_faces(
        phisical_properties.gravity_vector,
        biphasic_data['mob_w_internal_faces'],
        biphasic_data['mob_o_internal_faces'],
        biphasic.properties.rho_w,
        biphasic.properties.rho_o,
        geom['hi'],
        rock_data['keq_faces'][elements.internal_faces]
    )

    g_total_velocity_internal_faces = biphasic_data['g_velocity_w_internal_faces'] + biphasic_data['g_velocity_o_internal_faces']

    biphasic_data['g_source_w_internal_faces'], biphasic_data['g_source_o_internal_faces'] = biphasic.get_g_source_w_o_internal_faces(
        geom['areas'][elements.internal_faces],
        geom['u_direction_internal_faces'],
        biphasic_data['g_velocity_w_internal_faces'],
        biphasic_data['g_velocity_o_internal_faces']
    )

    wells2.add_gravity_2(
        elements.volumes,
        phisical_properties.gravity_vector,
        geom['centroid_volumes'],
        elements.get('volumes_adj_volumes_by_faces'),
        geom['centroid_nodes'],
        biphasic_data['saturation'],
        biphasic.properties.rho_w,
        biphasic.properties.rho_o
    )

    g_source_total_internal_faces = biphasic_data['g_source_w_internal_faces'] + biphasic_data['g_source_o_internal_faces']
    g_source_total_volumes = phisical_properties.get_total_g_source_volumes(
        elements.volumes,
        elements.get('volumes_adj_internal_faces'),
        g_source_total_internal_faces
    )

    T = monophasic.mount_transmissibility_matrix(
        biphasic_data['transmissibility_faces'][elements.internal_faces],
        elements.internal_faces,
        elements.get('volumes_adj_internal_faces'),
        elements.volumes
    )

    # T_with_boundary, b = monophasic.get_linear_problem(
    #     wells2['ws_p'],
    #     wells2['ws_q'],
    #     wells2['values_p'],
    #     wells2['values_q'],
    #     g_source_total_volumes,
    #     T
    # )

    # import pdb; pdb.set_trace()

    T_with_boundary, b = get_linear_problem(well_models, T, g_source_total_volumes, elements.volumes)

    # import pdb; pdb.set_trace()
    #
    # import scipy.sparse as sp
    # norm_T = sp.linalg.norm(T_with_boundary, np.inf)
    # norm_invT = sp.linalg.norm(sp.linalg.inv(T_with_boundary), np.inf)
    # cond = norm_T*norm_invT

    # import pdb; pdb.set_trace()
    pressure1 = solver.direct_solver(T_with_boundary, b)
    pressure = pressure1[elements.volumes]
    all_pressures = np.vstack((all_pressures, pressure))

    total_velocity_internal_faces = biphasic.get_total_velocity_internal_faces(
        pressure,
        elements.get('volumes_adj_internal_faces'),
        rock_data['keq_faces'][elements.internal_faces],
        biphasic_data['mob_w_internal_faces'],
        biphasic_data['mob_o_internal_faces'],
        geom['abs_u_normal_faces'][elements.internal_faces],
        biphasic_data['g_velocity_w_internal_faces'] + biphasic_data['g_velocity_o_internal_faces']
    )
    # vif = total_velocity_internal_faces

    # import pdb; pdb.set_trace()

    data_impress['velocity_faces'][elements.internal_faces] = total_velocity_internal_faces
    data_impress['pressure'] = pressure

    velocity_w_internal_faces, velocity_o_internal_faces = biphasic.get_velocity_w_and_o_internal_faces(
        total_velocity_internal_faces,
        phisical_properties.gravity_vector,
        rock_data['keq_faces'][elements.internal_faces],
        biphasic.properties.rho_w,
        biphasic.properties.rho_o,
        biphasic_data['mob_w_internal_faces'],
        biphasic_data['mob_o_internal_faces'],
        geom['hi'],
        geom['abs_u_normal_faces'][elements.internal_faces]
    )

    data_impress['velocity_w_faces'][elements.internal_faces] = velocity_w_internal_faces
    data_impress['velocity_o_faces'][elements.internal_faces] = velocity_o_internal_faces

    flux_w_internal_faces = biphasic.test_flux(
        velocity_w_internal_faces,
        geom['abs_u_normal_faces'][elements.internal_faces],
        geom['areas'][elements.internal_faces]
    )

    flux_o_internal_faces = biphasic.test_flux(
        velocity_o_internal_faces,
        geom['abs_u_normal_faces'][elements.internal_faces],
        geom['areas'][elements.internal_faces]
    )

    total_flux_internal_faces = flux_w_internal_faces + flux_o_internal_faces

    if test_coarse_flux is True:
        ml_data = M.multilevel_data
        flux_coarse_volumes_pf, total_flux_internal_faces_pf, velocity_internal_faces_pf, local_pressure_pf = conservation_test.conservation_with_gravity(
            elements.volumes,
            data_impress['GID_1'],
            pressure,
            T,
            ml_data['coarse_faces_level_1'],
            ml_data['coarse_intersect_faces_level_1'],
            geom['areas'],
            ml_data['coarse_primal_id_level_1'],
            elements.get('volumes_adj_internal_faces'),
            ml_data['coarse_internal_faces_level_1'],
            elements.get('map_internal_faces'),
            geom['abs_u_normal_faces'],
            phisical_properties.gravity_vector,
            rock_data['keq_faces'],
            biphasic.properties.rho_w,
            biphasic.properties.rho_o,
            geom['hi'],
            g_source_total_volumes,
            ml_data['vertex_level_1'],
            g_source_total_internal_faces,
            g_total_velocity_internal_faces,
            wells2['ws_p'],
            wells2['values_p'],
            wells2['ws_q'],
            wells2['values_q']
        )
        print()
        print(np.allclose(total_flux_internal_faces_pf, total_flux_internal_faces))
        print(np.allclose(local_pressure_pf, pressure))
        print()
        pdb.set_trace()

    if loop == 0:
        loop += 1
        continue

    p1 = pressure

    flux_total_volumes = monophasic.get_total_flux_volumes(
        total_flux_internal_faces,
        elements.volumes,
        elements.get('volumes_adj_internal_faces')
    )

    # flux_w_volumes = biphasic.get_flux_phase_volumes(
    #     flux_w_internal_faces,
    #     wells['all_wells'],
    #     mob_w/(mob_w+mob_o),
    #     elements.get('volumes_adj_internal_faces'),
    #     elements.volumes,
    #     flux_total_volumes
    # )

    flux_w_volumes = biphasic.get_flux_phase_volumes_with_wells(
        flux_w_internal_faces,
        wells['all_wells'],
        mob_w / (mob_w + mob_o),
        elements.get('volumes_adj_internal_faces'),
        elements.volumes,
        flux_total_volumes
    )

    # flux_o_volumes = biphasic.get_flux_phase_volumes(
    #     flux_o_internal_faces,
    #     wells['all_wells'],
    #     mob_o/(mob_w+mob_o),
    #     elements.get('volumes_adj_internal_faces'),
    #     elements.volumes,
    #     flux_total_volumes
    # )

    flux_o_volumes = biphasic.get_flux_phase_volumes_with_wells(
        flux_o_internal_faces,
        wells['all_wells'],
        mob_o / (mob_w + mob_o),
        elements.get('volumes_adj_internal_faces'),
        elements.volumes,
        flux_total_volumes
    )

    ###########
    flux_phases = AllWells.update_flux_phases(well_models, flux_phase_0=flux_w_volumes, flux_phase_1=flux_o_volumes)
    ## flux_phases[0] = flux_w_volumes
    ## flux_phases[1] = flux_o_volumes
    #############
    flux_w_volumes[:] = flux_phases[0][:] * biphasic.properties.rho_w
    flux_o_volumes[:] = flux_phases[1][:] * biphasic.properties.rho_o

    biphasic_data['saturation_last'] = biphasic_data['saturation'].copy()

    # delta_t, biphasic_data['saturation'] = biphasic.update_saturation(
    #     flux_w_volumes,
    #     rock_data['porosity'],
    #     geom['volume'],
    #     flux_total_volumes,
    #     total_velocity_internal_faces,
    #     mob_w/(mob_w + mob_o),
    #     biphasic_data['saturation_last'],
    #     elements.get('volumes_adj_internal_faces'),
    #     geom['hi'],
    #     geom['abs_u_normal_faces'][elements.internal_faces]
    # )

    delta_t, biphasic_data['saturation'] = biphasic.update_saturation(
        flux_w_volumes,
        rock_data['porosity'],
        geom['volume'],
        flux_total_volumes,
        total_velocity_internal_faces,
        mob_w / (mob_w + mob_o),
        biphasic_data['saturation_last'],
        elements.get('volumes_adj_internal_faces'),
        geom['hi'],
        geom['abs_u_normal_faces'][elements.internal_faces],
        dt=dT,
        loop=loop
    )

    sat = biphasic_data['saturation']
    sata = biphasic_data['saturation_last']
    # print(sat)
    # print()
    #
    import pdb; pdb.set_trace()

    all_saturations = np.vstack((all_saturations, biphasic_data['saturation']))

    prod_w, prod_o, wor, dvpi = biphasic.get_production_w_o(
        flux_total_volumes,
        mob_w,
        mob_o,
        wells2['ws_inj'],
        wells2['ws_prod'],
        delta_t,
        geom['volume'],
        rock_data['porosity']
    )

    loop += 1

    current_data['current'] = np.array([(prod_o, prod_w, wor, delta_t, dvpi, loop, accumulate.global_identifier)], dtype=dty)
    accumulate.insert_data(current_data['current'])

    data_impress['saturation_last'][:] = biphasic_data['saturation_last']
    data_impress['saturation'][:] = biphasic_data['saturation']
    data_impress['flux_volumes'][:] = flux_total_volumes
    data_impress['pressure'][:] = pressure
    data_impress['flux_w_volumes'][:] = flux_w_volumes
    data_impress['flux_o_volumes'][:] = flux_o_volumes

    name = 'results/test_volumes_' + str(loop) + '.vtk'
    # print_test_volumes(name)
    # M.core.print(folder='results', file='test_volumes_'+str(loop), extension='.vtk', config_input='input_cards/print_settings0.yml')

    name2 = 'results/test_faces_' + str(loop) + '.vtk'
    # M.core.print(folder='results', file='test_faces_'+str(loop), extension='.vtk', config_input='input_cards/print_settings2.yml')
    # print_test_faces(name2)

    w1 = well_models[0]
    w2 = well_models[1]
    p1 = pressure[w1.volumes_ids]
    p2 = pressure[w2.volumes_ids]
    AllWells.update_wells_pbh(well_models, pressure1, elements.volumes)
    AllWells.update_wells_flow_rate(well_models, flux_total_volumes)

    '''
        update: wells flow_rate
        genererate wells report
    '''

    pwells = pressure1[elements.volumes.max() + 1:]
    q1 = flux_total_volumes[w1.volumes_ids]
    q2 = flux_total_volumes[w2.volumes_ids]

    # if loop >= 2:
    #     import pdb; pdb.set_trace()
    # import pdb; pdb.set_trace()

    # if loop % 100 == 0:
    #     accumulate.export(local_key_identifier='loop')
    #     current_data.export_to_npz()
    #     current_data['current']['global_identifier'] = accumulate.global_identifier
    #
    # if loop % 200 == 0:
    #     all_datas = accumulate.load_all_datas()
    #     pdb.set_trace()
    #
    # if loop % 40 == 0:
    #     import pdb; pdb.set_trace()

    t += dT

    if t >= Tmax:
        loop = loop_max + 1

    # if loop == 94:
    #     loop = loop_max + 1



# np.save(os.path.join('data', 'pressures_p7_g.npy'), all_pressures)
# np.save(os.path.join('data', 'pressures_p7_g_3.npy'), all_pressures)
# np.save(os.path.join('data', 'centroids.npy'), centroids_save)
# np.save(os.path.join('data', 'saturations_p7_g_3.npy'), all_saturations)

import pdb; pdb.set_trace()





















    # pdb.set_trace()












    # name = 'results/test_volumes_' + str(loop) + '.vtk'
    # print_test_volumes(name)




    # pdb.set_trace()

pdb.set_trace()
