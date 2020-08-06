from packs.tpfa.biphasic import TpfaBiphasicCons
from packs.tpfa.monophasic import TpfaMonophasic
from packs.data_class import GeometricData, ElementsData, DataManager, SparseDataManager, RockData, BiphasicData, SimulationData, WellsData, AccumulativeBiphasicData, CurrentBiphasicData
from packs.preprocess import TpfaPreprocess
from packs.running.initial_mesh_properties import initial_mesh
from packs.properties import PhisicalProperties
from packs.solvers.solvers_scipy.solver_sp import SolverSp
from packs.multiscale.operators.prolongation.AMS.ams_tpfa_new1 import AMSTpfa
from packs.multiscale.operators.prolongation.AMS.paralel_ams_new1 import MasterOP
from copy import deepcopy
import numpy as np
import os
import pdb
import time

M, _, data_impress, wells = initial_mesh()
ml_data = M.multilevel_data
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

print_test_volumes('test1.vtk')

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

    biphasic_data = BiphasicData()
    biphasic_data['saturation'] = M.saturation[:].flatten()

    simulation_data = SimulationData()
    simulation_data['nkga_internal_faces'] = phisical_properties.get_nkga(rock_data['keq_faces'][elements.internal_faces], geom['u_direction_internal_faces'], geom['areas'][elements.internal_faces])


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

# wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate = preprocessar()
wells2, elements, geom, rock_data, biphasic_data, simulation_data, current_data, accumulate = carregar()
# pdb.set_trace()




###############################################
monophasic = TpfaMonophasic()
biphasic = TpfaBiphasicCons()
###############################################

solver = SolverSp()

ams_tpfa_level1 = AMSTpfa(elements.volumes, data_impress['GID_1'], data_impress['DUAL_1'])
###########################################

loop_max = 100000
loop = current_data['current']['loop'][0]

while loop <= loop_max:

    biphasic_data['krw'], biphasic_data['kro'] = biphasic.get_krw_and_kro(biphasic_data['saturation'])
    mob_w, mob_o = biphasic.get_mobilities_w_o(biphasic_data['krw'], biphasic_data['kro'])

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
            total_flux_internal_faces
        )

    ###
    biphasic_data['mob_w_internal_faces'] = mob_w[elements.get('volumes_adj_internal_faces')[biphasic_data['upwind_w']]]
    biphasic_data['mob_o_internal_faces'] = mob_o[elements.get('volumes_adj_internal_faces')[biphasic_data['upwind_o']]]
    biphasic_data['transmissibility_faces'] = biphasic.get_transmissibility_faces(
        geom['areas'],
        elements.internal_faces,
        elements.boundary_faces,
        elements.get('volumes_adj_internal_faces'),
        elements.get('volumes_adj_boundary_faces'),
        biphasic_data['upwind_w'],
        biphasic_data['upwind_o'],
        mob_w,
        mob_o,
        rock_data['keq_faces']
    )

    biphasic_data['g_source_w_internal_faces'], biphasic_data['g_source_o_internal_faces'] = phisical_properties.get_g_source_w_o_internal_faces(
        simulation_data['nkga_internal_faces'],
        biphasic_data['mob_w_internal_faces'],
        biphasic_data['mob_o_internal_faces'],
        biphasic.properties.rho_w,
        biphasic.properties.rho_o
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

    T_with_boundary, b = monophasic.get_linear_problem(
        wells2['ws_p'],
        wells2['ws_q'],
        wells2['values_p'],
        wells2['values_q'],
        g_source_total_volumes,
        T
    )

    # OP1, pcorr1 = ams_tpfa_level1.run(T, total_source_term=g_source_total_volumes)
    master = MasterOP(T, ml_data['dual_structure_level_1'], 1)
    OP2, pcorr2 = master.run()
    del master
    OP1, pcorr1 = ams_tpfa_level1.run(T)

    pressure = solver.direct_solver(T_with_boundary, b)
    total_velocity_internal_faces = biphasic.get_total_velocity_internal_faces(
        pressure,
        elements.internal_faces,
        phisical_properties.gravity_vector,
        elements.get('volumes_adj_internal_faces'),
        rock_data['keq_faces'][elements.internal_faces],
        biphasic_data['mob_w_internal_faces'],
        biphasic_data['mob_o_internal_faces'],
        geom['centroid_volumes'],
        biphasic.properties.rho_w,
        biphasic.properties.rho_o
    )

    # total_flux_internal_faces = biphasic.get_total_flux_internal_faces(
    #     total_velocity_internal_faces,
    #     geom['u_direction_internal_faces'],
    #     geom['areas'][elements.internal_faces]
    # )

    # flux_w_internal_faces, flux_o_internal_faces = biphasic.get_flux_w_and_o_internal_faces(
    #     total_flux_internal_faces,
    #     biphasic.properties.rho_w,
    #     biphasic.properties.rho_o,
    #     biphasic_data['mob_w_internal_faces'],
    #     biphasic_data['mob_o_internal_faces'],
    #     simulation_data['nkga_internal_faces']
    # )

    velocity_w_internal_faces, velocity_o_internal_faces = biphasic.get_velocity_w_and_o_internal_faces(
        total_velocity_internal_faces,
        phisical_properties.gravity_vector,
        rock_data['keq_faces'][elements.internal_faces],
        biphasic.properties.rho_w,
        biphasic.properties.rho_o,
        biphasic_data['mob_w_internal_faces'],
        biphasic_data['mob_o_internal_faces']
    )

    flux_w_internal_faces = biphasic.test_flux(
        velocity_w_internal_faces,
        geom['u_direction_internal_faces'],
        geom['areas'][elements.internal_faces]
    )

    flux_o_internal_faces = biphasic.test_flux(
        velocity_o_internal_faces,
        geom['u_direction_internal_faces'],
        geom['areas'][elements.internal_faces]
    )

    total_flux_internal_faces = flux_w_internal_faces + flux_o_internal_faces

    flux_total_volumes = monophasic.get_total_flux_volumes(
        total_flux_internal_faces,
        elements.volumes,
        elements.get('volumes_adj_internal_faces')
    )

    flux_w_volumes = biphasic.get_flux_phase_volumes(
        flux_w_internal_faces,
        wells['all_wells'],
        mob_w/(mob_w+mob_o),
        elements.get('volumes_adj_internal_faces'),
        elements.volumes,
        flux_total_volumes
    )

    flux_o_volumes = biphasic.get_flux_phase_volumes(
        flux_o_internal_faces,
        wells['all_wells'],
        mob_o/(mob_w+mob_o),
        elements.get('volumes_adj_internal_faces'),
        elements.volumes,
        flux_total_volumes
    )

    biphasic_data['saturation_last'] = biphasic_data['saturation'].copy()

    delta_t, biphasic_data['saturation'] = biphasic.update_saturation(
        flux_w_volumes,
        rock_data['porosity'],
        geom['volume'],
        flux_total_volumes,
        total_velocity_internal_faces,
        mob_w/(mob_w + mob_o),
        biphasic_data['saturation_last'],
        elements.get('volumes_adj_internal_faces'),
        geom['hi'],
        geom['abs_u_normal_faces'][elements.internal_faces]
    )

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

    # data_to_export = {
    #     'prod_o': np.array([prod_o]),
    #     'prod_w': np.array([prod_w]),
    #     'wor': np.array([wor]),
    #     'delta_t': np.array([delta_t]),
    #     'dvpi': np.array([dvpi]),
    #     'loop': np.array([loop])
    # }
    loop += 1

    current_data['current'] = np.array([(prod_o, prod_w, wor, delta_t, dvpi, loop, accumulate.global_identifier)], dtype=dty)
    accumulate.insert_data(current_data['current'])

    data_impress['saturation'] = biphasic_data['saturation_last']
    data_impress['flux_volumes'] = flux_total_volumes
    data_impress['pressure'] = pressure
    data_impress['flux_w_volumes'] = flux_w_volumes
    data_impress['flux_o_volumes'] = flux_o_volumes

    name = 'results/test_volumes_' + str(loop) + '.vtk'
    print_test_volumes(name)

    pdb.set_trace()

    if loop % 5 == 0:
        accumulate.export(local_key_identifier='loop')
        current_data.export_to_npz()
        current_data['current']['global_identifier'] = accumulate.global_identifier

    if loop % 20 == 0:
        all_datas = accumulate.load_all_datas()
        pdb.set_trace()

























    # pdb.set_trace()












    # name = 'results/test_volumes_' + str(loop) + '.vtk'
    # print_test_volumes(name)




    # pdb.set_trace()

pdb.set_trace()
