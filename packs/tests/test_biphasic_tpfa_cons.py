from packs.tpfa.biphasic import TpfaBiphasicCons
from packs.tpfa.monophasic import TpfaMonophasic
from packs.data_class import GeometricData, Elements, RockData, SparseDataManager, DataManager, BiphasicData, SimulationData
from packs.preprocess import TpfaPreprocess
from packs.running.initial_mesh_properties import initial_mesh
from packs.properties import PhisicalProperties
from packs.solvers.solvers_scipy.solver_sp import SolverSp
import numpy as np
import pdb





###########################################
## preprocessamento
M, _, data_impress, wells = initial_mesh()
meshset_volumes = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset_volumes, M.core.all_volumes)
meshset_faces = M.core.mb.create_meshset()
M.core.mb.add_entities(meshset_faces, M.core.all_faces)

def print_test_volumes(file_name):
    data_impress.update_variables_to_mesh()
    M.core.mb.write_file(name, [meshset_volumes])



perms = M.permeability[:]
perms = perms.reshape(len(perms), 3, 3)

elements = Elements()
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
elements.insert_value('volumes_adj_internal_faces', elements.faces_to_volumes(elements.internal_faces), 'array')
elements.insert_value('volumes_adj_boundary_faces', elements.faces_to_volumes(elements.boundary_faces).flatten(), 'array')
elements.insert_value('volumes_adj_volumes_by_faces', M.volumes.bridge_adjacencies(elements.volumes, 2, 3), 'array')
map_internal_faces = np.repeat(-1, len(elements.faces))
map_internal_faces[elements.internal_faces] = np.arange(len(elements.internal_faces))
elements.insert_value('map_internal_faces', map_internal_faces, 'array')

n_vols = len(elements.volumes)

preprocess = TpfaPreprocess()
geom = GeometricData()
geom['abs_u_normal_faces'] = M.faces.normal
geom['block_dimension'] = np.repeat(1.0, n_vols*3).reshape(n_vols, 3)
geom['volume'] = np.prod(geom['block_dimension'], axis=1)
geom['centroid_volumes'] = M.volumes.center[:]
geom['centroid_faces'] = M.faces.center[:]
geom['centroid_edges'] = M.edges.center[:]
geom['centroid_nodes'] = M.nodes.center[:]
nodes_of_faces = elements.faces_to_nodes(elements.faces)
points_faces = geom['centroid_nodes'][nodes_of_faces]
geom['areas'] = preprocess.get_areas_faces(elements.faces, points_faces)
geom['u_direction_internal_faces'] = preprocess.get_u_normal_internal_faces(elements.get_value('volumes_adj_internal_faces'), geom['abs_u_normal_faces'][elements.internal_faces], geom['centroid_volumes'])

rock_data = RockData()
rock_data['permeability'] = perms.copy()
rock_data['permeability_volumes_internal_faces_direction'] = preprocess.get_k_volumes_internal_faces_direction(elements.get_value('volumes_adj_internal_faces'), perms, geom['abs_u_normal_faces'][elements.internal_faces])
rock_data['porosity'] = M.poro[:].flatten()
rock_data['keq_faces'] = preprocess.get_equivalent_permeability_faces_from_diagonal_permeability(
    elements.volumes,
    elements.faces,
    elements.get_value('volumes_adj_internal_faces'),
    elements.get_value('volumes_adj_boundary_faces'),
    geom['abs_u_normal_faces'],
    rock_data['permeability'],
    elements.internal_faces,
    geom['block_dimension']
)

phisical_properties = PhisicalProperties()

biphasic_data = BiphasicData()
biphasic_data['saturation'] = M.saturation[:].flatten()

monophasic = TpfaMonophasic()
biphasic = TpfaBiphasicCons()
biphasic_data['krw'], biphasic_data['kro'] = biphasic.get_krw_and_kro(biphasic_data['saturation'])
# krw, kro = biphasic.get_krw_and_kro(biphasic_data['saturation'])
mob_w, mob_o = biphasic.get_mobilities_w_o(biphasic_data['krw'], biphasic_data['kro'])
###
##initialize
biphasic_data['upwind_w'], biphasic_data['upwind_o'] = biphasic.set_initial_upwind_internal_faces(
    biphasic_data['saturation'],
    elements.get_value('volumes_adj_internal_faces'),
    wells['ws_inj'],
    elements.internal_faces,
    elements.volumes_to_faces(elements.volumes),
    elements.boundary_faces,
    elements.get_value('map_internal_faces')
)
###
biphasic_data['mob_w_internal_faces'] = mob_w[elements.get_value('volumes_adj_internal_faces')[biphasic_data['upwind_w']]]
biphasic_data['mob_o_internal_faces'] = mob_o[elements.get_value('volumes_adj_internal_faces')[biphasic_data['upwind_o']]]
biphasic_data['transmissibility_faces'] = biphasic.get_transmissibility_faces(
    geom['areas'],
    elements.internal_faces,
    elements.boundary_faces,
    elements.get_value('volumes_adj_internal_faces'),
    elements.get_value('volumes_adj_boundary_faces'),
    biphasic_data['upwind_w'],
    biphasic_data['upwind_o'],
    mob_w,
    mob_o,
    rock_data['keq_faces']
)

simulation_data = SimulationData()
simulation_data['nkga_internal_faces'] = phisical_properties.get_nkga(rock_data['keq_faces'][elements.internal_faces], geom['u_direction_internal_faces'], geom['areas'][elements.internal_faces])

biphasic_data['g_source_w_internal_faces'], biphasic_data['g_source_o_internal_faces'] = phisical_properties.get_g_source_w_o_internal_faces(
    simulation_data['nkga_internal_faces'],
    biphasic_data['mob_w_internal_faces'],
    biphasic_data['mob_o_internal_faces'],
    biphasic.properties.rho_w,
    biphasic.properties.rho_o
)

wells.add_gravity_2(
    elements.volumes,
    phisical_properties.gravity_vector,
    geom['centroid_volumes'],
    elements.get_value('volumes_adj_volumes_by_faces'),
    geom['centroid_nodes'],
    biphasic_data['saturation'],
    biphasic.properties.rho_w,
    biphasic.properties.rho_o
)



solver = SolverSp()
###########################################

loop_max = 30
loop = 0

while loop <= loop_max:

    g_source_total_internal_faces = biphasic_data['g_source_w_internal_faces'] + biphasic_data['g_source_o_internal_faces']
    g_source_total_volumes = phisical_properties.get_total_g_source_volumes(
        elements.volumes,
        elements.get_value('volumes_adj_internal_faces'),
        g_source_total_internal_faces
    )
    
    data_impress['flux_grav_volumes'] = g_source_total_volumes
    data_impress['saturation'] = saturation

    T = monophasic.mount_transmissibility_matrix(
        biphasic_data['transmissibility_faces'][elements.internal_faces],
        elements.internal_faces,
        elements.get_value('volumes_adj_internal_faces'),
        elements.volumes
    )

    T_with_boundary, b = monophasic.get_linear_problem(
        wells['ws_p'],
        wells['ws_q'],
        wells['values_p'],
        wells['values_q'],
        g_source_total_volumes,
        T
    )

    data_impress['pressure'] = solver.direct_solver(T_with_boundary, b)
    total_velocity_internal_faces = biphasic.get_total_velocity_internal_faces(
        data_impress['pressure'],
        elements.internal_faces,
        phisical_properties.gravity_vector,
        elements.get_value('volumes_adj_internal_faces'),
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
        elements.get_value('volumes_adj_internal_faces')
    )

    flux_w_volumes = biphasic.get_flux_phase_volumes(
        flux_w_internal_faces,
        wells['all_wells'],
        mob_w,
        elements.get_value('volumes_adj_internal_faces'),
        elements.volumes,
        flux_total_volumes
    )

    flux_o_volumes = biphasic.get_flux_phase_volumes(
        flux_o_internal_faces,
        wells['all_wells'],
        mob_o,
        elements.get_value('volumes_adj_internal_faces'),
        elements.volumes,
        flux_total_volumes
    )
















    pdb.set_trace()












    name = 'results/test_volumes_' + str(loop) + '.vtk'
    print_test_volumes(name)




    pdb.set_trace()











    loop += 1

pdb.set_trace()
