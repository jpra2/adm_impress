from packs.biphasic.relative_perm.brooks_and_corey import BrooksAndCorey
from packs.biphasic.mobility import BiphasicMobility
from packs.biphasic.unstructured.mobility_mesh_elements import direct_edges_mobility
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths
from packs.manager import MeshProperty, MeshData, BoundaryConditions, SimulationData
from packs.multiscale.unstructured.test.test_cross import set_weights_nodes, set_fine_transmissibility
# from packs.multiscale.unstructured.test.test_brazil import define_faces_in_losangle, set_permeability
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.mpfa_methods.weight_interpolation.lpew import get_lpew2_weights
from packs.mpfa_methods.flux_calculation.diamond_method import get_xi_params_ds_flux

from packs.examples.biphasic_mpfa import initial_loop, while_loop, update_data

import os
import numpy as np
from typing import Tuple
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from packs.utils.permfields import chueh_perm_artur_paper, random_permeability_chueh, random_permeability_chueh_v2
from packs.utils.utils_old import is_point_inside_circle

def get_properties():

    fine_mesh_properties_name = 'finescale4_ameba'
    fine_mesh_path = defpaths.ameba_fine4

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)

    return fine_properties, fine_mesh_path

def set_boundary_conditions(fine_properties: MeshProperty):
    bc = BoundaryConditions()

    initial_saturation = np.zeros(fine_properties['faces'].shape[0])
    faces_of_nodes = fine_properties['faces_of_nodes']
    nodes_centroids = fine_properties['nodes_centroids']

    nodes_presc = fine_properties['physical_vertex_201']
    centroids_nodes_presc = nodes_centroids[nodes_presc]

    x_node1 = 0.3
    # x_node2 = 0.9
    x_node3 = 1.7

    dists = np.absolute(centroids_nodes_presc[:, 0] - x_node1) 

    node1 = nodes_presc[dists <= dists.min()]

    # dists[:] = np.absolute(centroids_nodes_presc[:, 0] - x_node2)
    # node2 = nodes_presc[dists <= dists.min()]

    dists[:] = np.absolute(centroids_nodes_presc[:, 0] - x_node3)
    node3 = nodes_presc[dists <= dists.min()]

    node_p0 = node3
    nodes_q01 = node1

    faces_node_p0 = faces_of_nodes[node_p0[0]]

    faces_nodes_q01 = faces_of_nodes[node1[0]]    

    faces_pressure = faces_node_p0
    pressure_presc = 0*np.ones(faces_pressure.shape[0])

    areas_faces_q01 = fine_properties['areas'][faces_nodes_q01]
    neummann_presc_faces_q01 = 1.0*areas_faces_q01/areas_faces_q01.sum()

    faces_neumann = faces_nodes_q01
    neummann_presc_faces = neummann_presc_faces_q01

    bc.set_boundary('dirichlet_volumes', faces_pressure, pressure_presc)

    bc.set_boundary('neumann_volumes', faces_neumann, neummann_presc_faces)

    walls_edges = fine_properties.boundary_edges

    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    fine_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    bc.set_boundary('water_saturation_volumes', faces_neumann, np.ones(faces_neumann.shape[0]))
    initial_saturation[faces_neumann] = 1.0

    bc.set_boundary('injectors', faces_neumann, np.array([True]))
    bc.set_boundary('producers', faces_pressure, np.array([True]))
    bc.set_boundary('initial_saturation', fine_properties['faces'], initial_saturation)

    bc.update_zero_bcs()

    return bc

def get_R(theta):
    R = np.array([
        np.array([np.cos(theta),  np.sin(theta)]),
        np.array([-np.sin(theta), np.cos(theta)])
    ])
    return R

def set_permeability(fine_mesh_path, fine_properties: MeshProperty, export_permfield=True, update_permfield=True, **kwargs):


    tag_preprocess = 'permeability'
    if fine_properties.verify_name_in_data_names(tag_preprocess) and update_permfield is False:
        return

    faces = fine_properties['faces']
    faces_centroids = fine_properties['faces_centroids']
    # N = 150
    # state = 15
    # N = 130
    # state = 92

    N = 30
    # state = 95
    state = 97


    aditional_points = np.array([
        [0.6, 0.3],
        [1.1, 0.5],
        [0.75, 0.3]
    ])

    aditional_ids = []
    for point in aditional_points:
        dist2 = np.linalg.norm(faces_centroids - point, axis=1)
        aditional_ids.append(faces[dist2<=dist2.min()][0])
    
    aditional_ids = np.array(aditional_ids)


    permeability = np.zeros((faces.shape[0], 2, 2))
    # perm = random_permeability_chueh(faces_centroids, N, state, aditional_ids)
    perm = random_permeability_chueh_v2(faces_centroids, N, state, aditional_ids)
    # perm = chueh_perm_artur_paper(faces_centroids)

    permeability[:, 0, 0] = perm
    permeability[:, 1, 1] = perm
    
    fine_properties.insert_or_update_data({
        tag_preprocess: permeability
    })

    if export_permfield:
        mesh_data = MeshData(mesh_path=fine_mesh_path)
        mesh_data.create_tag('permeability_xx')
        mesh_data.create_tag('permeability_xy')
        mesh_data.create_tag('permeability_yx')
        mesh_data.create_tag('permeability_yy')
        mesh_data.insert_tag_data(
            'permeability_xx',
            fine_properties['permeability'][:, 0, 0],
            elements_type='faces'
        )
        mesh_data.insert_tag_data(
            'permeability_xy',
            fine_properties['permeability'][:, 0, 1],
            elements_type='faces'
        )
        mesh_data.insert_tag_data(
            'permeability_yx',
            fine_properties['permeability'][:, 1, 0],
            elements_type='faces'
        )
        mesh_data.insert_tag_data(
            'permeability_yy',
            fine_properties['permeability'][:, 1, 1],
            elements_type='faces'
        )
        mesh_data.export_all_elements_type_to_vtk('permfield', element_type='faces')

def initial_funcs(
        fp: MeshProperty,
        fine_mesh_path: str
):
    
    set_permeability(fine_mesh_path, fp, export_permfield=True, update_permfield=True)
    set_weights_nodes(fp, update=True)

    if fp.verify_name_in_data_names('nodes_org'):
        pass
    else:
        nodes_org = fp.get_nodes_org_from_faces_of_nodes_object()
        fp.insert_or_update_data(nodes_org)

    # get_xi_params_ds_flux(fp, update=True)

    # fp.insert_or_update_data({
    #     'xi_params': fp['xi_params_ds']
    # })

    # get_lpew2_weights(fp, update=True)

    fp.backup_data('xi_params', 'xi_params_backup')
    fp.export_data()

def load_or_update_initial_loop(
        load: bool, 
        fp: MeshProperty, 
        fine_mesh_path: str,
        pressure: np.ndarray,
        newS: np.ndarray,
        vpi: float,
        cumulative_oil: float,
        cumulative_water: float,
        relative_perm: BrooksAndCorey,
        biphasic_mobility: BiphasicMobility,
        saturation: np.ndarray,
        bc: BoundaryConditions,
        lsds: LsdsFluxCalculation,
        dt: float,
        porosity: np.ndarray,
        total_area_reservoir: float,
        mesh_data: MeshData,
        simulation_data: SimulationData,
        loop: int,
        cfl
    ):
    
    if load is False:
        initial_funcs(fp, fine_mesh_path)
        pressure[:], newS[:], vpi, cumulative_oil, cumulative_water, faces_flux, water_faces_flux, water_flux, oil_flux = initial_loop(
            relative_perm,
            biphasic_mobility,
            saturation,
            fp,
            bc,
            lsds,
            dt,
            porosity,
            total_area_reservoir,
            vpi,
            cumulative_oil,
            cumulative_water,
            cfl
        )
        mesh_data.insert_tag_data('pressure', pressure, 'faces')
        mesh_data.insert_tag_data('faces_flux', faces_flux, 'faces')
        mesh_data.insert_tag_data('water_faces_flux', water_faces_flux, 'faces')
        mesh_data.insert_tag_data('saturation', saturation, 'faces')
        mesh_data.export_all_elements_type_to_vtk('pressure_faces_' + str(loop), 'faces')
        saturation[:] = newS
        simulation_data.insert_or_update_data({
            'all_loops': np.array([0]),
            'all_vpi': np.array([0.0]),
            'all_cumulative_oil': np.array([0.0]),
            'all_cumulative_water': np.array([0.0]),
            'pressure_' + str(loop): pressure,
            'saturation_' + str(loop): saturation,
            'water_flux': np.array([water_flux]),
            'oil_flux': np.array([oil_flux])
        })
        saturation[:] = newS
    elif load is True:
        # import pdb; pdb.set_trace()
        simulation_data.load_data()
        loop = simulation_data['all_loops'][-1]
        vpi = simulation_data['all_vpi'][-1]
        cumulative_oil = simulation_data['all_cumulative_oil'][-1]
        cumulative_water = simulation_data['all_cumulative_water'][-1]
        saturation[:] = simulation_data['saturation_' + str(loop)]
        pressure[:] = simulation_data['pressure_' + str(loop)]
    
    return loop, cumulative_oil, cumulative_water, vpi

def update_while_loop(
        loop_intervals: int,
        loop: int,
        pressure: np.ndarray,
        newS: np.ndarray,
        vpi: float,
        cumulative_oil: float,
        cumulative_water: float,
        relative_perm,
        biphasic_mobility,
        saturation: np.ndarray,
        fp: MeshProperty,
        bc: BoundaryConditions,
        lsds: LsdsFluxCalculation,
        porosity: np.ndarray,
        total_area_reservoir: float,
        saturation_plot: np.ndarray,
        simulation_data: SimulationData,
        mesh_data: MeshData,
        cfl: float,
        **kwargs
):
    
    for i in range(loop_intervals):
        loop += 1
        pressure[:], newS[:], vpi, cumulative_oil, cumulative_water, faces_flux, water_faces_flux, dt, water_flux, oil_flux = while_loop(
            relative_perm,
            biphasic_mobility,
            saturation,
            fp,
            bc,
            lsds,
            porosity,
            vpi,
            cumulative_oil,
            cumulative_water,
            total_area_reservoir,
            cfl=cfl
        )
        saturation_plot[:] = saturation
        saturation[:] = newS
        
        print()
        print('##########################')
        print(f'VPI: {vpi}')
        print(f'Cum oil: {cumulative_oil}')
        print(f'Cum water: {cumulative_water}')
        print(f'Loop: {loop}')
        print(f'Dt: {dt}')
        print('##########################')
        print()
    
    update_data(
        simulation_data,
        vpi,
        cumulative_oil,
        cumulative_water,
        loop,
        pressure,
        saturation,
        water_flux,
        oil_flux
    )
    
    mesh_data.insert_tag_data('pressure', pressure, 'faces')
    mesh_data.insert_tag_data('faces_flux', faces_flux, 'faces')
    mesh_data.insert_tag_data('water_faces_flux', water_faces_flux, 'faces')
    mesh_data.insert_tag_data('saturation', saturation_plot, 'faces')
    mesh_data.export_all_elements_type_to_vtk('pressure_faces_' + str(loop), 'faces')

    return loop, cumulative_oil, cumulative_water, vpi


def run5():

    dt = 0.00005
    max_vpi = 1.3
    loop = 0
    max_loop = np.inf
    load = False
    loop_intervals = 5
    cfl = 0.9

    cumulative_oil = 0.0
    cumulative_water = 0.0
    vpi = 0.0

    relative_perm = BrooksAndCorey(Sor=0.0, Swc=0.0, debug=True)
    biphasic_mobility = BiphasicMobility(mio=4)
    lsds = LsdsFluxCalculation()
    simulation_data = SimulationData('biphasic_ameba_finescale4')
    simulation_data.insert_or_update_data({'label': np.array(['finescale4'])})

    fp, fine_mesh_path = get_properties()
    bc = set_boundary_conditions(fp)

    porosity = np.repeat(0.2, len(fp['faces']))
    total_area_reservoir = porosity.dot(fp['areas'])
    saturation: np.ndarray = bc['initial_saturation']['value'].copy() 
    fp.insert_or_update_data({'sat_for_weight': saturation.copy()})
    pressure = np.repeat(0.0, fp['faces'].shape[0])
    newS = saturation.copy()
    saturation_plot = saturation.copy()

    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.create_tag('pressure')
    mesh_data.create_tag('faces_flux')
    mesh_data.create_tag('water_faces_flux')
    mesh_data.create_tag('saturation')

    loop, cumulative_oil, cumulative_water, vpi = load_or_update_initial_loop(
        load,
        fp,
        fine_mesh_path,
        pressure,
        newS,
        vpi,
        cumulative_oil,
        cumulative_water,
        relative_perm,
        biphasic_mobility,
        saturation,
        bc,
        lsds,
        dt,
        porosity,
        total_area_reservoir,
        mesh_data,
        simulation_data,
        loop,
        cfl
    )

    import pdb; pdb.set_trace()

    while vpi < max_vpi and loop < max_loop:
        
        loop, cumulative_oil, cumulative_water, vpi = update_while_loop(
            loop_intervals,
            loop,
            pressure,
            newS,
            vpi,
            cumulative_oil,
            cumulative_water,
            relative_perm,
            biphasic_mobility,
            saturation,
            fp,
            bc,
            lsds,
            porosity,
            total_area_reservoir,
            saturation_plot,
            simulation_data,
            mesh_data,
            cfl
        )

    import pdb; pdb.set_trace()


    # nodes_org, faces_of_nodes_org, n_nodes_org = fp.get_internal_nodes_org_from_faces_of_nodes_object()






    import pdb; pdb.set_trace()




    pass

