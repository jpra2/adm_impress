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

def get_properties():
    fine_mesh_properties_name = 'finescale_barreira'
    fine_mesh_path = defpaths.barreira_mesh

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)

    return fine_properties, fine_mesh_path

def set_boundary_conditions_tri(fine_properties: MeshProperty):
    bc = BoundaryConditions()

    initial_saturation = np.zeros(fine_properties['faces'].shape[0])

    pl_in = fine_properties['physical_line_201']
    pl_out = fine_properties['physical_line_202']
    pl_wall = fine_properties['physical_line_203']

    p0 = 5e5
    p1 = 1e5

    nodes_p0 = np.unique(fine_properties['nodes_of_edges'][pl_in].flatten())
    nodes_p1 = np.unique(fine_properties['nodes_of_edges'][pl_out].flatten())

    nodes_pressure = np.concatenate([nodes_p0, nodes_p1])
    nodes_pressure_value = np.concatenate([np.repeat(p0, len(nodes_p0)), np.repeat(p1, len(nodes_p1))])

    bc.set_boundary('dirichlet_nodes', nodes_pressure, nodes_pressure_value)
    bc.set_boundary('water_saturation_edges', pl_in, np.repeat(1.0, len(pl_in)))

    walls_edges = pl_wall
    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    fine_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    bc.set_boundary('initial_saturation', fine_properties['faces'], initial_saturation)

    bc.update_zero_bcs()

    return bc

def set_permeability_barreira(fine_mesh_path, fine_properties: MeshProperty, export_permfield=True, update_permfield=True, **kwargs):
    

    faces1 = fine_properties['physical_triangle_1']
    faces2 = fine_properties['physical_triangle_2']
    faces3 = fine_properties['physical_triangle_3']

    k1 = np.eye(2)*1e-10
    k2 = np.eye(2)*1e-4

    # k1 = np.eye(2)*1
    # k2 = np.eye(2)*1e6

    tag_preprocess = 'permeability'
    if fine_properties.verify_name_in_data_names(tag_preprocess) and update_permfield is False:
        return
    
    faces_k1 = np.concatenate([fine_properties['physical_triangle_2'], fine_properties['physical_triangle_3']])
    faces_k2 = fine_properties['physical_triangle_1']

    faces = fine_properties['faces']

    permeability = np.zeros((faces.shape[0], 2, 2))
    permeability[faces_k1] = k1
    permeability[faces_k2] = k2
    
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

    pass

def initial_funcs(
        fp: MeshProperty,
        fine_mesh_path: str
):
    
    set_permeability_barreira(fine_mesh_path, fp, export_permfield=True, update_permfield=True)
    # set_weights_nodes(fp, update=True)
    fine_properties = fp

    lsds = LsdsFluxCalculation()
    fine_properties.insert_or_update_data(
        lsds.get_all_edges_flux_params(**fine_properties.get_all_data())
    )

    fine_properties.insert_or_update_data({
        'nodes_to_calculate': fine_properties.nodes
    })

    weights = get_gls_nodes_weights(**fine_properties.get_all_data())
    fine_properties.insert_or_update_data(weights)

    fine_properties.export_data()

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
        loop: int
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
            cumulative_water
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
    loop_intervals = 1
    cfl = 0.9

    cumulative_oil = 0.0
    cumulative_water = 0.0
    vpi = 0.0

    relative_perm = BrooksAndCorey(Sor=0.1, Swc=0.1)
    biphasic_mobility = BiphasicMobility(miw=0.001, mio=0.004)
    lsds = LsdsFluxCalculation()
    simulation_data = SimulationData('biphasic_barreira')

    fp, fine_mesh_path = get_properties()
    bc = set_boundary_conditions_tri(fp)

    porosity = np.repeat(0.2, len(fp['faces']))
    porosity[fp['physical_triangle_1']] = 0.35
    print(len(fp['faces']))
    import pdb; pdb.set_trace()

    total_area_reservoir = porosity.dot(fp['areas'])
    saturation: np.ndarray = bc['initial_saturation']['value'].copy() 
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
        loop
    )

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





    

    import pdb; pdb.set_trace()