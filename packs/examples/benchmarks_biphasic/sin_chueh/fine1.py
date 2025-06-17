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
from packs.utils.permfields import chueh_1

# from packs.examples.biphasic_mpfa import initial_loop, while_loop, update_data

from packs.examples.benchmarks_biphasic.ameba.mpfa_fine4 import (
    load_or_update_initial_loop,
    update_while_loop,
    create_path_mesh_data
)

import os
import numpy as np
from typing import Tuple
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

def get_properties():

    fine_mesh_properties_name = 'finescale1_sin'
    fine_mesh_path = defpaths.sin1_fine
    fine_mesh_path_v4 = defpaths.sin1_fine_v4

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name, fine_mesh_path_v4)

    return fine_properties, fine_mesh_path

def set_boundary_conditions(fine_properties: MeshProperty):
   
   
    bc = BoundaryConditions()

    initial_saturation = np.zeros(fine_properties['faces'].shape[0])

    pl_in = fine_properties['physical_line_201']
    pl_out = fine_properties['physical_line_202']
    pl_wall = fine_properties['physical_line_203']

    p0 = 1
    p1 = 0

    nodes_p0 = np.unique(fine_properties['nodes_of_edges'][pl_in].flatten())
    nodes_p1 = np.unique(fine_properties['nodes_of_edges'][pl_out].flatten())

    nodes_pressure = np.concatenate([nodes_p0, nodes_p1])
    nodes_pressure_value = np.concatenate([np.repeat(p0, len(nodes_p0)), np.repeat(p1, len(nodes_p1))])

    bc.set_boundary('dirichlet_nodes', nodes_pressure, nodes_pressure_value)
    bc.set_boundary('water_saturation_edges', pl_in, np.repeat(1.0, len(pl_in)))
    bc.set_boundary('edges_injector', pl_in, np.full(len(pl_in), True, bool))
    bc.set_boundary('edges_producer', pl_out, np.full(len(pl_in), True, bool))

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

def set_boundary_conditions_linear(fine_properties: MeshProperty):
       
    bc = BoundaryConditions()

    initial_saturation = np.zeros(fine_properties['faces'].shape[0])

    pl_in = fine_properties['physical_line_201']
    pl_out = fine_properties['physical_line_202']
    pl_wall = fine_properties['physical_line_203']

    p0 = 1
    p1 = 0

    nodes_p0 = np.unique(fine_properties['nodes_of_edges'][pl_in].flatten())
    nodes_p1 = np.unique(fine_properties['nodes_of_edges'][pl_out].flatten())

    nodes_pressure = np.concatenate([nodes_p0, nodes_p1])
    nodes_pressure_value = np.concatenate([np.repeat(p0, len(nodes_p0)), np.repeat(p1, len(nodes_p1))])

    bc.set_boundary('dirichlet_nodes', nodes_pressure, nodes_pressure_value)
    bc.set_boundary('water_saturation_edges', pl_in, np.repeat(1.0, len(pl_in)))
    bc.set_boundary('edges_injector', pl_in, np.full(len(pl_in), True, bool))
    bc.set_boundary('edges_producer', pl_out, np.full(len(pl_in), True, bool))

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

def get_R(theta):
    R = np.array([
        np.array([np.cos(theta),  np.sin(theta)]),
        np.array([-np.sin(theta), np.cos(theta)])
    ])
    return R

def set_permeability(fine_mesh_path, fine_properties: MeshProperty, simulation_data: SimulationData, export_permfield=True, update_permfield=True, **kwargs):



    tag_preprocess = 'permeability'
    if fine_properties.verify_name_in_data_names(tag_preprocess) and update_permfield is False:
        return
    
    faces = fine_properties['faces']
    faces_centroids = fine_properties['faces_centroids']

    permeability = np.zeros((faces.shape[0], 2, 2))
    perm = chueh_1(faces_centroids)
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
        name_export = os.path.join(simulation_data.name, 'permfield')
        mesh_data.export_all_elements_type_to_vtk(name_export, element_type='faces')

        mesh_data.export_all_elements_type_to_vtk('permfield', element_type='faces')

def run5():

    dt = 0.00005
    max_vpi = 0.3
    loop = 0
    max_loop = np.inf
    load = False
    loop_intervals = 20
    cfl = 0.9
    vpis_to_plot = np.linspace(0, 0.1155, 16)[1:]

    max_cum_water = 0.1

    cumulative_oil = 0.0
    cumulative_water = 0.0
    vpi = 0.0

    relative_perm = BrooksAndCorey(Sor=0.0, Swc=0.0)
    biphasic_mobility = BiphasicMobility(mio=4)
    lsds = LsdsFluxCalculation()
    simulation_data = SimulationData('biphasic_sin_chueh_fine1')
    simulation_data.insert_or_update_data({'label': np.array(['finescale1'])})
    create_path_mesh_data(simulation_data)

    fp, fine_mesh_path = get_properties()
    set_permeability(fine_mesh_path, fp, simulation_data)
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
        cfl,
        vpis_to_plot
    )

    while vpi < max_vpi and loop < max_loop and cumulative_water < max_cum_water:
        
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
            cfl,
            vpis_to_plot
        )

    import pdb; pdb.set_trace()


    # nodes_org, faces_of_nodes_org, n_nodes_org = fp.get_internal_nodes_org_from_faces_of_nodes_object()






    import pdb; pdb.set_trace()




    pass

