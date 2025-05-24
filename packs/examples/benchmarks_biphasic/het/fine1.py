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

    fine_mesh_properties_name = 'finescale1_het'
    fine_mesh_path = defpaths.het_fine1_mesh

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)

    return fine_properties, fine_mesh_path

def set_boundary_conditions(fine_properties: MeshProperty):
    bc = BoundaryConditions()

    initial_saturation = np.zeros(fine_properties['faces'].shape[0])
    faces_of_nodes = fine_properties['faces_of_nodes']
    nodes_centroids = fine_properties['nodes_centroids']
    nodes = fine_properties['nodes']

    xmin, ymin = nodes_centroids.min(axis=0)
    xmax, ymax = nodes_centroids.max(axis=0)

    c_q0 = np.array([xmin, ymin])
    c_p0 = np.array([xmax, ymax])

    dists1 = np.linalg.norm(nodes_centroids - c_q0, axis=1)
    dists2 = np.linalg.norm(nodes_centroids - c_p0, axis=1)

    node_q0 = nodes[dists1 <= dists1.min()]
    node_p0 = nodes[dists2 <= dists2.min()]

    faces_q0 = faces_of_nodes[node_q0[0]]
    faces_p0 = faces_of_nodes[node_p0[0]]

    faces_pressure = faces_p0
    pressure_presc = 0*np.array([1.0, 1.0])

    faces_neumann = faces_q0
    areas_faces_neumann = fine_properties['areas'][faces_neumann]
    neummann_presc_faces = 1.0*areas_faces_neumann/areas_faces_neumann.sum()

    bc.set_boundary('dirichlet_volumes', faces_pressure, pressure_presc)

    bc.set_boundary('neumann_volumes', faces_neumann, neummann_presc_faces)

    walls_edges = fine_properties.boundary_edges

    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    fine_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    bc.set_boundary('water_saturation_volumes', faces_neumann, np.array([1.0, 1.0]))
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

def set_permeability(fine_mesh_path, fine_properties: MeshProperty, simulation_data: SimulationData, export_permfield=True, update_permfield=True, **kwargs):

    x_points = fine_properties['faces_centroids'][:, 0]
    y_points = fine_properties['faces_centroids'][:, 1]
    faces = fine_properties['faces']

    factor = 2*np.cos(6*np.pi*x_points)*np.cos(6*np.pi*y_points)
    # perm = np.exp(factor)
    perm = np.power(5, factor)

    tag_preprocess = 'permeability'
    if fine_properties.verify_name_in_data_names(tag_preprocess) and update_permfield is False:
        return

    permeability = np.zeros((faces.shape[0], 2, 2))
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
    max_vpi = 1.3
    loop = 0
    max_loop = np.inf
    load = True
    loop_intervals = 5
    cfl = 0.9

    cumulative_oil = 0.0
    cumulative_water = 0.0
    vpi = 0.0

    relative_perm = BrooksAndCorey(Sor=0.0, Swc=0.0)
    biphasic_mobility = BiphasicMobility(mio=4)
    lsds = LsdsFluxCalculation()
    simulation_data = SimulationData('biphasic_het1_finescale')
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

