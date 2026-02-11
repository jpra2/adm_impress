from packs.biphasic.relative_perm.brooks_and_corey import BrooksAndCorey
from packs.biphasic.mobility import BiphasicMobility
from packs.biphasic.unstructured.mobility_mesh_elements import direct_edges_mobility
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths
from packs.manager import MeshProperty, MeshData, BoundaryConditions, SimulationData
from packs.multiscale.unstructured.test.test_cross import set_weights_nodes, set_fine_transmissibility, set_permeability
# from packs.multiscale.unstructured.test.test_brazil import define_faces_in_losangle, set_permeability
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.manager.predef_names import TimeProfile

from packs.mpfa_methods.weight_interpolation.lpew import get_lpew2_weights
from packs.mpfa_methods.flux_calculation.diamond_method import get_xi_params_ds_flux
from packs.examples.same_functions import(
    update_simulation_data
)
from packs.utils import utils_old

import os
import numpy as np
from typing import Tuple
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import pint
import time



def get_properties() -> Tuple[MeshProperty, str]:
    # rel_path = os.path.join(
    #     defpaths.unstructured_coarse_test_mesh_folder,
    #     'brazil'
    # )
    # fine_mesh_path = os.path.join(rel_path, 'brazilf.msh')
    # fine_mesh_properties_name = 'brazilf' 
    # fine_mesh_path_v4 = os.path.join(rel_path, 'brazilf_v4.msh')

    rel_path = os.path.join(
        defpaths.unstructured_coarse_test_mesh_folder,
        'cross'
    )

    fine_mesh_path = os.path.join(rel_path, 'crossf.msh')
    fine_mesh_properties_name = 'crossf' 
    fine_mesh_path_v4 = os.path.join(rel_path, 'crossf.msh')

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)

    return fine_properties, fine_mesh_path

# def set_boundary_conditions(fine_properties: MeshProperty):
#     bc = BoundaryConditions()

#     nodes_centroids = fine_properties['nodes_centroids']
#     faces = fine_properties['faces']
#     faces_centroids = fine_properties['faces_centroids']

#     xmin, ymin = nodes_centroids.min(axis=0)
#     xmax, ymax = nodes_centroids.max(axis=0)

#     c_p1 = np.array([xmin, ymax])
#     c_p0 = np.array([xmax, ymin])

#     dists = np.linalg.norm(faces_centroids - c_p1, axis=1)
#     face_p1 = faces[dists <= dists.min()][0]
#     dists[:] = np.linalg.norm(faces_centroids - c_p0, axis=1)
#     face_p0 = faces[dists <= dists.min()][0]
#     faces_pressure = np.array([face_p1, face_p0])
#     pressure_presc = np.array([1.0, 0.0])

#     bc.set_boundary('dirichlet_volumes', faces_pressure, pressure_presc)
#     bc.set_boundary('dirichlet_nodes', np.array([]), np.array([]))

#     walls_edges = fine_properties['edges'][fine_properties['bool_boundary_edges']]

#     edges_values = np.repeat(0.0, walls_edges.shape[0])
#     bc.set_boundary('neumann_edges', walls_edges, edges_values)

#     fine_properties.insert_or_update_data({
#         'neumann_edges': bc['neumann_edges']['id'],
#         'neumann_edges_value': bc['neumann_edges']['value']
#     })

#     bc.set_boundary('water_saturation_volumes', np.array([face_p1]), np.array([1.0]))
#     bc.set_boundary('water_saturation_edges', np.array([]), np.array([]))

#     bc.set_boundary('injectors', np.array([face_p1]), np.array([True]))
#     bc.set_boundary('producers', np.array([face_p0]), np.array([True]))

#     bc.update_zero_bcs()

#     return bc

def set_boundary_conditions(fine_properties: MeshProperty):
    bc = BoundaryConditions()

    nodes_centroids = fine_properties['nodes_centroids']
    faces = fine_properties['faces']
    faces_centroids = fine_properties['faces_centroids']

    xmin, ymin = nodes_centroids.min(axis=0)
    xmax, ymax = nodes_centroids.max(axis=0)

    c_p1 = np.array([xmin, ymax])
    c_p0 = np.array([xmax, ymin])

    dists = np.linalg.norm(faces_centroids - c_p1, axis=1)
    face_p1 = faces[dists <= dists.min()][0]
    dists[:] = np.linalg.norm(faces_centroids - c_p0, axis=1)
    face_p0 = faces[dists <= dists.min()][0]
    
    faces_pressure = np.array([face_p0])
    pressure_presc = np.array([0.0])

    faces_neumann = np.array([face_p1])
    neummann_presc_faces = np.array([1.0])

    bc.set_boundary('dirichlet_volumes', faces_pressure, pressure_presc)

    bc.set_boundary('neumann_volumes', faces_neumann, neummann_presc_faces)

    walls_edges = fine_properties['edges'][fine_properties['bool_boundary_edges']]

    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    fine_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    bc.set_boundary('water_saturation_volumes', np.array([face_p1]), np.array([1.0]))
    bc.set_boundary('water_saturation_edges', np.array([]), np.array([]))

    bc.set_boundary('injectors', np.array([face_p1]), np.array([True]))
    bc.set_boundary('producers', np.array([face_p0]), np.array([True]))

    bc.update_zero_bcs()

    return bc

def update_xi_params(xi_params, total_mobility_edges):
    xi_params_new = xi_params.copy()
    xi_params_new[:] = xi_params*total_mobility_edges[:, np.newaxis]
    return xi_params_new
    
def update_water_faces_flux(water_faces_flux: np.ndarray, bc: BoundaryConditions, total_faces_flux: np.ndarray, relative_perm: BrooksAndCorey, biphasic_mobility: BiphasicMobility, fw_faces: np.ndarray):
    volumes_sat_prescription = bc['water_saturation_volumes']['id']
    if len(volumes_sat_prescription) > 0:
        sat_volumes_value = bc['water_saturation_volumes']['value']
        krw, kro = relative_perm.calculate(sat_volumes_value)
        mobw, mobo = biphasic_mobility.calculate(krw, kro)
        fw_presc = biphasic_mobility.get_fw(mobw, mobo)
        water_faces_flux[volumes_sat_prescription] += total_faces_flux[volumes_sat_prescription]*fw_presc
    
    producers = bc['producers']['id']
    if len(producers) > 0:
        water_flux = total_faces_flux[producers]*fw_faces[producers]
        water_faces_flux[producers] += water_flux


def calculate_dt(faces_centroids: np.ndarray, adjacencies: np.ndarray, bool_boundary_edges: np.ndarray, total_flux_edges: np.ndarray, edges_dim: np.ndarray, fw_faces: np.ndarray, saturation: np.ndarray, porosity: np.ndarray, areas: np.ndarray, faces_flux: np.ndarray, dist_centroids: np.ndarray, fw_edges: np.ndarray, edges_saturation: np.ndarray, cfl: float=0.9):
    
    dt = 1e20
    bool_internal_edges = ~bool_boundary_edges
    velocity_edges = total_flux_edges/edges_dim

    dist_internal_edges = dist_centroids[bool_internal_edges]

    v_internal_edges = np.abs(velocity_edges[bool_internal_edges])
    dfw = np.absolute(fw_faces[adjacencies[bool_internal_edges, 0]] - fw_faces[adjacencies[bool_internal_edges, 1]])
    ds = np.absolute(saturation[adjacencies[bool_internal_edges, 0]] - saturation[adjacencies[bool_internal_edges, 1]])
    adj_phi = porosity[adjacencies[bool_internal_edges]]

    test = ds != 0

    if test.sum() > 0:
        dfw = dfw[test]
        ds = ds[test]

        dfds = dfw/ds
        all_dt1: np.ndarray = cfl*dist_internal_edges[test]*adj_phi[test,0]/(v_internal_edges[test]*dfds)
        all_dt2: np.ndarray = cfl*dist_internal_edges[test]*adj_phi[test,1]/(v_internal_edges[test]*dfds)
        dt = min([all_dt1.min(), all_dt2.min()])

    velocity_bedges = np.abs(velocity_edges[bool_boundary_edges])
    dfw_bedges = np.absolute(fw_faces[adjacencies[bool_boundary_edges, 0]] - fw_edges[bool_boundary_edges])
    ds_bedges = np.absolute(saturation[adjacencies[bool_boundary_edges, 0]] - edges_saturation[bool_boundary_edges])
    phis = porosity[adjacencies[bool_boundary_edges, 0]]

    test = ds_bedges != 0

    if test.sum() > 0:
        dist_bedges = dist_centroids[bool_boundary_edges]
        dfw_bedges = dfw_bedges[test]
        ds_bedges = ds_bedges[test]

        dfds_bedges = dfw_bedges/ds_bedges
        all_dt3 = cfl*dist_bedges[test]*phis[test]/(velocity_bedges[test]*dfds_bedges)
        dt3 = all_dt3.min()

        if dt3 == 0:
            pass
        else:
            dt = min([dt, dt3])

    return dt
 

def update_saturation_dep0(water_faces_flux, areas, dt, porosity, saturation):

    ds = dt*water_faces_flux/(porosity*areas)
    newS = saturation + ds
    return newS

# def update_simulation_data(faces_flux: np.ndarray, injectors: np.ndarray, producers: np.ndarray, vpi: list, cum_oil: list, cum_water: list, fw_faces: np.ndarray, total_area_reservoir: float, dt: float) -> None:
#     total_volume_injected = faces_flux[injectors].sum()*dt
#     water_flux = (faces_flux[producers]*fw_faces[producers]).sum()
#     total_volume_water_produced = water_flux*dt
#     fo_faces = 1-fw_faces
#     oil_flux = (faces_flux[producers]*fo_faces[producers]).sum()
#     total_volume_oil_produced = oil_flux*dt

#     vpi += total_volume_injected/total_area_reservoir
#     cum_oil += total_volume_oil_produced
#     cum_water += total_volume_water_produced

#     return vpi, cum_oil, cum_water, water_flux, oil_flux

def set_fine_transmissibility_biphasic(fine_mesh_properties: MeshProperty, bc: BoundaryConditions, lsds: LsdsFluxCalculation):
    resp = lsds.mount_transmissibility_matrix(
        bc,
        **fine_mesh_properties.get_all_data()
    )
    return resp

def initial_funcs(
        fp: MeshProperty,
        fine_mesh_path: str,
        type_k: str
):
    # define_faces_in_losangle(fp)
    set_permeability(fine_mesh_path, fp, typek=type_k, export_permfield=True, update_permfield=True)
    set_weights_nodes(fp, update=True)
    nodes_org_dict = fp.get_nodes_org_from_faces_of_nodes_object()
    fp.insert_or_update_data(nodes_org_dict)
    fp.backup_data('xi_params', 'xi_params_backup')
    fp.export_data()

@utils_old.time_func
def initial_loop(
        relative_perm: BrooksAndCorey,
        biphasic_mobility: BiphasicMobility,
        saturation: np.ndarray,
        fp: MeshProperty,
        bc: BoundaryConditions,
        lsds: LsdsFluxCalculation,
        dt: float,
        porosity: np.ndarray,
        total_area_reservoir: float,
        vpi: float,
        cumulative_oil: float,
        cumulative_water: float,
        cfl: float,
        vpis_to_plot=[],
        **kwargs
):
    krw_faces, kro_faces = relative_perm.calculate(saturation)
    mobw_faces, mobo_faces = biphasic_mobility.calculate(krw_faces, kro_faces)
    total_mobility_faces = biphasic_mobility.get_total_mobility(mobw_faces, mobo_faces)
    # fp.insert_or_update_data({'faces_multiplier': total_mobility_faces})
    fw_faces = biphasic_mobility.get_fw(mobw_faces, mobo_faces)
    total_mobility_edges = direct_edges_mobility(
        total_mobility_faces,
        fp['areas'],
        fp['faces_of_nodes'],
        fp['nodes_of_edges'],
        bc,
        biphasic_mobility,
        relative_perm
    )

    fp.insert_or_update_data({
        'edges_multiplier': total_mobility_edges
    })

    fp.insert_or_update_data({
        'xi_params': update_xi_params(fp['xi_params_backup'], total_mobility_edges)
    })

    weights = get_gls_nodes_weights(**fp)
    fp.insert_or_update_data(weights)

    # get_lpew2_weights(fp, update=True)

    resp = set_fine_transmissibility_biphasic(
        fp,
        bc,
        lsds
    )

    pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])

    edges_flux = lsds.get_edges_flux(
        bc,
        pressure,
        fp['xi_params'],
        fp['nodes_weights'],
        fp['nodes_of_edges'],
        fp['adjacencies'],
        fp['neumann_weights']
    )

    edges_flux2, nodes_pressure = lsds.get_edges_flux_and_nodes_pressure(
        bc,
        pressure,
        fp['xi_params'],
        fp['nodes_weights'],
        fp['nodes_of_edges'],
        fp['adjacencies'],
        fp['neumann_weights']
    )

    # gradient_faces_dif = lsds.get_gradient_faces_dif(
    #     fp['matrix_for_gradient'],
    #     pressure,
    #     nodes_pressure,
    #     fp['nodes_of_edges'],
    #     fp['adjacencies'],
    #     fp.internal_edges,
    #     fp['Gkl']
    # )

    # estimator1 = lsds.get_estimator_1(
    #     gradient_faces_dif,
    #     fp.edges_dim,
    #     fp.internal_edges
    # )

    faces_flux = lsds.get_faces_flux(
        edges_flux,
        fp['adjacencies'],
        fp['bool_boundary_edges']
    )

    edges_saturation = edges_flux.copy()

    edges_saturation[:] = biphasic_mobility.update_edges_saturation_foum(
            saturation,
            edges_flux,
            fp['adjacencies'],
            bc,
            edges_saturation,
            fp['bool_boundary_edges']
        )
    
    # mobw_edges = direct_edges_mobility(
    #     mobw_faces,
    #     fp['areas'],
    #     fp['faces_of_nodes'],
    #     fp['nodes_of_edges'],
    #     bc,
    #     biphasic_mobility,
    #     relative_perm
    # )

    # mobo_edges = direct_edges_mobility(
    #     mobo_faces,
    #     fp['areas'],
    #     fp['faces_of_nodes'],
    #     fp['nodes_of_edges'],
    #     bc,
    #     biphasic_mobility,
    #     relative_perm
    # )



    krw_edges, kro_edges = relative_perm.calculate(edges_saturation)
    mobw_edges, mobo_edges = biphasic_mobility.calculate(krw_edges, kro_edges)
    fw_edges = biphasic_mobility.get_fw(mobw_edges, mobo_edges)

    water_edges_flux = -fw_edges*edges_flux

    water_faces_flux = lsds.get_faces_flux(
        water_edges_flux,
        fp['adjacencies'],
        fp['bool_boundary_edges']
    )

    update_water_faces_flux(
        water_faces_flux,
        bc,
        faces_flux,
        relative_perm,
        biphasic_mobility,
        fw_faces
    )

    dt = calculate_dt(
        fp['faces_centroids'],
        fp['adjacencies'],
        fp['bool_boundary_edges'],
        edges_flux,
        fp.edges_dim,
        fw_faces,
        saturation,
        porosity,
        fp['areas'],
        faces_flux,
        fp.dist_centroids,
        fw_edges,
        edges_saturation,
        cfl=cfl
    )

    newS, dt = update_saturation(water_faces_flux, fp['areas'], dt, porosity, saturation, relative_perm)
    relative_perm._test_saturations(newS)

    new_vpi, new_cumulative_oil, new_cumulative_water, water_flux, oil_flux, dt, plot_vpi = update_simulation_data(
        faces_flux,
        bc['injectors']['id'],
        bc['producers']['id'],
        vpi,
        cumulative_oil,
        cumulative_water,
        fw_faces,
        total_area_reservoir,
        dt,
        vpis_to_plot,
        edges_flux,
        bc,
        fw_edges
    )
    newS, dt = update_saturation(water_faces_flux, fp['areas'], dt, porosity, saturation, relative_perm)
    
    fp.insert_or_update_data({'dt1': np.array([dt])})
    fp.insert_or_update_data({'edges_flux0': edges_flux})
    
    return pressure, newS, new_vpi, new_cumulative_oil, new_cumulative_water, faces_flux, water_faces_flux, water_flux, oil_flux

def update_saturation(water_faces_flux, areas, dt, porosity, saturation, relative_perm: BrooksAndCorey, ratio=0.5):
    
    ds = dt*water_faces_flux/(porosity*areas)
    newS = saturation + ds
    verify = relative_perm.is_saturations_max_bound(newS)

    while verify == True:
        print('####################')
        print('dt updated')
        print('####################')
        dt = ratio*dt
        ds = dt*water_faces_flux/(porosity*areas)
        newS[:] = saturation + ds
        verify = relative_perm.is_saturations_max_bound(newS)   
    
    return newS, dt



def define_nodes_for_weight_from_delta_sat(fp: MeshProperty, saturation: np.ndarray, delta_sat_for_weight, **kwargs):

    sat_for_weight = fp['sat_for_weight']
    
    nodes_org = fp['nodes_org']
    faces_of_nodes_org = fp['faces_of_nodes_org']

    my_nodes = []
    faces_to_update_sat_for_weight = []

    for i, nodes in enumerate(nodes_org):
        faces_nodes = faces_of_nodes_org[i]
        sat1 = sat_for_weight[faces_nodes]
        sat2 = saturation[faces_nodes]
        dsat = np.absolute(sat1 - sat2)
        test = dsat >= delta_sat_for_weight
        test = np.any(test, axis=1)
        if np.any(test):
            my_nodes.append(nodes[test])
            faces_to_update_sat_for_weight.append(np.unique(faces_nodes[test].flatten()))
    if len(my_nodes) > 0:
        my_nodes = np.concatenate(my_nodes)
        faces_to_update_sat_for_weight = np.unique(np.concatenate(faces_to_update_sat_for_weight))
        sat_for_weight[faces_to_update_sat_for_weight] = saturation[faces_to_update_sat_for_weight]
        fp.insert_or_update_data({'sat_for_weight': sat_for_weight})
    else:
        my_nodes = np.array([])
    fp.insert_or_update_data({'nodes_to_calculate': my_nodes})


def update_weight_new_function(fp: MeshProperty, saturation: np.ndarray, delta_sat_for_weight=0.1, **kwargs):

    define_nodes_for_weight_from_delta_sat(fp, saturation, delta_sat_for_weight)

    weights = get_gls_nodes_weights(**fp)
    
    nodes_updated = np.unique(weights['nodes_weights']['node_id'])
    neumann_nodes_updated = weights['neumann_weights']['node_id']

    mesh_node_weight = fp['nodes_weights']
    mesh_node_neumman_weights = fp['neumann_weights']
    
    if nodes_updated.shape[0] > 0:
        test = np.isin(mesh_node_weight['node_id'], nodes_updated)
        test = ~test
        new_weight = mesh_node_weight[test].copy()
        new_weight = np.hstack([new_weight, weights['nodes_weights']])
        fp.insert_or_update_data({'nodes_weights': new_weight})
    
    if neumann_nodes_updated.shape[0] > 0:
        test2 = np.isin(mesh_node_neumman_weights['node_id'], neumann_nodes_updated)
        test2 = ~test2
        new_neumann_weight = mesh_node_neumman_weights[test2]
        new_neumann_weight = np.hstack([new_neumann_weight, weights['neumann_weights']])
        fp.insert_or_update_data({'neumann_weights': new_neumann_weight})



@utils_old.time_func_cum(export_time=True)
def while_loop(
        relative_perm: BrooksAndCorey,
        biphasic_mobility: BiphasicMobility,
        saturation: np.ndarray,
        fp: MeshProperty,
        bc: BoundaryConditions,
        lsds: LsdsFluxCalculation,
        porosity: np.ndarray,
        vpi: float,
        cumulative_oil: float,
        cumulative_water: float,
        total_area_reservoir: float,
        cfl: float,
        vpis_to_plot=[],
        **kwargs
):
    
    krw_faces, kro_faces = relative_perm.calculate(saturation)
    mobw_faces, mobo_faces = biphasic_mobility.calculate(krw_faces, kro_faces)
    total_mobility_faces = biphasic_mobility.get_total_mobility(mobw_faces, mobo_faces)
    # fp.insert_or_update_data({'faces_multiplier': total_mobility_faces})
    fw_faces = biphasic_mobility.get_fw(mobw_faces, mobo_faces)

    total_mobility_edges = direct_edges_mobility(
        total_mobility_faces,
        fp['areas'],
        fp['faces_of_nodes'],
        fp['nodes_of_edges'],
        bc,
        biphasic_mobility,
        relative_perm
    )

    fp.insert_or_update_data({
        'edges_multiplier': total_mobility_edges
    })

    fp.insert_or_update_data({
        'xi_params': update_xi_params(fp['xi_params_backup'], total_mobility_edges)
    })


    # weights = get_gls_nodes_weights(**fp)
    # fp.insert_or_update_data(weights)

    update_weight_new_function(fp, saturation)


    # get_lpew2_weights(fp, update=True)

    resp = set_fine_transmissibility_biphasic(
        fp,
        bc,
        lsds
    )

    pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])

    edges_flux = lsds.get_edges_flux(
        bc,
        pressure,
        fp['xi_params'],
        fp['nodes_weights'],
        fp['nodes_of_edges'],
        fp['adjacencies'],
        fp['neumann_weights']
    )

    faces_flux = lsds.get_faces_flux(
        edges_flux,
        fp['adjacencies'],
        fp['bool_boundary_edges']
    )
    edges_saturation = edges_flux.copy()
    edges_saturation[:] = biphasic_mobility.update_edges_saturation_foum(
        saturation,
        edges_flux,
        fp['adjacencies'],
        bc,
        edges_saturation,
        fp['bool_boundary_edges']
    )

    # mobw_edges = direct_edges_mobility(
    #     mobw_faces,
    #     fp['areas'],
    #     fp['faces_of_nodes'],
    #     fp['nodes_of_edges'],
    #     bc,
    #     biphasic_mobility,
    #     relative_perm
    # )

    # mobo_edges = direct_edges_mobility(
    #     mobo_faces,
    #     fp['areas'],
    #     fp['faces_of_nodes'],
    #     fp['nodes_of_edges'],
    #     bc,
    #     biphasic_mobility,
    #     relative_perm
    # )

    krw_edges, kro_edges = relative_perm.calculate(edges_saturation)
    mobw_edges, mobo_edges = biphasic_mobility.calculate(krw_edges, kro_edges)
    fw_edges = biphasic_mobility.get_fw(mobw_edges, mobo_edges)

    water_edges_flux = -fw_edges*edges_flux

    water_faces_flux = lsds.get_faces_flux(
        water_edges_flux,
        fp['adjacencies'],
        fp['bool_boundary_edges']
    )

    update_water_faces_flux(
        water_faces_flux,
        bc,
        faces_flux,
        relative_perm,
        biphasic_mobility,
        fw_faces
    )

    dt = calculate_dt(
        fp['faces_centroids'],
        fp['adjacencies'],
        fp['bool_boundary_edges'],
        edges_flux,
        fp.edges_dim,
        fw_faces,
        saturation,
        porosity,
        fp['areas'],
        faces_flux,
        fp.dist_centroids,
        fw_edges,
        edges_saturation,
        cfl=cfl
    )


    newS, dt = update_saturation(water_faces_flux, fp['areas'], dt, porosity, saturation, relative_perm)

    new_vpi, new_cumulative_oil, new_cumulative_water, water_flux, oil_flux, dt, plot_vpi = update_simulation_data(
        faces_flux,
        bc['injectors']['id'],
        bc['producers']['id'],
        vpi,
        cumulative_oil,
        cumulative_water,
        fw_faces,
        total_area_reservoir,
        dt,
        vpis_to_plot,
        edges_flux,
        bc,
        fw_edges
    )
    newS, dt = update_saturation(water_faces_flux, fp['areas'], dt, porosity, saturation, relative_perm)
    
    fp.insert_or_update_data({'dt1': np.array([dt])})
    fp.insert_or_update_data({'edges_flux1': edges_flux})

    return pressure, newS, new_vpi, new_cumulative_oil, new_cumulative_water, faces_flux, water_faces_flux, dt, water_flux, oil_flux, plot_vpi


def update_pressure_only(
        relative_perm: BrooksAndCorey,
        biphasic_mobility: BiphasicMobility,
        saturation: np.ndarray,
        fp: MeshProperty,
        bc: BoundaryConditions,
        lsds: LsdsFluxCalculation,
        **kwargs
):
    krw_faces, kro_faces = relative_perm.calculate(saturation)
    mobw_faces, mobo_faces = biphasic_mobility.calculate(krw_faces, kro_faces)
    total_mobility_faces = biphasic_mobility.get_total_mobility(mobw_faces, mobo_faces)
    # fp.insert_or_update_data({'faces_multiplier': total_mobility_faces})
    fw_faces = biphasic_mobility.get_fw(mobw_faces, mobo_faces)

    total_mobility_edges = direct_edges_mobility(
        total_mobility_faces,
        fp['areas'],
        fp['faces_of_nodes'],
        fp['nodes_of_edges'],
        bc,
        biphasic_mobility,
        relative_perm
    )

    fp.insert_or_update_data({
        'edges_multiplier': total_mobility_edges
    })

    fp.insert_or_update_data({
        'xi_params': update_xi_params(fp['xi_params_backup'], total_mobility_edges)
    })


    # weights = get_gls_nodes_weights(**fp)
    # fp.insert_or_update_data(weights)

    update_weight_new_function(fp, saturation)


    # get_lpew2_weights(fp, update=True)

    t0 = time.perf_counter()
    resp = set_fine_transmissibility_biphasic(
        fp,
        bc,
        lsds
    )
    t1 = time.perf_counter()
    TimeProfile.dt_set_finescale_problem = t1 - t0

    t0 = time.perf_counter()
    pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])
    t1 = time.perf_counter()
    TimeProfile.dt_solution_fs = t1 - t0

    edges_flux = lsds.get_edges_flux(
        bc,
        pressure,
        fp['xi_params'],
        fp['nodes_weights'],
        fp['nodes_of_edges'],
        fp['adjacencies'],
        fp['neumann_weights']
    )
    
    return pressure, edges_flux, fw_faces


def update_saturation_only(
    edges_flux: np.ndarray,
    relative_perm: BrooksAndCorey,
    biphasic_mobility: BiphasicMobility,
    saturation: np.ndarray,
    fp: MeshProperty,
    bc: BoundaryConditions,
    lsds: LsdsFluxCalculation,
    porosity: np.ndarray,
    total_area_reservoir: float,
    vpi: float,
    cumulative_oil: float,
    cumulative_water: float,
    dtmax: float,
    cfl: float,
    vpis_to_plot=[]
):
    
    krw_faces, kro_faces = relative_perm.calculate(saturation)
    mobw_faces, mobo_faces = biphasic_mobility.calculate(krw_faces, kro_faces)
    total_mobility_faces = biphasic_mobility.get_total_mobility(mobw_faces, mobo_faces)
    # fp.insert_or_update_data({'faces_multiplier': total_mobility_faces})
    fw_faces = biphasic_mobility.get_fw(mobw_faces, mobo_faces)
    
    edges_saturation = edges_flux.copy()
    edges_saturation[:] = biphasic_mobility.update_edges_saturation_foum(
        saturation,
        edges_flux,
        fp['adjacencies'],
        bc,
        edges_saturation,
        fp['bool_boundary_edges']
    )

    krw_edges, kro_edges = relative_perm.calculate(edges_saturation)
    mobw_edges, mobo_edges = biphasic_mobility.calculate(krw_edges, kro_edges)
    fw_edges = biphasic_mobility.get_fw(mobw_edges, mobo_edges)

    water_edges_flux = -fw_edges*edges_flux

    water_faces_flux = lsds.get_faces_flux(
        water_edges_flux,
        fp['adjacencies'],
        fp['bool_boundary_edges']
    )
    
    faces_flux = lsds.get_faces_flux(
        edges_flux,
        fp['adjacencies'],
        fp['bool_boundary_edges']
    )

    update_water_faces_flux(
        water_faces_flux,
        bc,
        faces_flux,
        relative_perm,
        biphasic_mobility,
        fw_faces
    )
    
    dt = calculate_dt(
        fp['faces_centroids'],
        fp['adjacencies'],
        fp['bool_boundary_edges'],
        edges_flux,
        fp.edges_dim,
        fw_faces,
        saturation,
        porosity,
        fp['areas'],
        faces_flux,
        fp.dist_centroids,
        fw_edges,
        edges_saturation,
        cfl=cfl
    )
    
    if dt > dtmax:
        dt = dtmax
        
    newS, dt = update_saturation(water_faces_flux, fp['areas'], dt, porosity, saturation, relative_perm)

    new_vpi, new_cumulative_oil, new_cumulative_water, water_flux, oil_flux, dt, plot_vpi = update_simulation_data(
        faces_flux,
        bc['injectors']['id'],
        bc['producers']['id'],
        vpi,
        cumulative_oil,
        cumulative_water,
        fw_faces,
        total_area_reservoir,
        dt,
        vpis_to_plot,
        edges_flux,
        bc,
        fw_edges
    )
    newS, dt = update_saturation(water_faces_flux, fp['areas'], dt, porosity, saturation, relative_perm)
    
    return newS, dt, plot_vpi, new_vpi, new_cumulative_oil, new_cumulative_water, water_flux, oil_flux, faces_flux, water_faces_flux


def update_data(
        simulation_data: SimulationData,
        vpi: float,
        cumulative_oil: float,
        cumulative_water: float,
        loop: int,
        pressure: np.ndarray,
        saturation: np.ndarray,
        water_flux,
        oil_flux
):
    
    all_loops = simulation_data['all_loops']
    all_loops = np.append(all_loops, [loop])

    all_vpi = simulation_data['all_vpi']
    all_vpi = np.append(all_vpi, [vpi])

    all_cum_oil = simulation_data[simulation_data.my_data_names[2]]
    all_cum_oil = np.append(all_cum_oil, [cumulative_oil])

    all_cum_wat = simulation_data[simulation_data.my_data_names[3]]
    all_cum_wat = np.append(all_cum_wat, [cumulative_water])
    
    all_water_flux = simulation_data[simulation_data.my_data_names[6]]
    all_water_flux = np.append(all_water_flux, [water_flux])

    all_oil_flux = simulation_data[simulation_data.my_data_names[7]]
    all_oil_flux = np.append(all_oil_flux, [oil_flux])
    
    simulation_data.insert_or_update_data({
        simulation_data.my_data_names[0]: all_loops,
        simulation_data.my_data_names[1]: all_vpi,
        simulation_data.my_data_names[2]: all_cum_oil,
        simulation_data.my_data_names[3]: all_cum_wat,
        # simulation_data.my_data_names[4] + str(loop): pressure,
        # simulation_data.my_data_names[5] + str(loop): saturation,
        simulation_data.my_data_names[6]: all_water_flux,
        simulation_data.my_data_names[7]: all_oil_flux
    })

    simulation_data.export_data()

    print()
    print(f'Simulation data updated at loop {loop}')
    print()

def plot_graph():
    fig_path = os.path.join(defpaths.results, 'cumulative_graph.svg')
    simulation_data = SimulationData('biphasic')
    simulation_data.load_data()

    all_cum_oil = simulation_data[simulation_data.my_data_names[2]]
    all_cum_wat = simulation_data[simulation_data.my_data_names[3]]
    all_vpi = simulation_data[simulation_data.my_data_names[1]]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(all_vpi, -all_cum_oil, label='Cumulative Oil')
    ax.plot(all_vpi, -all_cum_wat, label='Cumulative water')
    ax.set_xlabel('VPI')
    ax.set_ylabel('Cumulative production')
    ax.set_title('Cumulative production X VPI')

    ax.legend()
    fig.savefig(fig_path)

    



def run():

    type_k = 'channel'
    dt = 0.00005
    max_vpi = 1.3
    loop = 0
    max_loop = np.inf
    load = False
    loop_intervals = 10

    cumulative_oil = 0.0
    cumulative_water = 0.0
    vpi = 0.0
    cfl = 1.0

    relative_perm = BrooksAndCorey()
    biphasic_mobility = BiphasicMobility()
    lsds = LsdsFluxCalculation()
    simulation_data = SimulationData('biphasic_fine_mesh')

    fp, fine_mesh_path = get_properties()
    bc = set_boundary_conditions(fp)
    
    porosity = np.repeat(0.2, len(fp['faces']))
    total_area_reservoir = porosity.dot(fp['areas'])
    saturation = np.repeat(0.2, fp['faces'].shape[0])
    saturation[bc['injectors']['id']] = 0.9
    pressure = np.repeat(0.0, fp['faces'].shape[0])
    newS = saturation.copy()
    saturation_plot = saturation.copy()

    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.create_tag('pressure')
    mesh_data.create_tag('faces_flux')
    mesh_data.create_tag('water_faces_flux')
    mesh_data.create_tag('saturation')
    # mesh_data.export_all_elements_type_to_vtk('pressure_faces_' + str(loop), 'faces')
   
    if load is False:
        initial_funcs(fp, fine_mesh_path, type_k)
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
            'oil_flux': np.ndarray([oil_flux])
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

    while vpi < max_vpi and loop < max_loop:
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

        if loop % 100 == 0 and abs(cumulative_water) > 1e-4:
            import pdb; pdb.set_trace()

    import pdb; pdb.set_trace()



    

    

    




    


    




    















    print('fim')
