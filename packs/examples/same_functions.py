from packs.manager import MeshProperty, MeshData, BoundaryConditions
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.manager.generic_data import PrimalCoarseData
from packs.multiscale.unstructured.create_primal_dual.primal_coarse_volumes_2d import get_coarse_structure, load_coarse_structure
from packs import defnames
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.multiscale.unstructured.operators.prolongation.ams import Unstructured2DAmsOperator
from packs.utils import utils_old
from packs.manager.predef_names import TimeProfile

import numpy as np
from typing import Sequence
from shapely import geometry
from scipy.sparse.linalg import spsolve
import scipy.sparse as sp
import os
import shutil
import copy
from functools import reduce
import time
import timeit


def get_perm_diag(value):
    return np.array([[value, 0], [0, value]])

def set_permeability_brazil(fine_mesh_path, fine_properties: MeshProperty, typek='barrier', export_permfield=True, update_permfield=True) -> None:
    typeks = ['barrier', 'channel']
    assert typek in typeks

    tag_preprocess = 'permeability'
    if fine_properties.verify_name_in_data_names(tag_preprocess) and update_permfield is False:
        return

    k1 = 1.0
    k2 = 1e-3
    k3 = 1e3

    faces_in_losangle = fine_properties['faces_in_losangle']
    faces = fine_properties['faces']

    permeability = np.zeros((faces.shape[0], 2, 2))
    permeability[:] = get_perm_diag(k1)
    if typek == typeks[0]:
        permeability[faces_in_losangle] = permeability[faces_in_losangle]*k2
    else:
        permeability[faces_in_losangle] = permeability[faces_in_losangle]*k3
    
    fine_properties.insert_or_update_data({
        tag_preprocess: permeability
    })

    if export_permfield:
        mesh_data = MeshData(mesh_path=fine_mesh_path)
        mesh_data.create_tag('permeability')
        mesh_data.insert_tag_data(
            'permeability',
            fine_properties['permeability'][:, 1, 1],
            elements_type='faces'
        )
        mesh_data.export_all_elements_type_to_vtk('permfield', element_type='faces')

def define_faces_in_losangle(fine_properties: MeshProperty) -> None:
    tag_preprocess = 'faces_in_losangle'
    if fine_properties.verify_name_in_data_names(tag_preprocess):
        return
    
    Lx = 1.5
    Ly = 1

    d1 = 1.2
    d2 = 0.75

    x1 = (Lx-d1)/2
    x2 = x1 + d1/2
    x3 = x1 + d1
    x4 = x2

    y1 = Ly/2
    y2 = (Ly-d2)/2
    y3 = y1
    y4 = y2 + d2

    losangle = geometry.Polygon([
        (x1, y1),
        (x2, y2),
        (x3, y3),
        (x4, y4)
    ])

    poly = losangle
    faces_centroids = fine_properties['faces_centroids']
    points_list = geometry.MultiPoint(faces_centroids)
    test = np.array([poly.contains(i) for i in points_list.geoms])
    
    faces_in_losangle = fine_properties['faces'][test]

    fine_properties.insert_or_update_data({
        tag_preprocess: faces_in_losangle
    })

def define_coarse_structure(fine_mesh_properties: MeshProperty, lsds: LsdsFluxCalculation, level=1, update=True) -> Sequence[PrimalCoarseData]:
    fp = fine_mesh_properties
    if update is True:
        coarse_struct = get_coarse_structure(
            1,
            fp[defnames.get_primal_id_name_by_level(level)],
            fp['faces'],
            fp['adjacencies'],
            fp['edges'],
            fp['nodes_of_edges'],
            fp['bool_boundary_edges'],
            fp['bool_boundary_nodes'],
            fp['nodes'],
            fp['nodes_weights'],
            fp['nodes_of_nodes'],
            fp['edges_of_nodes'],
            fp['faces_of_nodes'],
            fp['nodes_centroids'],
            fp['faces_centroids'],
            fp['permeability'],
            fp['unitary_normal_edges'],
            fp[defnames.get_dual_id_name_by_level(1)],
            fp.edges_dim,
            lsds
        )
    else:
        coarse_struct = load_coarse_structure(
            level,
            fp[defnames.get_primal_id_name_by_level(level)]
        )
    
    return coarse_struct

def set_weights_nodes_cstruct(cstruct: PrimalCoarseData):

    fine_properties = cstruct
    weights = get_gls_nodes_weights(**fine_properties.get_all_data())
    return weights

def update_xi_params(xi_params, total_mobility_edges):
    xi_params_new = xi_params.copy()
    xi_params_new[:] = xi_params*total_mobility_edges[:, np.newaxis]
    return xi_params_new

def set_fine_transmissibility_biphasic(fine_mesh_properties: MeshProperty, bc: BoundaryConditions, lsds: LsdsFluxCalculation, only_internal_nodes=False):
    resp = lsds.mount_transmissibility_matrix(
        bc,
        **fine_mesh_properties.get_all_data(),
        only_internal_nodes=only_internal_nodes
    )
    return resp


def set_fine_transmissibility_biphasic_local(fine_mesh_properties: MeshProperty, bc: BoundaryConditions, lsds: LsdsFluxCalculation, only_internal_nodes=False):
    resp = lsds.mount_transmissibility_matrix(
        bc,
        **fine_mesh_properties.get_all_data(),
        only_internal_nodes=only_internal_nodes
    )
    return resp


def update_fine_flux(
        coarse_struct: Sequence[PrimalCoarseData],
        edges_flux: np.ndarray,
        ms_pressure: np.ndarray,
        total_mobility_edges: np.ndarray,
        lsds: LsdsFluxCalculation,
        nodes_pressure: np.ndarray,
        global_bc: BoundaryConditions,
        global_nodes_of_edges: np.ndarray,
        edges_dim: np.ndarray,
        fine_ids: np.ndarray,
        saturation: np.ndarray
    ):

    global_dirichlet_nodes = global_bc['dirichlet_nodes']['id']

    for cstruct in coarse_struct:
        local_edges_flux, edges_update = _update_fine_flux_aux(
            cstruct,
            edges_flux,
            ms_pressure,
            total_mobility_edges,
            lsds,
            nodes_pressure,
            global_bc,
            global_nodes_of_edges,
            edges_dim,
            global_dirichlet_nodes,
            fine_ids,
            saturation
        )

        edges_flux[edges_update] = local_edges_flux

def get_local_nodes_to_calculate(local_boundary_edges_flux, saturation, cstruct: PrimalCoarseData):
    
    bool_boundary_edges = cstruct['bool_boundary_edges']
    test_faces = np.isin(cstruct['faces'], cstruct['adjacencies'][bool_boundary_edges, 0])
    local_boundary_faces = cstruct['faces'][test_faces]
    global_boundary_faces = cstruct['map_faces'][test_faces]
    
    boundary_edges_flux = local_boundary_edges_flux
    saturation_local_boundary_faces = saturation[global_boundary_faces]
    
    local_saturation_test = cstruct['local_saturation_test']
    local_boundary_edge_flux_test = cstruct['local_boundary_edge_flux_for_test']
    
    if np.all(local_saturation_test == True): 
        ## se for primeiro loop
        cstruct.insert_or_update_data({
            'local_boundary_edge_flux_for_test': copy.deepcopy(boundary_edges_flux),
            'local_saturation_test': copy.deepcopy(saturation_local_boundary_faces)
        })
        nodes_to_calculate = cstruct['nodes_to_calculate'] ## mantem todos os nos do contorno do primal
    else:
        delta_sat = np.absolute(local_saturation_test - saturation_local_boundary_faces)
        test_sat = delta_sat >= cstruct['max_delta_sat']
        if test_sat.sum() > 0:
            local_saturation_test[test_sat] = saturation_local_boundary_faces[test_sat]
            cstruct.insert_or_update_data({
                'local_saturation_test': local_saturation_test
            })
        
        delta_flux_percent = cstruct['percent_var_flux']*np.absolute(local_boundary_edge_flux_test)
        delta_flux = np.absolute(local_boundary_edge_flux_test - boundary_edges_flux)
        test_max_flux = delta_flux >= delta_flux_percent
        if test_max_flux.sum() > 0:
            local_boundary_edge_flux_test[test_max_flux] = boundary_edges_flux[test_max_flux]
            cstruct.insert_or_update_data({
                'local_boundary_edge_flux_for_test': local_boundary_edge_flux_test
            })
        
        test_edges_to_update_nodes = np.isin(cstruct['adjacencies'][bool_boundary_edges, 0], local_boundary_faces[test_sat]) | test_max_flux
        if test_edges_to_update_nodes.sum() > 0:
            nodes_of_boundary_edges = cstruct['nodes_of_edges'][bool_boundary_edges]
            nodes_to_calculate = np.unique(nodes_of_boundary_edges[test_edges_to_update_nodes].flatten())
        else:
            nodes_to_calculate = np.array([])
        
    return nodes_to_calculate.astype(np.int64)

def _update_fine_flux_aux(
        cstruct: PrimalCoarseData,
        edges_flux: np.ndarray,
        ms_pressure: np.ndarray,
        total_mobility_edges: np.ndarray,
        lsds: LsdsFluxCalculation,
        nodes_pressure: np.ndarray,
        global_bc: BoundaryConditions,
        global_nodes_of_edges: np.ndarray,
        edges_dim: np.ndarray,
        global_dirichlet_nodes: np.ndarray,
        fine_ids: np.ndarray,
        saturation: np.ndarray
):
    
    global_edges = cstruct['map_edges']
    bool_boundary_edges = cstruct['bool_boundary_edges']
    # local_flux_presc = edges_flux[global_edges[bool_boundary_edges]]/(edges_dim[global_edges[bool_boundary_edges]])
    local_flux_presc = edges_flux[global_edges].copy()
    local_flux_presc[cstruct['other_side_flux']] *= -1
    local_flux_presc = -1*local_flux_presc[bool_boundary_edges]

    global_dirichlet_faces = np.intersect1d(global_bc['dirichlet_volumes']['id'], cstruct['map_faces'])
    global_neumman_faces = np.intersect1d(global_bc['neumann_volumes']['id'], cstruct['map_faces'])

    # local_flux_presc = edges_flux[global_edges[bool_boundary_edges]]
    neumann_edges = cstruct['edges'][bool_boundary_edges]
    bc = BoundaryConditions()

    # ##################
    # ## verificar edges com pressao prescrita:
    # nodes_of_global_edges = global_nodes_of_edges[global_edges]
    # test1 = np.isin(nodes_of_global_edges[:,0], global_dirichlet_nodes)
    # test2 = np.isin(nodes_of_global_edges[:,1], global_dirichlet_nodes)
    # test3 = test1 & test2
    # edges_to_remove = global_edges[test3]
    # test4 = np.isin(global_edges[bool_boundary_edges], edges_to_remove)
    # test4 = ~test4
    # #######################3
    test4 = np.full(neumann_edges.shape[0], True, dtype=bool)

    ###################
    ## verificar nos com pressao prescrita
    with_nodes_pressure = False
    # local_global_dirichlet_nodes = np.intersect1d(global_dirichlet_nodes, cstruct['map_nodes'][cstruct['bool_boundary_nodes']])
    # if local_global_dirichlet_nodes.shape[0] > 0:
    #     with_nodes_pressure = True
    #     my_nodes = cstruct['map_nodes']
    #     local_dirichlet_nodes = cstruct['nodes'][np.isin(my_nodes, local_global_dirichlet_nodes)]
    #     test5 = np.isin(global_dirichlet_nodes, local_global_dirichlet_nodes)
    #     bc.set_boundary('dirichlet_nodes', local_dirichlet_nodes, global_bc['dirichlet_nodes']['value'][test5])
    #     bc.set_boundary('dirichlet_volumes', np.array([]), np.array([]))
    ##############

    ##########
    # ## remover edges com fluxo prescrito e modificar para volumes com fluxo prescrito
    # local_adjacencies = cstruct['adjacencies']
    # neumann_volumes = local_adjacencies[neumann_edges[test4], 0]
    # neumann_value = local_flux_presc[test4]
    # bc.set_boundary('neumann_volumes', neumann_volumes, neumann_value)
    # bc.set_boundary('neumann_edges', np.array([]), np.array([]))
    # #############

    bc.set_boundary('neumann_edges', neumann_edges[test4], local_flux_presc[test4])
    # bc.set_boundary('neumann_volumes', np.array([]), np.array([]))
    # bc.set_boundary('neumann_edges', np.array([]), np.array([]))

    # if with_nodes_pressure is True:
    #     pass
    if global_dirichlet_faces.shape[0] > 0:
        all_values = global_bc['dirichlet_volumes']['value']
        all_gids =  global_bc['dirichlet_volumes']['id']
        values = []
        lids = []
        for i in global_dirichlet_faces:
            values.append(all_values[all_gids==i][0])
            lids.append(cstruct['faces'][cstruct['map_faces']==i][0])
        values = np.array(values)
        lids = np.array(lids)
        bc.set_boundary('dirichlet_volumes', lids, values)

    else:
        local_vertice = cstruct['faces'][cstruct['dual_id']==defnames.dual_ids('vertice_id')]
        global_local_vertice = cstruct['map_faces'][cstruct['dual_id']==defnames.dual_ids('vertice_id')]
        bc.set_boundary('dirichlet_volumes', local_vertice, ms_pressure[global_local_vertice])
        # bc.set_boundary('dirichlet_volumes', np.array([]), np.array([]))
        # bc.set_boundary('dirichlet_nodes', np.array([]), np.array([]))
    
    if global_neumman_faces.shape[0] > 0:
        all_values = global_bc['neumann_volumes']['value']
        all_gids =  global_bc['neumann_volumes']['id']
        values = []
        lids = []
        for i in global_neumman_faces:
            values.append(all_values[all_gids==i][0])
            lids.append(cstruct['faces'][cstruct['map_faces']==i][0])
        values = np.array(values)
        lids = np.array(lids)
        bc.set_boundary('neumann_volumes', lids, values)

    bc.update_zero_bcs()

    # bool_boundary_nodes = cstruct['bool_boundary_nodes']
    # mapbnodes = cstruct['map_nodes'][bool_boundary_nodes]
    # bnodes = cstruct['nodes'][bool_boundary_nodes]
    # bc.set_boundary('dirichlet_nodes', bnodes, nodes_pressure[mapbnodes])

    edges_multiplier = total_mobility_edges[global_edges].copy()
    ########################
    # edges_multiplier[:] = 1
    #########################
    cstruct.update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value'],
        'edges_multiplier': edges_multiplier,
    })
    
    t0 = time.perf_counter()
    nodes_to_calculate = get_local_nodes_to_calculate(
        local_flux_presc,
        saturation,
        cstruct  
    )
    
    cstruct.insert_or_update_data({'nodes_to_calculate': nodes_to_calculate})
    boundary_nodes_weights_all = set_weights_nodes_cstruct(cstruct)
    cstruct.update_all_nodes_weights(nodes_to_calculate, boundary_nodes_weights_all)
    t1 = time.perf_counter()
    TimeProfile.dt_update_boundary_nodes_weights_neumann += t1 - t0
    
    cstruct.insert_or_update_data({
        'xi_params': update_xi_params(cstruct['xi_params_backup'], total_mobility_edges[global_edges]),
    })

    t0 = time.perf_counter()
    # lt = set_fine_transmissibility_biphasic(cstruct, bc, lsds, only_internal_nodes=True)
    lt = set_fine_transmissibility_biphasic_local(cstruct, bc, lsds)
    t1 = time.perf_counter()
    TimeProfile.dt_set_local_neumann_problem += t1 - t0

    # #### segunda modificação
    # local_edges_presc_neumann = neumann_edges[test4]
    # local_presc_neumann_value = local_flux_presc[test4]
    # local_boundary_faces = cstruct['adjacencies'][local_edges_presc_neumann, 0]
    # diag = lt['transmissibility'].diagonal()
    # local_boundary_faces_diag = diag[local_boundary_faces]
    # lt['transmissibility'][local_boundary_faces] = 0
    # lt['source'][local_boundary_faces] = 0
    # lt['transmissibility'][local_boundary_faces, local_boundary_faces] = local_boundary_faces_diag
    # lt['source'][local_boundary_faces] = local_presc_neumann_value
    # lt['transmissibility'].eliminate_zeros()
    # ##################

    local_pressure = spsolve(lt['transmissibility'], lt['source'])
    
    local_edges_flux = lsds.get_edges_flux(
        bc,
        local_pressure,
        cstruct['xi_params'],
        cstruct['nodes_weights'],
        cstruct['nodes_of_edges'],
        cstruct['adjacencies'],
        cstruct['neumann_weights']
    )


    # biedges = ~bool_boundary_edges
    # local_adj = cstruct['adjacencies']
    # local_balance = np.bincount(
    #     np.concatenate([
    #         local_adj[biedges, 0],
    #         local_adj[biedges, 1],
    #         local_adj[bool_boundary_edges, 0]
    #     ]),
    #     weights=np.concatenate([
    #         local_edges_flux[biedges],
    #         -local_edges_flux[biedges],
    #         local_edges_flux[bool_boundary_edges]
    #     ])
    # )

    # print(local_balance)
    # import pdb; pdb.set_trace()

    # if cstruct['coarse_id'][0] == 1:
    #     import pdb; pdb.set_trace()

    cstruct.insert_or_update_data({
        'edges_flux': local_edges_flux
    })

    cstruct.insert_or_update_data({
        'local_pressure': local_pressure
    })

    bool_internal_edges = ~bool_boundary_edges
    # local_edges_flux[cstruct['other_side_flux']] *= -1
    # local_edges_flux[bool_boundary_edges] *= -1
    # edges_flux[cstruct['map_edges'][bool_internal_edges]] = local_edges_flux[bool_internal_edges]
    # edges_flux[cstruct['map_edges']] = local_edges_flux

    # l1 = local_edges_flux[bool_boundary_edges]
    # l2 = edges_flux[global_edges].copy()
    # l2[cstruct['other_side_flux']] *= -1
    # l2 = l2[bool_boundary_edges]

    return local_edges_flux[bool_internal_edges], global_edges[bool_internal_edges]
    
def export_op(mesh_path, OP_AMS, op_name):

    # flying_mesh_path = _create_flying_mesh(mesh_path)
    flying_mesh_path = mesh_path


    mesh_data = MeshData(mesh_path=flying_mesh_path)
    all_data = sp.find(OP_AMS)
    lines = all_data[0]
    cols = all_data[1]
    data = all_data[2]

    elements = []
    data_array = []

    cids = np.unique(cols)

    for cid in cids:
        test = cols == cid
        elements.append(lines[test])
        data_array.append(data[test])

    mesh_data.insert_array_tag_data(
        'OP',
        data_array,
        'faces',
        elements
    )

    mesh_data.export_all_elements_type_to_vtk(
        op_name,
        'faces'
    )

def define_fine_ids_from_saturation(saturation, adjacencies, internal_edges, delta_sat_lim=0.1):
    adj_sat = saturation[adjacencies[internal_edges]]
    delta_sat = np.absolute(adj_sat[:, 1] - adj_sat[:, 0])
    test = delta_sat >= delta_sat_lim

    fine_ids = adjacencies[internal_edges]
    fine_ids = np.unique(fine_ids[test])
    return fine_ids

def define_fine_ids_from_saturation_internal_nodes_org(internal_nodes_org: np.ndarray, internal_faces_of_nodes_org: np.ndarray, saturation: np.ndarray, delta_sat_min=0.1, **kwargs):
    # 'internal_nodes_org': all_nodes_org,
    # 'internal_faces_of_nodes_org': faces_of_nodes_org,
    # 'internal_n_nodes': new_n_nodes

    all_faces = []
    
    for i, nodes in enumerate(internal_nodes_org):
        faces_of_nodes = internal_faces_of_nodes_org[i]
        saturations_faces = saturation[faces_of_nodes]
        max_sat = saturations_faces.max(axis=1)
        min_sat = saturations_faces.min(axis=1)
        delta_sat = max_sat - min_sat
        test = delta_sat > delta_sat_min
        all_faces.append(np.unique(faces_of_nodes[test].flatten()))
    
    all_faces = np.unique(np.concatenate(all_faces))

    return all_faces

def get_intersect_edges(coarse_struct, fine_properties):
    fp = fine_properties
    intersect_edges = []
    for cs in coarse_struct:
        int_edges = cs['map_edges'][cs['bool_boundary_edges']]
        intersect_edges.append(int_edges)

    intersect_edges = np.unique(np.concatenate(intersect_edges))
    intersect_edges = np.intersect1d(intersect_edges, fp.internal_edges)

    return intersect_edges

def get_OR_AMS(fine_mesh_properties: MeshProperty):
    ams_prolongation = Unstructured2DAmsOperator()
    OR_AMS = ams_prolongation.get_finite_volume_restriction_operator(
        fine_mesh_properties['faces'],
        fine_mesh_properties[defnames.get_primal_id_name_by_level(1)]
    )
    return OR_AMS

def define_initial_fine_volumes(fp: MeshProperty, bc: BoundaryConditions):
    dirichlet_vols = bc['dirichlet_volumes']['id']
    values = bc['dirichlet_volumes']['id']

    cids = fp[defnames.get_primal_id_name_by_level(1)][dirichlet_vols]
    test = np.isin(fp[defnames.get_primal_id_name_by_level(1)], cids)
    fine_vols = fp['faces'][test]
    return fine_vols

def define_new_fine_levels_v1(
        fine_mesh_properties: MeshProperty,
        bc: BoundaryConditions
) -> np.ndarray:
    
    """
    Apenas as duais com prescricao na malha fina usando a dual tipo 1
    """

    faces_of_nodes = fine_mesh_properties['faces_of_nodes']
    dual_volumes = fine_mesh_properties['dual_volumes_level1']
    dual_id = fine_mesh_properties[defnames.get_dual_id_name_by_level(1)]
    
    nodes_pressure_presc = bc['dirichlet_nodes']['id']
    faces_pressure_presc = bc['dirichlet_faces']['id']
    faces_of_nodes_presc = np.unique(
        np.concatenate(faces_of_nodes[nodes_pressure_presc])
    )
    faces_of_nodes_presc = faces_of_nodes_presc[dual_id[faces_of_nodes_presc] == defnames.dual_ids('face_id')]
    
    boundary_faces = np.unique(np.concatenate([faces_of_nodes_presc, faces_pressure_presc]))
    dual_in_boundary = []
    for dual in dual_volumes:
        if np.any(np.isin(dual, boundary_faces)):
            dual_in_boundary.append(dual)
    
    dual_in_boundary = np.unique(np.concatenate(dual_in_boundary))
    return dual_in_boundary

def set_permeability_barreira_for_weights(fine_properties: MeshProperty, **kwargs):
    

    faces1 = fine_properties['physical_triangle_1']
    faces2 = fine_properties['physical_triangle_2']

    k1 = np.eye(2)*1
    k2 = np.eye(2)*1e-3

    # k1 = np.eye(2)*1
    # k2 = np.eye(2)*1e6

    tag_preprocess = 'permeability'
    
    faces_k1 = faces1
    faces_k2 = faces2

    faces = fine_properties['faces']

    permeability = np.zeros((faces.shape[0], 2, 2))
    permeability[faces_k1] = k1
    permeability[faces_k2] = k2
    
    fine_properties.insert_or_update_data({
        tag_preprocess: permeability
    })

def set_permeability_barreira_for_simulation(fine_properties: MeshProperty, **kwargs):
    

    faces1 = fine_properties['physical_triangle_1']
    faces2 = fine_properties['physical_triangle_2']

    k1 = np.eye(2)*1
    k2 = np.eye(2)*1e-3

    # k1 = np.eye(2)*1
    # k2 = np.eye(2)*1e6

    tag_preprocess = 'permeability'
    
    faces_k1 = faces1
    faces_k2 = faces2

    faces = fine_properties['faces']

    permeability = np.zeros((faces.shape[0], 2, 2))
    permeability[faces_k1] = k1
    permeability[faces_k2] = k2
    
    fine_properties.insert_or_update_data({
        tag_preprocess: permeability
    })

def calculate_data_from_dt(faces_flux: np.ndarray, injectors: np.ndarray, producers: np.ndarray, vpi: list, cum_oil: list, cum_water: list, fw_faces: np.ndarray, total_area_reservoir: float, dt: float, edges_flux: np.ndarray, bc: BoundaryConditions, fw_edges: np.ndarray):
    edges_injector = bc['edges_injector']['id']
    edges_producer = bc['edges_producer']['id']
    
    total_flux_injected = faces_flux[injectors].sum() + -1*edges_flux[edges_injector].sum()
    
    total_volume_injected = total_flux_injected*dt
    
    water_flux = (faces_flux[producers]*fw_faces[producers]).sum() + (-1*edges_flux[edges_producer]*fw_edges[edges_producer]).sum() 
    total_volume_water_produced = water_flux*dt
    fo_faces = 1-fw_faces
    fo_edges = 1-fw_edges
    oil_flux = (faces_flux[producers]*fo_faces[producers]).sum() + (-1*edges_flux[edges_producer]*fo_edges[edges_producer]).sum()
    total_volume_oil_produced = oil_flux*dt
    
    old_vpi = vpi
    
    dvpi = total_volume_injected/total_area_reservoir

    # vpi += total_volume_injected/total_area_reservoir
    vpi = dvpi + old_vpi
    cum_oil += total_volume_oil_produced
    cum_water += total_volume_water_produced
    
    max_delta = 1e-9
    if vpi - old_vpi <= max_delta:
        import pdb; pdb.set_trace()
    
    return vpi, cum_oil, cum_water, water_flux, oil_flux

def calculate_data_from_vpi(faces_flux: np.ndarray, injectors: np.ndarray, producers: np.ndarray, vpi: list, cum_oil: list, cum_water: list, fw_faces: np.ndarray, total_area_reservoir: float, dt: float, new_vpi:float, edges_flux: np.ndarray, bc: BoundaryConditions, fw_edges: np.ndarray):
    max_delta = 1e-9
    edges_injector = bc['edges_injector']['id']
    edges_producer = bc['edges_producer']['id']
    
    # total_volume_injected = faces_flux[injectors].sum()*dt
    total_flux_injected = faces_flux[injectors].sum() + -1*edges_flux[edges_injector].sum()
    # dvpi = total_volume_injected/total_volume_injected
    
    dvpi = new_vpi - vpi
    dt = dvpi*total_area_reservoir/total_flux_injected
    # total_volume_injected = faces_flux[injectors].sum()*dt  
    total_volume_injected = total_flux_injected*dt
    
    
    water_flux = (faces_flux[producers]*fw_faces[producers]).sum() + (-1*edges_flux[edges_producer]*fw_edges[edges_producer]).sum() 
    total_volume_water_produced = water_flux*dt
    fo_faces = 1-fw_faces
    fo_edges = 1-fw_edges
    oil_flux = (faces_flux[producers]*fo_faces[producers]).sum() + (-1*edges_flux[edges_producer]*fo_edges[edges_producer]).sum()
    total_volume_oil_produced = oil_flux*dt

    vpi += total_volume_injected/total_area_reservoir
    cum_oil += total_volume_oil_produced
    cum_water += total_volume_water_produced
    
    test = (new_vpi - vpi) <= max_delta
    if test == True:
        pass
    else:
        raise ValueError 
    
    return vpi, cum_oil, cum_water, water_flux, oil_flux, dt

def update_simulation_data(
    faces_flux: np.ndarray, 
    injectors: np.ndarray, 
    producers: np.ndarray, 
    vpi: list, 
    cum_oil: list, 
    cum_water: list, 
    fw_faces: np.ndarray, 
    total_area_reservoir: float, 
    dt: float, 
    vpis_to_plot,
    edges_flux: np.ndarray, 
    bc: BoundaryConditions,
    fw_edges: np.ndarray
) -> None:
    
    delta_max = 1e-9
    old_vpi = vpi
    old_cum_oil = cum_oil
    old_cum_water = cum_water
    plot_vpi = False
    
    # total_volume_injected = faces_flux[injectors].sum()*dt
    # water_flux = (faces_flux[producers]*fw_faces[producers]).sum()
    # total_volume_water_produced = water_flux*dt
    # fo_faces = 1-fw_faces
    # oil_flux = (faces_flux[producers]*fo_faces[producers]).sum()
    # total_volume_oil_produced = oil_flux*dt

    # vpi += total_volume_injected/total_area_reservoir
    # cum_oil += total_volume_oil_produced
    # cum_water += total_volume_water_produced
    
    vpi, cum_oil, cum_water, water_flux, oil_flux = calculate_data_from_dt(
        faces_flux,
        injectors,
        producers,
        old_vpi,
        old_cum_oil,
        old_cum_water,
        fw_faces,
        total_area_reservoir,
        dt,
        edges_flux,
        bc,
        fw_edges
    )
    
    dvpi = vpi - old_vpi
    if dvpi <= delta_max:
        import pdb; pdb.set_trace() 
    
    n1 = len(vpis_to_plot)
    if n1 > 0:
        calculate = True
        ids_vpis_to_plot = np.arange(n1)
        max_id = max(ids_vpis_to_plot)
        max_vpis_to_plot = max(vpis_to_plot)
        
        err1 = vpis_to_plot - old_vpi
        test1 = np.absolute(err1) <= delta_max
        if np.any(test1):
            id_target = ids_vpis_to_plot[test1] + 1
            if id_target > max_id:
                calculate = False
        elif old_vpi >= max_vpis_to_plot-delta_max:
            calculate = False
        else:
            test2 = vpis_to_plot > old_vpi
            id_target = ids_vpis_to_plot[test2][0]
        
        if calculate == True:
            vpi_target = vpis_to_plot[id_target]
            if vpi > vpi_target:
                vpi, cum_oil, cum_water, water_flux, oil_flux, dt = calculate_data_from_vpi(
                    faces_flux,
                    injectors,
                    producers,
                    old_vpi,
                    old_cum_oil,
                    old_cum_water,
                    fw_faces,
                    total_area_reservoir,
                    dt,
                    vpi_target,
                    edges_flux,
                    bc,
                    fw_edges
                )
                plot_vpi = True
 

    return vpi, cum_oil, cum_water, water_flux, oil_flux, dt, plot_vpi

def create_folders_pressure_results(pressure_folder, saturation_folder):
    
    try:
        os.makedirs(pressure_folder)
        os.makedirs(saturation_folder)
    except FileExistsError:
        shutil.rmtree(pressure_folder)
        shutil.rmtree(saturation_folder)
        os.makedirs(pressure_folder)
        os.makedirs(saturation_folder)
    

def export_ps_results(loop, pressure, saturation, pressure_folder, saturation_folder):
    ext = '.npy'
    p_str = os.path.join(pressure_folder, 'pressure_' + str(loop) + ext)
    sat_str = os.path.join(saturation_folder, 'saturation_' + str(loop) + ext)
    
    np.save(p_str, pressure)
    np.save(sat_str, saturation)

def load_ps_results(loop,  pressure_folder, saturation_folder):
    ext = '.npy'
    p_str = os.path.join(pressure_folder, 'pressure_' + str(loop) + ext)
    sat_str = os.path.join(saturation_folder, 'saturation_' + str(loop) + ext)
    
    pressure = np.load(p_str)
    saturation = np.load(sat_str)
    
    return pressure, saturation