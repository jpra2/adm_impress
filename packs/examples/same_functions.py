from packs.manager import MeshProperty, MeshData, BoundaryConditions
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.manager.generic_data import PrimalCoarseData
from packs.multiscale.unstructured.create_primal_dual.primal_coarse_volumes_2d import get_coarse_structure, load_coarse_structure
from packs import defnames
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.multiscale.unstructured.operators.prolongation.ams import Unstructured2DAmsOperator

import numpy as np
from typing import Sequence
from shapely import geometry
from scipy.sparse.linalg import spsolve
import scipy.sparse as sp


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
        fine_ids: np.ndarray
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
            fine_ids
        )

        edges_flux[edges_update] = local_edges_flux

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
        fine_ids: np.ndarray
):
    
    global_edges = cstruct['map_edges']
    bool_boundary_edges = cstruct['bool_boundary_edges']
    # local_flux_presc = edges_flux[global_edges[bool_boundary_edges]]/(edges_dim[global_edges[bool_boundary_edges]])
    local_flux_presc = edges_flux[global_edges].copy()
    local_flux_presc[cstruct['other_side_flux']] *= -1
    local_flux_presc = -1*local_flux_presc[bool_boundary_edges]

    global_dirichlet_faces = np.intersect1d(global_bc['dirichlet_volumes']['id'], cstruct['map_faces'])

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

    if with_nodes_pressure is True:
        pass
    elif global_dirichlet_faces.shape[0] > 0:
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
    
    boundary_nodes_weights = set_weights_nodes_cstruct(cstruct)
    local_nodes_weights = cstruct['nodes_weights_internal'].copy()
    local_nodes_weights = np.hstack([local_nodes_weights, boundary_nodes_weights['nodes_weights']])

    cstruct.insert_or_update_data({
        'xi_params': update_xi_params(cstruct['xi_params_backup'], total_mobility_edges[global_edges]),
        'nodes_weights': local_nodes_weights,
        'neumann_weights': boundary_nodes_weights['neumann_weights']
    })

    # lt = set_fine_transmissibility_biphasic(cstruct, bc, lsds, only_internal_nodes=True)
    lt = set_fine_transmissibility_biphasic(cstruct, bc, lsds)

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
        local_nodes_weights,
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