from packs.biphasic.relative_perm.brooks_and_corey import BrooksAndCorey
from packs.biphasic.mobility import BiphasicMobility
from packs.biphasic.unstructured.mobility_mesh_elements import direct_edges_mobility
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths, defnames
from packs.manager import MeshProperty, MeshData, BoundaryConditions, SimulationData
from packs.manager.predef_names import TimeProfile
from packs.multiscale.unstructured.test.test_cross import set_weights_nodes, set_fine_transmissibility, get_properties
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.examples.same_functions import (
    define_faces_in_losangle, 
    set_permeability_brazil as set_permeability,
    define_coarse_structure,
    update_xi_params,
    set_fine_transmissibility_biphasic,
    set_fine_transmissibility_biphasic_local,
    update_fine_flux,
    define_fine_ids_from_saturation,
    define_fine_ids_from_saturation_internal_nodes_org,
    update_simulation_data
)
from packs.examples.biphasic_mpfa import initial_funcs, set_boundary_conditions, update_saturation, calculate_dt, update_weight_new_function
from packs.examples.diss_test1 import define_new_fine_levels_v1

from packs.multiscale.unstructured.test.test_uns_ams_prolongation import export_adm_levels
from packs.utils.multiscale_methods import print_adm_interfaces_2d

from packs.multiscale.unstructured.test.test_brazil import create_primal_ids, create_dual_ids, export_primal_ids, export_dual_ids
from packs.multiscale.unstructured.test.test_brazil import get_OR_AMS
from packs.multiscale.unstructured.operators.prolongation.dual_interaction_region import create_dual_interaction_regions
from packs.multiscale.unstructured.operators.prolongation.get_op_from_amsu import update_global_op_from_amsu
from packs.utils import utils_old
from packs.adm.non_uniform import fine_level_from_alpha
from packs.fim_nu_adm.packs.processor import nu_adm_funcs
from packs.manager.generic_data import PrimalCoarseData

from packs.multiscale.unstructured.operators.precond.enhanced import Enhanced
from packs.multiscale.unstructured.operators.prolongation.msrsb_klevtsov import MsRSB

import os
import numpy as np
from typing import Tuple, Sequence
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from typing import Sequence
import time

# ### iterative ms classic
from packs.multiscale.ms_solvers.iterative_solver import iterative_ms_ilu0_bicgstab, iterative_ms_ilu0_bicgstab_fv_amsu

from packs.multiscale.unstructured.test.test_brazil import (
    f_linf_percent, 
    f_l2_error_paper_artur, 
    f_linf_percent_paper_artur, 
    f_l2_error_percent_paper_artur
)

# def get_properties():
#     rel_path = os.path.join(
#         defpaths.unstructured_coarse_test_mesh_folder,
#         'brazil'
#     )
#     fine_mesh_path = os.path.join(rel_path, 'brazilf.msh')
#     fine_mesh_properties_name = 'brazilf' 
#     fine_mesh_path_v4 = os.path.join(rel_path, 'brazilf_v4.msh')

#     coarse_mesh_path = os.path.join(rel_path, 'brazilC.msh')
#     coarse_mesh_properties_name = 'brazilC1'

#     # coarse_mesh_path = os.path.join(rel_path, 'brazilC2.msh')
#     # coarse_mesh_properties_name = 'brazilC2'


#     fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)
#     coarse_properties = preprocess_mesh(coarse_mesh_path, coarse_mesh_properties_name)

#     return fine_properties, coarse_properties, fine_mesh_path, coarse_mesh_path

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
    
#     faces_pressure = np.array([face_p0])
#     pressure_presc = np.array([0.0])

#     faces_neumann = np.array([face_p1])
#     neummann_presc_faces = np.array([1.0])

#     bc.set_boundary('dirichlet_volumes', faces_pressure, pressure_presc)

#     bc.set_boundary('neumann_volumes', faces_neumann, neummann_presc_faces)

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

def calculate_dt_v1_ms(faces_centroids: np.ndarray, adjacencies: np.ndarray, bool_boundary_edges: np.ndarray, total_flux_edges: np.ndarray, edges_dim: np.ndarray, fw_faces: np.ndarray, saturation: np.ndarray, porosity: np.ndarray, cfl: float=1.0):
    bool_internal_edges = ~bool_boundary_edges
    velocity_edges = total_flux_edges/edges_dim

    dist_internal_edges = np.linalg.norm(
        faces_centroids[adjacencies[bool_internal_edges, 0]] - faces_centroids[adjacencies[bool_internal_edges, 1]],
        axis=1
    )
    v_internal_edges = np.abs(velocity_edges[bool_internal_edges])
    dfw = np.absolute(fw_faces[adjacencies[bool_internal_edges, 0]] - fw_faces[adjacencies[bool_internal_edges, 1]])
    ds = np.absolute(saturation[adjacencies[bool_internal_edges, 0]] - saturation[adjacencies[bool_internal_edges, 1]])
    adj_phi = porosity[adjacencies[bool_internal_edges]]

    test = ds != 0

    dfw = dfw[test]
    ds = ds[test]

    dfds = dfw/ds
    all_dt1: np.ndarray = cfl*dist_internal_edges[test]*adj_phi[test,0]/(v_internal_edges[test]*dfds)
    all_dt2: np.ndarray = cfl*dist_internal_edges[test]*adj_phi[test,1]/(v_internal_edges[test]*dfds)
    dt = min([all_dt1.min(), all_dt2.min()])
    return dt
    



def get_OP_matrix(fine_mesh_properties: MeshProperty):
    primal_ids = fine_mesh_properties[defnames.get_primal_id_name_by_level(1)]
    cids = np.unique(primal_ids)
    ncids = cids.shape[0]
    nfids = fine_mesh_properties['faces'].shape[0]
    OP = sp.lil_matrix((nfids, ncids))
    return OP

def get_op_amsu(fine_mesh_properties: MeshProperty, lsds: LsdsFluxCalculation, OP: sp.lil_matrix, fine_transm_without_bc: sp.csc_matrix):
    T = fine_transm_without_bc
    level_str = '_level1'

    interaction_regions = create_dual_interaction_regions(
        fine_mesh_properties[defnames.get_dual_interation_region_name_by_level(1)],
        fine_mesh_properties[defnames.vertices_selected + level_str],
        fine_mesh_properties[defnames.get_primal_id_name_by_level(1)][fine_mesh_properties[defnames.vertices_selected + level_str]],
        fine_mesh_properties[defnames.internal_dual_path + level_str],
        fine_mesh_properties[defnames.boundary_dual_interaction + level_str],
        fine_mesh_properties[defnames.dual_initial_ccs + level_str],
        global_transmissibility=T,
        global_diagonal_term=np.zeros(fine_mesh_properties['faces'].shape[0]),
        dual_id=fine_mesh_properties[defnames.get_dual_id_name_by_level(1)]
    )

    OP = update_global_op_from_amsu(interaction_regions, OP)
    return OP

# def initial_funcs(
#         fp: MeshProperty,
#         fine_mesh_path: str,
#         type_k: str
# ):
#     define_faces_in_losangle(fp)
#     set_permeability(fine_mesh_path, fp, typek=type_k, export_permfield=True, update_permfield=True)
#     set_weights_nodes(fp, update=True)
#     fp.backup_data('xi_params', 'xi_params_backup')
#     fp.export_data()


def refine_by_gradient_v0(fp: MeshProperty, pressure: np.ndarray, max_grad: float=0.1, refine=True):
    # refine = False
    ## limite para refinar considerando os primais 0.5
    if refine is True: 
        # dist_centroids = fp['dist_centroids']
        dist_centroids = fp.dist_centroids
        nodes_edges = fp['nodes_of_edges']
        faces_nodes = fp['faces_of_nodes']
        primal_id = fp[defnames.get_primal_id_name_by_level(1)]
        faces_of_faces_by_nodes = fp.faces_of_faces_by_nodes

        dp = pressure[fp['adjacencies'][fp.internal_edges]]
        dp = np.absolute(dp[:, 1] - dp[:, 0])
        grad = dp/dist_centroids[fp.internal_edges]
        test = grad > max_grad
        nodes_selected = np.unique(nodes_edges[fp.internal_edges[test]].flatten())
        faces_selected = np.unique(np.concatenate(faces_nodes[nodes_selected]))
        # faces_selected = np.unique(np.concatenate(faces_of_faces_by_nodes[faces_selected])) 

        # primal_id_selected = np.unique(primal_id[faces_selected])
        # test1 = np.isin(primal_id, primal_id_selected)
        # faces_selected = fp['faces'][test1]
    elif refine is False:
        faces_selected  = np.array([])
    
    return faces_selected

def refine_by_pressure_lim(fp: MeshProperty, pressure: np.ndarray, min_pressure=1e5, max_pressure=5e5):
    f1 = fp['faces'][pressure > max_pressure]
    f2 = fp['faces'][pressure < min_pressure]
    all_faces = np.concatenate([f1, f2])

    all_faces = np.concatenate(fp.faces_of_faces_by_nodes[all_faces])

    return all_faces

def refine_by_estimator1_delta_pressure(
        fp: MeshProperty,
        lsds: LsdsFluxCalculation,
        bc: BoundaryConditions,
        pressure: np.ndarray,
        # max_value: float=700
        max_value: float=1200,
        delta_lim = 0.01
):
    
    """Estimador baseado no delta de pressao

    Returns:
        _type_: _description_
    """
    
    pressures_inj = pressure[bc['injectors']['id']]
    pressures_prod = pressure[bc['producers']['id']]
    all_pressures = np.concatenate([pressures_inj, pressures_prod])
    
    
    pmax = all_pressures.max()
    pmin = all_pressures.min()
    deltamax = pmax - pmin
    iadjacencies = fp['adjacencies'][fp.internal_edges]
    
    delta_p_iadjacencies = np.absolute(pressure[iadjacencies[:, 1]] - pressure[iadjacencies[:, 0]])/deltamax
    test = delta_p_iadjacencies >= delta_lim
    
    faces_selected = np.unique(iadjacencies[test].flatten())
    return faces_selected
    
    
    
    

def refine_by_estimator1(
        fp: MeshProperty,
        lsds: LsdsFluxCalculation,
        bc: BoundaryConditions,
        pressure: np.ndarray,
        # max_value: float=700
        max_value: float=1200
):
    """Estimaodr baseado no delta do grad

    Args:
        fp (MeshProperty): _description_
        lsds (LsdsFluxCalculation): _description_
        bc (BoundaryConditions): _description_
        pressure (np.ndarray): _description_
        max_value (float, optional): _description_. Defaults to 1200.

    Returns:
        _type_: _description_
    """
    
    injectors = bc['injectors']['id']
    producers = bc['producers']['id']
    wells = np.concatenate([injectors, producers])

    edges_flux, nodes_pressure = lsds.get_edges_flux_and_nodes_pressure(
        bc,
        pressure,
        fp['xi_params'],
        fp['nodes_weights'],
        fp['nodes_of_edges'],
        fp['adjacencies'],
        fp['neumann_weights']
    )

    gradient_faces_dif = lsds.get_gradient_faces_dif(
        fp['matrix_for_gradient'],
        pressure,
        nodes_pressure,
        fp['nodes_of_edges'],
        fp['adjacencies'],
        fp.internal_edges,
        fp['Gkl']
    )

    estimator1 = lsds.get_estimator_1_dkl(
        gradient_faces_dif,
        fp.dkl,
        fp.internal_edges
    )

    internal_edges = fp.internal_edges
    adjacencies = fp['adjacencies']
    faces_of_faces_by_nodes = fp.faces_of_faces_by_nodes

    estimator_internal_edges = estimator1[internal_edges]
    adj_internal_edges = adjacencies[internal_edges]
    
    # indicator1 = np.isin(adj_internal_edges[:, 0], wells)
    # indicator2 = np.isin(adj_internal_edges[:, 1], wells)
    # indicator = indicator1 | indicator2
    
    # ########
    # max_value_estimator_internal_edges = estimator_internal_edges[indicator].max()
    max_value_estimator_internal_edges = estimator_internal_edges.max()
    # print(f'Value Estimator: {max_value_estimator_internal_edges}')
    # max_value_estimator_internal_edges = estimator_internal_edges[indicator].min()
    # import pdb; pdb.set_trace()
    normalized_estimator = estimator_internal_edges/max_value_estimator_internal_edges
    # max_norm_estimator = normalized_estimator.max()
    # test2 = estimator_internal_edges >= max_value
    test2 = normalized_estimator >= max_value
    # new_param_test = normalized_estimator[test2]
    # ###################################

    # test = estimator_internal_edges >= max_value
    test = test2
    # faces_to_refine = np.unique(adj_internal_edges[test].flatten())
    # faces_to_refine = np.unique(np.concatenate(faces_of_faces_by_nodes[faces_to_refine]))
    nodes_selected = np.unique(np.concatenate(fp['nodes_of_edges'][internal_edges[test]]))
    faces_to_refine = np.unique(np.concatenate(fp['faces_of_nodes'][nodes_selected]))

    faces_selected = faces_to_refine
    # primal_id_selected = np.unique(primal_id[faces_selected])
    # test1 = np.isin(primal_id, primal_id_selected)
    # faces_selected = fine_faces[test1]

    return faces_selected




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
        matrices_path,
        op_name,
        initial_fine_volumes: np.ndarray,
        alpha_lim_finescale,
        beta_lim,
        etol_msrsb,
        maxit_msrsb,
        cfl: float,
        refine_by_grad_bool: bool,
        refine_by_estimator1_bool: bool,
        max_value_grad: float,
        max_value_estimator1: float,
        vpis_to_plot: np.ndarray,
        iterative_ms: bool,
        tol_iterative: float,
        p0: np.ndarray,
        **kwargs
):
    
    it = -1
    err = -1
    
    
    # refine_by_grad_bool = False
    # refine_by_estimator1_bool = True 

    krw_faces, kro_faces = relative_perm.calculate(saturation)
    mobw_faces, mobo_faces = biphasic_mobility.calculate(krw_faces, kro_faces)
    total_mobility_faces = biphasic_mobility.get_total_mobility(mobw_faces, mobo_faces)
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

    coarse_struct = define_coarse_structure(fp, lsds, level=1)

    intersect_edges = []
    for cs in coarse_struct:
        int_edges = cs['map_edges'][cs['bool_boundary_edges']]
        intersect_edges.append(int_edges)

    intersect_edges = np.unique(np.concatenate(intersect_edges))
    intersect_edges = np.intersect1d(intersect_edges, fp.internal_edges)

    cadj_fine = fp[defnames.get_primal_id_name_by_level(1)][fp['adjacencies']]
    cadj_fine[fp['adjacencies'] == -1] = -1

    cadj_intersect = cadj_fine[intersect_edges]

    resp = set_fine_transmissibility_biphasic(
        fp,
        bc,
        lsds
    )

    
    OR = get_OR_AMS(fp)
    
    if op_name == defnames.list_op_toget[1]:
        OP = get_OP_matrix(fp)
        fine_transm_without_bc = lsds.mount_transmissibility_matrix_without_bc(**fp.get_all_data())
        OP = get_op_amsu(fp, lsds, OP, fine_transm_without_bc['transmissibility_without_bc'])
    elif op_name == defnames.list_op_toget[2]:
        level_str = defnames.level_str(1)
        enhanced = Enhanced()
        monotone_transm = enhanced.get_enhanced_matrix(resp['transmissibility'])
        msrsb = MsRSB()
        OP = msrsb.get_OP(
            faces=fp['faces'],
            T=monotone_transm,
            diagonal_term=np.zeros(resp['source'].shape[0]),
            interation_regions=fp[defnames.get_dual_interation_region_name_by_level(1)],
            interation_boundaries=fp[defnames.boundary_dual_interaction + level_str],
            vertices=fp[defnames.vertices_selected + level_str],
            dual_edges=fp['faces'][fp[defnames.get_dual_id_name_by_level(1)]==defnames.dual_ids('edge_id')],
            dual_faces=fp['faces'][fp[defnames.get_dual_id_name_by_level(1)]==defnames.dual_ids('face_id')],
            coarse_ids=fp[defnames.get_primal_id_name_by_level(1)][fp[defnames.vertices_selected + level_str]],
            OR_fv=OR,
            etol=etol_msrsb,
            maxit=maxit_msrsb                
        )
    else:
        raise NameError

    utils_old.save_matrix(matrices_path, op_name, OP)

    fine_levels = np.full(fp['faces'].shape[0], -1)
    fine_levels[initial_fine_volumes] = 0
    # fine_ids_from_saturation = define_fine_ids_from_saturation(saturation, fp['adjacencies'], fp.internal_edges)
    fine_ids_from_saturation = define_fine_ids_from_saturation_internal_nodes_org(fp['internal_nodes_org'], fp['internal_faces_of_nodes_org'], saturation)


    fine_levels[fine_ids_from_saturation] = 0

    fine_ids_from_perm = refine_from_permeability_value(fp)
    if fine_ids_from_perm.shape[0] > 0:
        fine_levels[fine_ids_from_perm] = 0

    fine_ids_from_alpha = fine_level_from_alpha.define_fine_levels_from_alpha(
        OR,
        OP,
        resp['transmissibility'],
        fp[defnames.get_primal_id_name_by_level(1)],
        alpha_lim=alpha_lim_finescale
    )

    # fine_ids_from_alpha = fine_level_from_alpha.define_fine_levels_from_alpha(
    #     OR,
    #     OP,
    #     resp['transmissibility'],
    #     alpha_lim=alpha_lim_finescale
    # )

    fine_levels[fine_ids_from_alpha] = 0

    beta_groups, beta_ind, betas = nu_adm_funcs.get_beta_groups(
        fp['faces'],
        fp[defnames.get_primal_id_name_by_level(1)],
        sp.find(OP)[0:3],
        fp['adjacencies'][fp.internal_edges],
        beta_lim=beta_lim
    )

    finescale_faces = nu_adm_funcs.get_finescale_vols(
        fp['faces'][fine_levels==0],
        fine_ids_from_alpha,
        beta_ind,
        beta_groups
    )

    # finescale_faces = finescale_faces.astype(np.int)

    fine_levels[finescale_faces] = 0
    fine_levels[fine_levels==-1] = 1

    finescale_ids = fp['faces'][fine_levels==0]

    
    LEVEL_ID_1, ADM_COARSE_ID_LEVEL_1 = nu_adm_funcs.set_adm_mesh_non_nested(
        finescale_ids,
        fine_levels,
        fp['faces'],
        fp[defnames.get_primal_id_name_by_level(1)],
        fp[defnames.get_dual_id_name_by_level(1)]
    )

    OP_adm, OR_adm = nu_adm_funcs.organize(
        fine_levels,
        sp.find(OP)[0:3],
        fp['faces'],
        fp[defnames.get_primal_id_name_by_level(1)],
        LEVEL_ID_1,
        ADM_COARSE_ID_LEVEL_1,
        fp[defnames.get_dual_id_name_by_level(1)]
    )

    T_adm = OR_adm*(resp['transmissibility']*OP_adm)
    Q_adm = OR_adm*resp['source']
    P_adm = spsolve(T_adm.tocsc(), Q_adm)
    P_prol = OP_adm*P_adm

    ##############################################################
    ## refine by delta grad
    if refine_by_grad_bool is True:
        
        fine_faces_by_grad = refine_by_gradient_v0(fp, P_prol, max_grad=max_value_grad)
        fine_faces_by_grad = nu_adm_funcs.get_finescale_vols(
            fp['faces'][fine_levels==0],
            fine_faces_by_grad,
            beta_ind,
            beta_groups
        )
        fp.insert_or_update_data({'fine_faces_by_grad': fine_faces_by_grad})

        if fine_faces_by_grad.shape[0] > 0:
            fine_levels[fine_faces_by_grad] = 0
        finescale_ids = fp['faces'][fine_levels==0]

    # LEVEL_ID_1, ADM_COARSE_ID_LEVEL_1 = nu_adm_funcs.set_adm_mesh_non_nested(
    #     finescale_ids,
    #     fine_levels,
    #     fp['faces'],
    #     fp[defnames.get_primal_id_name_by_level(1)],
    #     fp[defnames.get_dual_id_name_by_level(1)]
    # )

    # OP_adm, OR_adm = nu_adm_funcs.organize(
    #     fine_levels,
    #     sp.find(OP)[0:3],
    #     fp['faces'],
    #     fp[defnames.get_primal_id_name_by_level(1)],
    #     LEVEL_ID_1,
    #     ADM_COARSE_ID_LEVEL_1,
    #     fp[defnames.get_dual_id_name_by_level(1)]
    # )

    # T_adm = OR_adm*(resp['transmissibility']*OP_adm)
    # Q_adm = OR_adm*resp['source']
    # P_adm = spsolve(T_adm.tocsc(), Q_adm)
    # P_prol = OP_adm*P_adm
    ##############################################################

    # fine_faces_by_pressure = refine_by_pressure_lim(fp, P_prol)
    # fine_levels[fine_faces_by_pressure] = 0
    # finescale_ids = fp['faces'][fine_levels==0]

    # LEVEL_ID_1, ADM_COARSE_ID_LEVEL_1 = nu_adm_funcs.set_adm_mesh_non_nested(
    #     finescale_ids,
    #     fine_levels,
    #     fp['faces'],
    #     fp[defnames.get_primal_id_name_by_level(1)],
    #     fp[defnames.get_dual_id_name_by_level(1)]
    # )

    # OP_adm, OR_adm = nu_adm_funcs.organize(
    #     fine_levels,
    #     sp.find(OP)[0:3],
    #     fp['faces'],
    #     fp[defnames.get_primal_id_name_by_level(1)],
    #     LEVEL_ID_1,
    #     ADM_COARSE_ID_LEVEL_1,
    #     fp[defnames.get_dual_id_name_by_level(1)]
    # )

    # T_adm = OR_adm*(resp['transmissibility']*OP_adm)
    # Q_adm = OR_adm*resp['source']
    # P_adm = spsolve(T_adm.tocsc(), Q_adm)
    # P_prol = OP_adm*P_adm

    ###############################################################
    ## refine by estimator 1 (gk - gl)/edge_dim
    if refine_by_estimator1_bool is True:
        
        fine_faces_by_estimator1 = refine_by_estimator1(
            fp,
            lsds,
            bc,
            P_prol,
            max_value=max_value_estimator1
        )

        if fine_faces_by_estimator1.shape[0] > 0:
            fine_levels[fine_faces_by_estimator1] = 0
        finescale_ids = fp['faces'][fine_levels==0]

        
    ###############################################################
    
    if refine_by_grad_bool is True or refine_by_estimator1_bool is True:
        LEVEL_ID_1, ADM_COARSE_ID_LEVEL_1 = nu_adm_funcs.set_adm_mesh_non_nested(
            finescale_ids,
            fine_levels,
            fp['faces'],
            fp[defnames.get_primal_id_name_by_level(1)],
            fp[defnames.get_dual_id_name_by_level(1)]
        )

        OP_adm, OR_adm = nu_adm_funcs.organize(
            fine_levels,
            sp.find(OP)[0:3],
            fp['faces'],
            fp[defnames.get_primal_id_name_by_level(1)],
            LEVEL_ID_1,
            ADM_COARSE_ID_LEVEL_1,
            fp[defnames.get_dual_id_name_by_level(1)]
        )

        T_adm = OR_adm*(resp['transmissibility']*OP_adm)
        Q_adm = OR_adm*resp['source']
        P_adm = spsolve(T_adm.tocsc(), Q_adm)
        P_prol = OP_adm*P_adm
    
    if iterative_ms:
        P_prol[:], it, err = iterative_ms_ilu0_bicgstab(
            resp['transmissibility'],
            resp['source'],
            p0,
            # P_prol,
            OP_adm,
            OR_adm,
            epsilon=tol_iterative
        )
        
        P_prol2, it2, err2 = iterative_ms_ilu0_bicgstab(
                resp['transmissibility'],
                resp['source'],
                p0,
                # P_prol,
                OP,
                OR,
                epsilon=1e-12
            )
        
        np.save('it_nuadm.npy', it)
        np.save('it2.npy', it2)
        np.save('err_nuadm.npy', err)
        np.save('err2.npy', err2)
    
    # pressure = spsolve(resp['transmissibility'], resp['source'])
    
    # plinf = f_linf_percent_paper_artur(pressure, P_prol)/100
    # pl2 = f_l2_error_paper_artur(pressure, P_prol)
    
    # plinf2 = f_linf_percent_paper_artur(pressure, P_prol2)/100
    # pl22 = f_l2_error_paper_artur(pressure, P_prol2)
    
    edges_flux, nodes_pressure = lsds.get_edges_flux_and_nodes_pressure(
        bc,
        P_prol,
        fp['xi_params'],
        fp['nodes_weights'],
        fp['nodes_of_edges'],
        fp['adjacencies'],
        fp['neumann_weights']
    )

    fine_edges_flux = edges_flux.copy()
    fine_faces_flux = lsds.get_faces_flux(
        fine_edges_flux,
        fp['adjacencies'],
        fp['bool_boundary_edges']
    )

    fid = bc['dirichlet_volumes']['id']

    intersect_flux = edges_flux[intersect_edges]
    bflux = edges_flux[fp.boundary_edges]

    v1 = np.concatenate([intersect_flux, -intersect_flux, bflux])
    v2 = np.concatenate([cadj_intersect[:, 0], cadj_intersect[:, 1], cadj_fine[fp.boundary_edges, 0]])
    coarse_face_flux = np.bincount(v2, weights=v1)
    cff = coarse_face_flux

    update_fine_flux(
        coarse_struct,
        edges_flux,
        P_prol,
        total_mobility_edges,
        lsds,
        nodes_pressure,
        bc,
        fp['nodes_of_edges'],
        fp.edges_dim,
        finescale_ids,
        saturation
    )

    faces_flux = lsds.get_faces_flux(
        edges_flux,
        fp['adjacencies'],
        fp['bool_boundary_edges']
    )

    # edges_flux2 = lsds.get_edges_flux(
    #     bc,
    #     pressure,
    #     fp['xi_params'],
    #     fp['nodes_weights'],
    #     fp['nodes_of_edges'],
    #     fp['adjacencies'],
    #     fp['neumann_weights']
    # )

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

    intitial_fine_faces = fp['faces'][fine_levels==0]
    # intitial_fine_faces = initial_fine_volumes
    fp.insert_or_update_data({'initial_fine_faces': intitial_fine_faces})
    fp.insert_or_update_data(
        {   
            'nuadm_vols': np.array([T_adm.shape[0]]),
            'it': np.array([it]),
            'err': np.array([err])
        }
    )
    
    fp.insert_or_update_data({'dt1': np.array([dt])})
    fp.insert_or_update_data({'edges_flux0': edges_flux})

    return P_prol, newS, new_vpi, new_cumulative_oil, new_cumulative_water, faces_flux, coarse_struct, OP, OR, fine_levels, water_flux, oil_flux

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
        matrices_path: str,
        op_name: str,
        initial_fine_volumes: np.ndarray,
        alpha_lim_finescale: float,
        beta_lim: float,
        OP: sp.csc_matrix,
        OR: sp.csc_matrix,
        coarse_struct: Sequence[PrimalCoarseData],
        cfl: float,
        vpis_to_plot: np.ndarray,
        iterative_ms: bool,
        tol_iterative: float,
        p0: np.ndarray,
        **kwargs
):
    
    it = -1
    err = -1
    
    krw_faces, kro_faces = relative_perm.calculate(saturation)
    mobw_faces, mobo_faces = biphasic_mobility.calculate(krw_faces, kro_faces)
    total_mobility_faces = biphasic_mobility.get_total_mobility(mobw_faces, mobo_faces)
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

    for coarse_data in coarse_struct:
        coarse_data.get_local_nodes_weights(
            fp['nodes_weights'],
            coarse_data['map_nodes'],
            coarse_data['nodes'],
            coarse_data['map_faces'],
            coarse_data['faces'],
            coarse_data['bool_boundary_nodes'],
            fp['nodes_to_calculate']
        )

    resp = set_fine_transmissibility_biphasic(
        fp,
        bc,
        lsds
    )

    fine_levels = np.full(fp['faces'].shape[0], -1)
    fine_levels[initial_fine_volumes] = 0
    fine_levels[fp['initial_fine_faces']] = 0
    # fine_ids_from_saturation = define_fine_ids_from_saturation(saturation, fp['adjacencies'], fp.internal_edges)
    fine_ids_from_saturation = define_fine_ids_from_saturation_internal_nodes_org(fp['internal_nodes_org'], fp['internal_faces_of_nodes_org'], saturation)

    fine_levels[fine_ids_from_saturation] = 0

    # fine_ids_from_perm = refine_from_permeability_value(fp)
    # if fine_ids_from_perm.shape[0] > 0:
    #     fine_levels[fine_ids_from_perm] = 0

    fine_ids_from_alpha = fine_level_from_alpha.define_fine_levels_from_alpha(
        OR,
        OP,
        resp['transmissibility'],
        fp[defnames.get_primal_id_name_by_level(1)],
        alpha_lim=alpha_lim_finescale
    )

    fine_levels[fine_ids_from_alpha] = 0

    beta_groups, beta_ind, betas = nu_adm_funcs.get_beta_groups(
        fp['faces'],
        fp[defnames.get_primal_id_name_by_level(1)],
        sp.find(OP)[0:3],
        fp['adjacencies'][fp.internal_edges],
        beta_lim=beta_lim
    )

    finescale_faces = nu_adm_funcs.get_finescale_vols(
        fp['faces'][fine_levels==0],
        fine_ids_from_alpha,
        beta_ind,
        beta_groups
    )

    finescale_faces = finescale_faces.astype(np.int)

    fine_levels[finescale_faces] = 0
    fine_levels[fine_levels==-1] = 1

    finescale_ids = fp['faces'][fine_levels==0]

    LEVEL_ID_1, ADM_COARSE_ID_LEVEL_1 = nu_adm_funcs.set_adm_mesh_non_nested(
        finescale_ids,
        fine_levels,
        fp['faces'],
        fp[defnames.get_primal_id_name_by_level(1)],
        fp[defnames.get_dual_id_name_by_level(1)]
    )

    OP_adm, OR_adm = nu_adm_funcs.organize(
        fine_levels,
        sp.find(OP)[0:3],
        fp['faces'],
        fp[defnames.get_primal_id_name_by_level(1)],
        LEVEL_ID_1,
        ADM_COARSE_ID_LEVEL_1,
        fp[defnames.get_dual_id_name_by_level(1)]
    )
    
    if iterative_ms:
        P_prol, it, err = iterative_ms_ilu0_bicgstab(
            resp['transmissibility'],
            resp['source'],
            p0,
            # P_prol,
            OP_adm,
            OR_adm,
            epsilon=tol_iterative
        )
    
    else:
        T_adm = OR_adm*(resp['transmissibility']*OP_adm)
        Q_adm = OR_adm*resp['source']
        P_adm = spsolve(T_adm.tocsc(), Q_adm)
        P_prol = OP_adm*P_adm
    
    

    # pmin = P_prol.min()
    # pmax = P_prol.max()
    # max_grad = 0.8*(pmax - pmin)/np.sqrt(2)

    # fine_faces_by_grad = refine_by_gradient_v0(fp, P_prol, max_grad=max_grad)
    # # fine_faces_by_grad = nu_adm_funcs.get_finescale_vols(
    # #     fp['faces'][fine_levels==0],
    # #     fine_faces_by_grad,
    # #     beta_ind,
    # #     beta_groups
    # # )
    # fp.insert_or_update_data({'fine_faces_by_grad': fine_faces_by_grad})

    # if fine_faces_by_grad.shape[0] > 0:
    #     fine_levels[fine_faces_by_grad] = 0
    # finescale_ids = fp['faces'][fine_levels==0]

    # LEVEL_ID_1, ADM_COARSE_ID_LEVEL_1 = nu_adm_funcs.set_adm_mesh_non_nested(
    #     finescale_ids,
    #     fine_levels,
    #     fp['faces'],
    #     fp[defnames.get_primal_id_name_by_level(1)],
    #     fp[defnames.get_dual_id_name_by_level(1)]
    # )

    # OP_adm, OR_adm = nu_adm_funcs.organize(
    #     fine_levels,
    #     sp.find(OP)[0:3],
    #     fp['faces'],
    #     fp[defnames.get_primal_id_name_by_level(1)],
    #     LEVEL_ID_1,
    #     ADM_COARSE_ID_LEVEL_1,
    #     fp[defnames.get_dual_id_name_by_level(1)]
    # )

    # T_adm = OR_adm*(resp['transmissibility']*OP_adm)
    # Q_adm = OR_adm*resp['source']
    # P_adm = spsolve(T_adm.tocsc(), Q_adm)
    # P_prol = OP_adm*P_adm

    edges_flux, nodes_pressure = lsds.get_edges_flux_and_nodes_pressure(
        bc,
        P_prol,
        fp['xi_params'],
        fp['nodes_weights'],
        fp['nodes_of_edges'],
        fp['adjacencies'],
        fp['neumann_weights']
    )

    
    
    update_fine_flux(
        coarse_struct,
        edges_flux,
        P_prol,
        total_mobility_edges,
        lsds,
        nodes_pressure,
        bc,
        fp['nodes_of_edges'],
        fp.edges_dim,
        finescale_ids,
        saturation
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
    
    fp.insert_or_update_data(
        {
            'nuadm_vols': np.array([OP_adm.shape[1]]),
            'it': np.array([it]),
            'err': np.array([err])
        }
    )

    return P_prol, newS, new_vpi, new_cumulative_oil, new_cumulative_water, faces_flux, fine_levels, water_faces_flux, dt, water_flux, oil_flux, plot_vpi

def update_data(
        simulation_data: SimulationData,
        vpi: float,
        cumulative_oil: float,
        cumulative_water: float,
        loop: int,
        pressure: np.ndarray,
        saturation: np.ndarray,
        water_flux: float,
        oil_flux: float,
        fp: MeshProperty
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

    all_nu_adm_vols = simulation_data['nuadm_vols']
    all_nu_adm_vols = np.append(all_nu_adm_vols, fp['nuadm_vols'])
    
    all_it = simulation_data['all_it']
    all_it = np.append(all_it, fp['it'])
    
    all_err = simulation_data['all_err']
    all_err = np.append(all_err, fp['err'])
    
    simulation_data.insert_or_update_data({
        simulation_data.my_data_names[0]: all_loops,
        simulation_data.my_data_names[1]: all_vpi,
        simulation_data.my_data_names[2]: all_cum_oil,
        simulation_data.my_data_names[3]: all_cum_wat,
        # simulation_data.my_data_names[4] + str(loop): pressure,
        # simulation_data.my_data_names[5] + str(loop): saturation,
        simulation_data.my_data_names[6]: all_water_flux,
        simulation_data.my_data_names[7]: all_oil_flux,
        'nuadm_vols': all_nu_adm_vols,
        'all_it': all_it,
        'all_err': all_err
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
    
def define_initial_fine_volumes(fp: MeshProperty, bc: BoundaryConditions):
    dirichlet_vols = bc['dirichlet_volumes']['id']
    values = bc['dirichlet_volumes']['id']

    cids = fp[defnames.get_primal_id_name_by_level(1)][dirichlet_vols]
    test = np.isin(fp[defnames.get_primal_id_name_by_level(1)], cids)
    fine_vols = fp['faces'][test]
    return fine_vols

def refine_from_permeability_contrast(fp: MeshProperty):
    perm = fp['permeabiity']
    permx = perm[:, 0, 0]

    adjacencies = fp['adjacencies']

    dperm = permx[adjacencies[fp.internal_edges]]

    ddperm = np.absolute(dperm[:, 0] - dperm[:, 1])
    test = ddperm > 0
    edges1 = fp.internal_edges[test]
    volumes = np.unique(adjacencies[edges1].flatten())
    return volumes

def refine_from_permeability_value_v0(fp: MeshProperty):
    perm = fp['permeability']
    permx = perm[:, 0, 0]

    adjacencies = fp['adjacencies']

    test = permx <= 1e-5

    vols1 = fp.faces[test]

    t1 = np.isin(adjacencies[:, 0], vols1)
    t2 = np.isin(adjacencies[:, 1], vols1)
    t3 = t1 | t2

    vols2 = np.unique(adjacencies[t3].flatten())

    return vols2

def refine_from_permeability_value_v2(fp: MeshProperty):
    perm = fp['permeability']
    primal_id = fp[defnames.get_primal_id_name_by_level(1)]

    permx = perm[:, 0, 0]

    test = permx <= 1e-5

    vols1 = fp.faces[test]

    primal_ids = np.unique(primal_id[vols1])

    test = np.isin(primal_id, primal_ids)

    vols2 = fp.faces[test]


    return vols2

def refine_from_permeability_value(fp: MeshProperty):
    return np.array([])


def update_pressure_only_ms(
        relative_perm: BrooksAndCorey,
        biphasic_mobility: BiphasicMobility,
        saturation: np.ndarray,
        fp: MeshProperty,
        bc: BoundaryConditions,
        lsds: LsdsFluxCalculation,
        initial_fine_volumes: np.ndarray,
        alpha_lim_finescale: float,
        beta_lim: float,
        OP: sp.csc_matrix,
        OR: sp.csc_matrix,
        coarse_struct: Sequence[PrimalCoarseData],
        iterative_ms: bool,
        tol_iterative: float,
        p0: np.ndarray,
        **kwargs
):
    
    it = -1
    err = -1
    
    krw_faces, kro_faces = relative_perm.calculate(saturation)
    mobw_faces, mobo_faces = biphasic_mobility.calculate(krw_faces, kro_faces)
    total_mobility_faces = biphasic_mobility.get_total_mobility(mobw_faces, mobo_faces)
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
    t0 = time.perf_counter()
    update_weight_new_function(fp, saturation)
    t1 = time.perf_counter()
    TimeProfile.dt_weight_finescale = t1 - t0

    t0 = time.perf_counter()
    for coarse_data in coarse_struct:
        coarse_data.get_local_nodes_weights(
            fp['nodes_weights'],
            coarse_data['map_nodes'],
            coarse_data['nodes'],
            coarse_data['map_faces'],
            coarse_data['faces'],
            coarse_data['bool_boundary_nodes'],
            fp['nodes_to_calculate']
        )
    t1 = time.perf_counter()
    TimeProfile.dt_update_coarse_weight = t1 - t0

    t0 = time.perf_counter()
    resp = set_fine_transmissibility_biphasic(
        fp,
        bc,
        lsds
    )
    t1 = time.perf_counter()
    TimeProfile.dt_set_finescale_problem = t1 - t0

    t0 = time.perf_counter()
    fine_levels = np.full(fp['faces'].shape[0], -1)
    fine_levels[initial_fine_volumes] = 0
    fine_levels[fp['initial_fine_faces']] = 0
    # fine_ids_from_saturation = define_fine_ids_from_saturation(saturation, fp['adjacencies'], fp.internal_edges)
    fine_ids_from_saturation = define_fine_ids_from_saturation_internal_nodes_org(fp['internal_nodes_org'], fp['internal_faces_of_nodes_org'], saturation)

    fine_levels[fine_ids_from_saturation] = 0

    # fine_ids_from_perm = refine_from_permeability_value(fp)
    # if fine_ids_from_perm.shape[0] > 0:
    #     fine_levels[fine_ids_from_perm] = 0

    fine_ids_from_alpha = fine_level_from_alpha.define_fine_levels_from_alpha(
        OR,
        OP,
        resp['transmissibility'],
        fp[defnames.get_primal_id_name_by_level(1)],
        alpha_lim=alpha_lim_finescale
    )

    fine_levels[fine_ids_from_alpha] = 0

    beta_groups, beta_ind, betas = nu_adm_funcs.get_beta_groups(
        fp['faces'],
        fp[defnames.get_primal_id_name_by_level(1)],
        sp.find(OP)[0:3],
        fp['adjacencies'][fp.internal_edges],
        beta_lim=beta_lim
    )

    finescale_faces = nu_adm_funcs.get_finescale_vols(
        fp['faces'][fine_levels==0],
        fine_ids_from_alpha,
        beta_ind,
        beta_groups
    )

    finescale_faces = finescale_faces.astype(np.int)

    fine_levels[finescale_faces] = 0
    fine_levels[fine_levels==-1] = 1

    finescale_ids = fp['faces'][fine_levels==0]

    LEVEL_ID_1, ADM_COARSE_ID_LEVEL_1 = nu_adm_funcs.set_adm_mesh_non_nested(
        finescale_ids,
        fine_levels,
        fp['faces'],
        fp[defnames.get_primal_id_name_by_level(1)],
        fp[defnames.get_dual_id_name_by_level(1)]
    )

    OP_adm, OR_adm = nu_adm_funcs.organize(
        fine_levels,
        sp.find(OP)[0:3],
        fp['faces'],
        fp[defnames.get_primal_id_name_by_level(1)],
        LEVEL_ID_1,
        ADM_COARSE_ID_LEVEL_1,
        fp[defnames.get_dual_id_name_by_level(1)]
    )
    t1 = time.perf_counter()
    TimeProfile.dt_set_adm_mesh = t1 - t0
    
    t0 = time.perf_counter()
    if iterative_ms:
        P_prol, it, err = iterative_ms_ilu0_bicgstab(
            resp['transmissibility'],
            resp['source'],
            p0,
            # P_prol,
            OP_adm,
            OR_adm,
            epsilon=tol_iterative
        )
    
    else:
        T_adm = OR_adm*(resp['transmissibility']*OP_adm)
        Q_adm = OR_adm*resp['source']
        P_adm = spsolve(T_adm.tocsc(), Q_adm)
        P_prol = OP_adm*P_adm
    t1 = time.perf_counter()
    TimeProfile.dt_solution_ms = t1 - t0

    edges_flux, nodes_pressure = lsds.get_edges_flux_and_nodes_pressure(
        bc,
        P_prol,
        fp['xi_params'],
        fp['nodes_weights'],
        fp['nodes_of_edges'],
        fp['adjacencies'],
        fp['neumann_weights']
    )
    
    fp.insert_or_update_data({'dt_update_boundary_nodes_weights_neumann': np.array([0])})
    
    t0 = time.perf_counter()
    update_fine_flux(
        coarse_struct,
        edges_flux,
        P_prol,
        total_mobility_edges,
        lsds,
        nodes_pressure,
        bc,
        fp['nodes_of_edges'],
        fp.edges_dim,
        finescale_ids,
        saturation
    )
    t1 = time.perf_counter()
    TimeProfile.dt_neumann = t1-t0
    
    fp.insert_or_update_data({
        'fine_levels': fine_levels
    })
    
    return P_prol, edges_flux


def run():

    matrices_path = 'matrices.h5'
    op_name = 'AMS-U'

    update_primal_mesh = True
    update_dual_mesh = True
    update_coarse_struct = True
    my_dual_type = 1
    cfl = 0.9

    alpha_lim_finescale = 0.05
    beta_lim = 3.0
    

    type_k = 'barrier'
    dt = 0.00005
    max_vpi = 1.3
    loop = 0
    max_loop = np.inf
    load = False
    loop_intervals = 1
    etol_msrsb = 0.01
    maxit_msrsb = 1000

    cumulative_oil = 0.0
    cumulative_water = 0.0
    vpi = 0.0

    relative_perm = BrooksAndCorey()
    biphasic_mobility = BiphasicMobility()
    lsds = LsdsFluxCalculation()
    simulation_data = SimulationData('biphasic_nu_adm')

    fp, cp, fine_mesh_path, coarse_mesh_path = get_properties()
    bc = set_boundary_conditions(fp)
    create_primal_ids(fp, cp, update=update_primal_mesh)
    export_primal_ids(fine_mesh_path, fp, coarse_mesh_path, export=update_primal_mesh)
    create_dual_ids(fp, cp, update=update_dual_mesh, dual_type=my_dual_type)
    export_dual_ids(fine_mesh_path, fp, export=update_dual_mesh)
    
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
    # # mesh_data.export_all_elements_type_to_vtk('pressure_faces_' + str(loop), 'faces')

    # initial_fine_vols = define_initial_fine_volumes(fp, bc)
    initial_fine_vols = define_new_fine_levels_v1(fp, bc)
   

    if load is False:
        initial_funcs(fp, fine_mesh_path, type_k)
        # coarse_struct = define_coarse_structure(fp, lsds, level=1, update=update_coarse_struct)
        pressure[:], newS[:], vpi, cumulative_oil, cumulative_water, faces_flux, coarse_struct, OP, OR, fine_levels, water_flux, oil_flux = initial_loop(
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
            matrices_path,
            op_name,
            initial_fine_vols,
            alpha_lim_finescale,
            beta_lim,
            etol_msrsb,
            maxit_msrsb
        )
        mesh_data.insert_tag_data('pressure', pressure, 'faces')
        mesh_data.insert_tag_data('saturation', saturation, 'faces')
        mesh_data.insert_tag_data('faces_flux', np.absolute(faces_flux), 'faces')
        mesh_data.export_all_elements_type_to_vtk('pressure_faces_' + str(loop), 'faces')
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
        saturation_plot[:] = saturation
        saturation[:] = newS
        adm_interfaces_name = 'adm_edges_' + str(loop)
        print_adm_interfaces_2d(
            fp,
            fine_mesh_path,
            fine_levels,
            adm_interfaces_name
        )
        fp.export_data()
    else:
        # import pdb; pdb.set_trace()
        simulation_data.load_data()
        loop = simulation_data['all_loops'][-1]
        vpi = simulation_data['all_vpi'][-1]
        cumulative_oil = simulation_data['all_cumulative_oil'][-1]
        cumulative_water = simulation_data['all_cumulative_water'][-1]
        saturation[:] = simulation_data['saturation_' + str(loop)]
        pressure[:] = simulation_data['pressure_' + str(loop)]
        coarse_struct = define_coarse_structure(fp, lsds, level=1, update=False)
        OP = utils_old.load_matrix(matrices_path, op_name)
        OR = get_OR_AMS(fp)

    while vpi < max_vpi and loop < max_loop:
        for i in range(loop_intervals):
            loop += 1
            pressure[:], newS[:], vpi, cumulative_oil, cumulative_water, faces_flux, fine_levels, water_faces_flux, dt, water_flux, oil_flux, plot_vpi = while_loop(
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
                matrices_path,
                op_name,
                initial_fine_vols,
                alpha_lim_finescale,
                beta_lim,
                OP,
                OR,
                coarse_struct,
                cfl
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
            
            if plot_vpi == True:
                break
        
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

        adm_interfaces_name = 'adm_edges_' + str(loop)
        print_adm_interfaces_2d(
            fp,
            fine_mesh_path,
            fine_levels,
            adm_interfaces_name
        )

        mesh_data.insert_tag_data('pressure', pressure, 'faces')
        mesh_data.insert_tag_data('faces_flux', faces_flux, 'faces')
        mesh_data.insert_tag_data('water_faces_flux', water_faces_flux, 'faces')
        mesh_data.insert_tag_data('saturation', saturation_plot, 'faces')
        mesh_data.export_all_elements_type_to_vtk('pressure_faces_' + str(loop), 'faces')
        
        

        # mesh_data.insert_tag_data('pressure', pressure, 'faces')
        # mesh_data.insert_tag_data('faces_flux', faces_flux, 'faces')
        # mesh_data.insert_tag_data('water_faces_flux', water_faces_flux, 'faces')
        # mesh_data.insert_tag_data('saturation', newS, 'faces')
        # mesh_data.export_all_elements_type_to_vtk('pressure_faces_' + str(loop), 'faces')


        # import pdb; pdb.set_trace()

    import pdb; pdb.set_trace()



    

    

    




    


    




    















    print('fim')
