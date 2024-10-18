from packs.biphasic.relative_perm.brooks_and_corey import BrooksAndCorey
from packs.biphasic.mobility import BiphasicMobility
from packs.biphasic.unstructured.mobility_mesh_elements import direct_edges_mobility
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths, defnames
from packs.manager import MeshProperty, MeshData, BoundaryConditions, SimulationData
from packs.multiscale.unstructured.test.test_cross import set_weights_nodes, set_fine_transmissibility
from packs.multiscale.unstructured.test.test_brazil import define_faces_in_losangle, set_permeability
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights

from packs.multiscale.unstructured.test.test_brazil import create_primal_ids, create_dual_ids, export_primal_ids, export_dual_ids
from packs.multiscale.unstructured.create_primal_dual.primal_coarse_volumes_2d import get_coarse_structure, load_coarse_structure
from packs.multiscale.unstructured.test.test_brazil import get_OR_AMS
from packs.multiscale.unstructured.operators.prolongation.dual_interaction_region import create_dual_interaction_regions
from packs.multiscale.unstructured.operators.prolongation.get_op_from_amsu import update_global_op_from_amsu
from packs.utils import utils_old
from packs.adm.non_uniform import fine_level_from_alpha
from packs.fim_nu_adm.packs.processor import nu_adm_funcs
from packs.manager.generic_data import PrimalCoarseData

import os
import numpy as np
from typing import Tuple, Sequence
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

def get_properties():
    rel_path = os.path.join(
        defpaths.unstructured_coarse_test_mesh_folder,
        'brazil'
    )
    fine_mesh_path = os.path.join(rel_path, 'brazilf.msh')
    fine_mesh_properties_name = 'brazilf' 
    fine_mesh_path_v4 = os.path.join(rel_path, 'brazilf_v4.msh')

    coarse_mesh_path = os.path.join(rel_path, 'brazilC.msh')
    coarse_mesh_properties_name = 'brazilC1'

    # coarse_mesh_path = os.path.join(rel_path, 'brazilC2.msh')
    # coarse_mesh_properties_name = 'brazilC2'


    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)
    coarse_properties = preprocess_mesh(coarse_mesh_path, coarse_mesh_properties_name)

    return fine_properties, coarse_properties, fine_mesh_path, coarse_mesh_path

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
    faces_pressure = np.array([face_p1, face_p0])
    pressure_presc = np.array([1.0, 0.0])

    bc.set_boundary('dirichlet_volumes', faces_pressure, pressure_presc)
    bc.set_boundary('dirichlet_nodes', np.array([]), np.array([]))

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

def calculate_dt(faces_centroids: np.ndarray, adjacencies: np.ndarray, bool_boundary_edges: np.ndarray, total_flux_edges: np.ndarray, edges_dim: np.ndarray, fw_faces: np.ndarray, saturation: np.ndarray, porosity: np.ndarray, cfl: float=1.0):
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
    
def update_saturation(water_faces_flux, areas, dt, porosity, saturation):

    ds = dt*water_faces_flux/(porosity*areas)
    newS = saturation + ds
    return newS

def update_simulation_data(faces_flux: np.ndarray, injectors: np.ndarray, producers: np.ndarray, vpi: list, cum_oil: list, cum_water: list, fw_faces: np.ndarray, total_area_reservoir: float, dt: float) -> None:
    total_volume_injected = faces_flux[injectors].sum()*dt
    water_flux = (faces_flux[producers]*fw_faces[producers]).sum()
    total_volume_water_produced = water_flux*dt
    fo_faces = 1-fw_faces
    oil_flux = (faces_flux[producers]*fo_faces[producers]).sum()
    total_volume_oil_produced = oil_flux*dt

    vpi += total_volume_injected/total_area_reservoir
    cum_oil += total_volume_oil_produced
    cum_water += total_volume_water_produced

    return vpi, cum_oil, cum_water

def set_fine_transmissibility_biphasic(fine_mesh_properties: MeshProperty, bc: BoundaryConditions, lsds: LsdsFluxCalculation):
    resp = lsds.mount_transmissibility_matrix(
        bc,
        **fine_mesh_properties.get_all_data()
    )
    return resp

def set_weights_nodes_cstruct(cstruct: PrimalCoarseData):

    fine_properties = cstruct
    weights = get_gls_nodes_weights(**fine_properties.get_all_data())
    return weights

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

def initial_funcs(
        fp: MeshProperty,
        fine_mesh_path: str,
        type_k: str
):
    define_faces_in_losangle(fp)
    set_permeability(fine_mesh_path, fp, typek=type_k, export_permfield=True, update_permfield=True)
    set_weights_nodes(fp, update=True)
    fp.backup_data('xi_params', 'xi_params_backup')
    fp.export_data()

def update_fine_flux(
        coarse_struct: Sequence[PrimalCoarseData],
        edges_flux: np.ndarray,
        ms_pressure: np.ndarray,
        total_mobility_edges: np.ndarray,
        lsds: LsdsFluxCalculation,
        nodes_pressure: np.ndarray
    ):

    for cstruct in coarse_struct:
        global_edges = cstruct['map_edges']
        bool_boundary_edges = cstruct['bool_boundary_edges']
        local_flux_presc = edges_flux[global_edges[bool_boundary_edges]]
        neumann_edges = cstruct['edges'][bool_boundary_edges]
        bc = BoundaryConditions()
        # bc.set_boundary('neumann_edges', neumann_edges, local_flux_presc)
        bc.set_boundary('neumann_edges', np.array([]), np.array([]))

        # local_vertice = cstruct['faces'][cstruct['dual_id']==defnames.dual_ids('vertice_id')]
        # bc.set_boundary('dirichlet_volumes', local_vertice, ms_pressure[local_vertice])
        bc.set_boundary('dirichlet_volumes', np.array([]), np.array([]))
        # bc.set_boundary('dirichlet_nodes', np.array([]), np.array([]))

        bool_boundary_nodes = cstruct['bool_boundary_nodes']
        mapbnodes = cstruct['map_nodes'][bool_boundary_nodes]
        bnodes = cstruct['nodes'][bool_boundary_nodes]
        bc.set_boundary('dirichlet_nodes', bnodes, nodes_pressure[mapbnodes])

        edges_multiplier = total_mobility_edges[global_edges].copy()
        ########################
        # edges_multiplier[:] = 1
        #########################
        cstruct.update_data({
            'neumann_edges': neumann_edges,
            'neumann_edges_value': local_flux_presc,
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

        lt = set_fine_transmissibility_biphasic(cstruct, bc, lsds)
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

        # if cstruct['coarse_id'][0] == 1:
        #     import pdb; pdb.set_trace()

        cstruct.insert_or_update_data({
            'edges_flux': local_edges_flux
        })

        bool_internal_edges = ~bool_boundary_edges
        edges_flux[cstruct['map_edges'][bool_internal_edges]] = local_edges_flux[bool_internal_edges]











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
        coarse_struct
):
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
        'xi_params': update_xi_params(fp['xi_params_backup'], total_mobility_edges)
    })

    resp = set_fine_transmissibility_biphasic(
        fp,
        bc,
        lsds
    )

    fine_transm_without_bc = lsds.mount_transmissibility_matrix_without_bc(**fp.get_all_data())

    OP = get_OP_matrix(fp)
    OR = get_OR_AMS(fp)
    OP = get_op_amsu(fp, lsds, OP, fine_transm_without_bc['transmissibility_without_bc'])
    utils_old.save_matrix(matrices_path, op_name, OP)

    fine_levels = np.full(fp['faces'].shape[0], -1)
    fine_levels[initial_fine_volumes] = 0

    fine_ids_from_alpha = fine_level_from_alpha.define_fine_levels_from_alpha(
        OR,
        OP,
        fine_transm_without_bc['transmissibility_without_bc'],
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
        nodes_pressure
    )

    # pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])

    # edges_flux2 = lsds.get_edges_flux(
    #     bc,
    #     pressure,
    #     fp['xi_params'],
    #     fp['nodes_weights'],
    #     fp['nodes_of_edges'],
    #     fp['adjacencies'],
    #     fp['neumann_weights']
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

    newS = update_saturation(water_faces_flux, fp['areas'], dt, porosity, saturation)
    relative_perm._test_saturations(newS)

    new_vpi, new_cumulative_oil, new_cumulative_water = update_simulation_data(
        faces_flux,
        bc['injectors']['id'],
        bc['producers']['id'],
        vpi,
        cumulative_oil,
        cumulative_water,
        fw_faces,
        total_area_reservoir,
        dt
    )

    return P_prol, newS, new_vpi, new_cumulative_oil, new_cumulative_water, faces_flux

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
        total_area_reservoir: float
):
    
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
        'xi_params': update_xi_params(fp['xi_params_backup'], total_mobility_edges)
    })

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
        cfl=0.5
    )

    newS = update_saturation(water_faces_flux, fp['areas'], dt, porosity, saturation)
    new_vpi, new_cumulative_oil, new_cumulative_water = update_simulation_data(
        faces_flux,
        bc['injectors']['id'],
        bc['producers']['id'],
        vpi,
        cumulative_oil,
        cumulative_water,
        fw_faces,
        total_area_reservoir,
        dt
    )

    return pressure, newS, new_vpi, new_cumulative_oil, new_cumulative_water

def update_data(
        simulation_data: SimulationData,
        vpi: float,
        cumulative_oil: float,
        cumulative_water: float,
        loop: int,
        pressure: np.ndarray,
        saturation: np.ndarray
):
    
    all_loops = simulation_data['all_loops']
    all_loops = np.append(all_loops, [loop])

    all_vpi = simulation_data['all_vpi']
    all_vpi = np.append(all_vpi, [vpi])

    all_cum_oil = simulation_data[simulation_data.my_data_names[2]]
    all_cum_oil = np.append(all_cum_oil, [cumulative_oil])

    all_cum_wat = simulation_data[simulation_data.my_data_names[3]]
    all_cum_wat = np.append(all_cum_wat, [cumulative_water])
    
    simulation_data.insert_or_update_data({
        simulation_data.my_data_names[0]: all_loops,
        simulation_data.my_data_names[1]: all_vpi,
        simulation_data.my_data_names[2]: all_cum_oil,
        simulation_data.my_data_names[3]: all_cum_wat,
        simulation_data.my_data_names[4] + str(loop): pressure,
        simulation_data.my_data_names[5] + str(loop): saturation
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

def define_coarse_structure(fine_mesh_properties: MeshProperty, lsds: LsdsFluxCalculation, level=1, update=True):
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
        
def define_initial_fine_volumes(fp: MeshProperty, bc: BoundaryConditions):
    dirichlet_vols = bc['dirichlet_volumes']['id']
    values = bc['dirichlet_volumes']['id']

    cids = fp[defnames.get_primal_id_name_by_level(1)][dirichlet_vols]
    test = np.isin(fp[defnames.get_primal_id_name_by_level(1)], cids)
    fine_vols = fp['faces'][test]
    return fine_vols




def run():

    matrices_path = 'matrices.h5'
    op_name = 'AMSU'

    update_primal_mesh = False
    update_coarse_struct = False
    update_dual_mesh = False
    my_dual_type = 1

    alpha_lim_finescale = 0.5
    beta_lim = 3.0
    

    type_k = 'barrier'
    dt = 0.00005
    max_vpi = 1.3
    loop = 0
    max_loop = np.inf
    load = False
    loop_intervals = 10

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

    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.create_tag('pressure')
    mesh_data.create_tag('faces_flux')
    # mesh_data.create_tag('water_faces_flux')
    mesh_data.create_tag('saturation')
    # # mesh_data.export_all_elements_type_to_vtk('pressure_faces_' + str(loop), 'faces')

    initial_fine_vols = define_initial_fine_volumes(fp, bc)
   

    if load is False:
        initial_funcs(fp, fine_mesh_path, type_k)
        coarse_struct = define_coarse_structure(fp, lsds, level=1, update=update_coarse_struct)
        pressure[:], newS[:], vpi, cumulative_oil, cumulative_water, faces_flux = initial_loop(
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
            coarse_struct
        )
        mesh_data.insert_tag_data('pressure', pressure, 'faces')
        mesh_data.insert_tag_data('saturation', saturation, 'faces')
        mesh_data.insert_tag_data('faces_flux', np.absolute(faces_flux), 'faces')
        mesh_data.export_all_elements_type_to_vtk('pressure_faces_' + str(loop), 'faces')
        import pdb; pdb.set_trace()
        simulation_data.insert_or_update_data({
            'all_loops': np.array([0]),
            'all_vpi': np.array([0.0]),
            'all_cumulative_oil': np.array([0.0]),
            'all_cumulative_water': np.array([0.0]),
            'pressure_' + str(loop): pressure,
            'saturation_' + str(loop): saturation
        })
        saturation[:] = newS
    else:
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
            pressure[:], newS[:], vpi, cumulative_oil, cumulative_water = while_loop(
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
                total_area_reservoir
            )
            saturation[:] = newS
            
            print()
            print('##########################')
            print(f'VPI: {vpi}')
            print(f'Cum oil: {cumulative_oil}')
            print(f'Cum water: {cumulative_water}')
            print(f'Loop: {loop}')
            print('##########################')
            print()
        
        import pdb; pdb.set_trace()
        
        update_data(
            simulation_data,
            vpi,
            cumulative_oil,
            cumulative_water,
            loop,
            pressure,
            saturation
        )
        
        

        # mesh_data.insert_tag_data('pressure', pressure, 'faces')
        # mesh_data.insert_tag_data('faces_flux', faces_flux, 'faces')
        # mesh_data.insert_tag_data('water_faces_flux', water_faces_flux, 'faces')
        # mesh_data.insert_tag_data('saturation', newS, 'faces')
        # mesh_data.export_all_elements_type_to_vtk('pressure_faces_' + str(loop), 'faces')


        # import pdb; pdb.set_trace()

    import pdb; pdb.set_trace()



    

    

    




    


    




    















    print('fim')
