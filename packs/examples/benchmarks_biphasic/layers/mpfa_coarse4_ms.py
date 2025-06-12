from packs.biphasic.relative_perm.brooks_and_corey import BrooksAndCorey
from packs.biphasic.mobility import BiphasicMobility
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths, defnames
from packs.manager import MeshProperty, MeshData, BoundaryConditions, SimulationData
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.examples.diss_test1 import define_new_fine_levels_v1

from packs.examples.same_functions import (
    define_faces_in_losangle, 
    define_coarse_structure,
    update_xi_params,
    set_fine_transmissibility_biphasic,
    update_fine_flux,
    define_fine_ids_from_saturation
)

from packs.examples.benchmarks_biphasic.layers.mpfa_fine4 import (
    get_properties as get_properties_finescale,
    set_boundary_conditions,
    initial_funcs
)

from packs.examples.biphasic_mpfa_nu_adm import (
    initial_loop,
    while_loop,
    update_data
)

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


import os
import numpy as np
from typing import Tuple, Sequence
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from typing import Sequence

def get_properties_coarse():
    coarse_mesh_properties_name = 'coarse4_layers'
    coarse_mesh_path = defpaths.coarse4_layers

    coarse_properties = preprocess_mesh(coarse_mesh_path, coarse_mesh_properties_name)

    return coarse_properties, coarse_mesh_path

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
        matrices_path: str,
        op_name: str,
        initial_fine_vols: np.ndarray,
        alpha_lim_finescale: float,
        beta_lim: float,
        saturation_plot: np.ndarray,
        etol_msrsb,
        maxit_msrsb,
        cfl
):
    
    if load is False:
        initial_funcs(fp, fine_mesh_path)
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
            maxit_msrsb,
            cfl
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
            'oil_flux': np.array([oil_flux])
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
    
    return loop, cumulative_oil, cumulative_water, vpi, OP, OR, coarse_struct

def update_while_loop_ms(
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
        matrices_path,
        op_name,
        initial_fine_vols,
        alpha_lim_finescale,
        beta_lim,
        OP,
        OR,
        coarse_struct,
        fine_mesh_path
):
    
    for i in range(loop_intervals):
        loop += 1
        pressure[:], newS[:], vpi, cumulative_oil, cumulative_water, faces_flux, fine_levels, water_faces_flux, dt, water_flux, oil_flux = while_loop(
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

    return loop, cumulative_oil, cumulative_water, vpi


def run6():
    matrices_path = 'matrices.h5'
    op_name = 'AMS-U'
    debug = False

    update_primal_mesh = True
    update_dual_mesh = True
    update_coarse_struct = True

    # update_primal_mesh = False
    # update_dual_mesh = False
    # update_coarse_struct = False


    my_dual_type = 1
    cfl = 0.9

    # alpha_lim_finescale = 0.1
    # beta_lim = 2

    alpha_lim_finescale = 1e3
    beta_lim = 1e3

    dt = 0.00005
    max_vpi = 1.3
    loop = 0
    max_loop = np.inf
    load = False
    loop_intervals = 10
    etol_msrsb = 0.01
    maxit_msrsb = 1000

    cumulative_oil = 0.0
    cumulative_water = 0.0
    vpi = 0.0

    relative_perm = BrooksAndCorey(Sor=0.0, Swc=0.0)
    biphasic_mobility = BiphasicMobility(mio=4)
    lsds = LsdsFluxCalculation()
    simulation_data = SimulationData('biphasic_layers_coarse4')
    
    fp, fine_mesh_path = get_properties_finescale()
    cp, coarse_mesh_path = get_properties_coarse()
    nodes_org = fp.get_internal_nodes_org_from_faces_of_nodes_object()
    fp.insert_or_update_data(nodes_org)
    bc = set_boundary_conditions(fp)
    create_primal_ids(fp, cp, update=update_primal_mesh)
    export_primal_ids(fine_mesh_path, fp, coarse_mesh_path, export=update_primal_mesh)
    create_dual_ids(fp, cp, update=update_dual_mesh, dual_type=my_dual_type)
    export_dual_ids(fine_mesh_path, fp, export=update_dual_mesh)

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

    initial_fine_vols = define_new_fine_levels_v1(fp, bc)

    loop, cumulative_oil, cumulative_water, vpi, OP, OR, coarse_struct = load_or_update_initial_loop(
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
        matrices_path,
        op_name,
        initial_fine_vols,
        alpha_lim_finescale,
        beta_lim,
        saturation_plot,
        etol_msrsb,
        maxit_msrsb,
        cfl
    )

    import pdb; pdb.set_trace()

    while vpi < max_vpi and loop < max_loop:

        loop, cumulative_oil, cumulative_water, vpi = update_while_loop_ms(
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
            matrices_path,
            op_name,
            initial_fine_vols,
            alpha_lim_finescale,
            beta_lim,
            OP,
            OR,
            coarse_struct,
            fine_mesh_path
        )

    



    import pdb; pdb.set_trace()



