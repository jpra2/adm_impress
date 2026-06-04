from packs.biphasic.relative_perm.brooks_and_corey import BrooksAndCorey
from packs.biphasic.mobility import BiphasicMobility
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths, defnames
from packs.manager import MeshProperty, MeshData, BoundaryConditions, SimulationData, configsim
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.examples.diss_test1 import define_new_fine_levels_v1

from packs.examples.same_functions import (
    define_coarse_structure,
    create_folders_pressure_results,
    export_ps_results,
    load_ps_results
)

from packs.examples.benchmarks_biphasic.coarse_test1.ftest1 import (
    get_properties as get_properties_finescale,
    set_boundary_conditions,
    initial_funcs,
    create_path_mesh_data,
    set_permeability
)

from packs.examples.benchmarks_biphasic.ameba.mpfa_coarse4 import (
    plot_results_ms
)

from packs.examples.biphasic_mpfa import (
    update_saturation_only
)

from packs.examples.benchmarks_biphasic.ameba.mpfa_fine4_mimpes import (
    set_dvtol,
    calculate_time_step_mimpes
)

from packs.examples.benchmarks_biphasic.ameba.mpfa_coarse4_mimpes import (
    load_or_update_initial_loop,
    update_while_loop_ms
)

from packs.examples.biphasic_mpfa_nu_adm import (
    initial_loop,
    while_loop,
    update_data,
    update_pressure_only_ms
)


from packs.utils.multiscale_methods import print_adm_interfaces_2d

from packs.multiscale.unstructured.test.test_brazil import create_primal_ids, create_dual_ids, export_primal_ids, export_dual_ids
from packs.multiscale.unstructured.test.test_brazil import get_OR_AMS
from packs.multiscale.unstructured.operators.prolongation.dual_interaction_region import create_dual_interaction_regions
from packs.multiscale.unstructured.operators.prolongation.get_op_from_amsu import update_global_op_from_amsu
from packs.utils import utils_old
from packs.adm.non_uniform import fine_level_from_alpha
from packs.fim_nu_adm.packs.processor import nu_adm_funcs
from packs.manager.generic_data import PrimalCoarseData
from packs.manager.predef_names import TimeProfile


import os
import numpy as np
from typing import Tuple, Sequence
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from typing import Sequence
import time

def get_properties_coarse():
    coarse_mesh_properties_name = 'ctest1_2'
    coarse_mesh_path = defpaths.ctest1_2

    # coarse_mesh_properties_name = 'coarse4f_ameba'
    # coarse_mesh_path = defpaths.ameba_coarse4f

    coarse_properties = preprocess_mesh(coarse_mesh_path, coarse_mesh_properties_name)

    return coarse_properties, coarse_mesh_path

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

    alpha_lim_finescale = 1e6
    beta_lim = 1e6

    dt = 0.00005
    max_vpi = 0.6
    loop = 0
    max_loop = np.inf
    load = False
    loop_intervals = 1
    etol_msrsb = 0.01
    maxit_msrsb = 1000
    vpis_to_plot = np.linspace(0, 0.6, 301)[1:]
    # dvtol = 1e-5
    Rdtmax = 1.3
    Rdtmin = 0.75
    saturation_intervals = 10
    iterative_ms = False
    tol_iterative = 1e-6

    refine_by_grad_bool = False
    refine_by_estimator1_bool = True
    max_value_grad = 100
    max_value_estimator1 =0.04
    
    gdict = {
        'funcname': '',
        'file_times': 'functions_times_ctest1_2.yaml',
        'load_simulation': load,
        'filename_export_times': os.path.join(defpaths.flying, 'pressure_times.csv')
    }

    cumulative_oil = 0.0
    cumulative_water = 0.0
    vpi = 0.0

    relative_perm = BrooksAndCorey(Sor=0.0, Swc=0.0)
    biphasic_mobility = BiphasicMobility(mio=4)
    lsds = LsdsFluxCalculation()
    # simulation_data = SimulationData('biphasic_layers_coarse4')
    simulation_data = SimulationData('bipdf = pd.DataFrame([data])hasic_ctest1_2')
    simulation_data.insert_or_update_data({'label': np.array(['ctest1_2'])})
    create_path_mesh_data(simulation_data)
    TimeProfile.config_sim = simulation_data['label'][0]

    
    fp, fine_mesh_path = get_properties_finescale()
    set_permeability(fine_mesh_path, fp, simulation_data)
    cp, coarse_mesh_path = get_properties_coarse()
    nodes_org = fp.get_internal_nodes_org_from_faces_of_nodes_object()
    fp.insert_or_update_data(nodes_org)
    bc = set_boundary_conditions(fp)
    gdict.update({'funcname': 'create_primal_ids'})
    create_primal_ids(fp, cp, update=update_primal_mesh,**gdict)
    export_primal_ids(fine_mesh_path, fp, coarse_mesh_path, export=update_primal_mesh)
    gdict.update({'funcname': 'create_dual_ids'})
    create_dual_ids(fp, cp, update=update_dual_mesh, dual_type=my_dual_type,**gdict)
    export_dual_ids(fine_mesh_path, fp, export=update_dual_mesh)

    porosity = np.repeat(0.2, len(fp['faces']))
    total_area_reservoir = porosity.dot(fp['areas'])
    saturation: np.ndarray = bc['initial_saturation']['value'].copy() 
    fp.insert_or_update_data({'sat_for_weight': saturation.copy()})
    pressure = np.repeat(0.0, fp['faces'].shape[0])
    newS = saturation.copy()
    saturation_plot = saturation.copy()
    
    dvtol = set_dvtol(fp, porosity, bc)

    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.create_tag('pressure')
    mesh_data.create_tag('faces_flux')
    mesh_data.create_tag('water_faces_flux')
    mesh_data.create_tag('saturation')
    
    pressures_update = []
    sat_update = []
    vpi_update = []

    initial_fine_vols = define_new_fine_levels_v1(fp, bc)
    
    gdict.update({'funcname': 'initial_loop'})

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
        cfl,
        refine_by_grad_bool,
        refine_by_estimator1_bool,
        max_value_grad,
        max_value_estimator1,
        vpis_to_plot,
        iterative_ms,
        tol_iterative,
        **gdict
    )
    
    TimeProfile.n_coarses = len(coarse_struct)
    TimeProfile.reset_data()

    print(f'LOOP: {loop} \n')
    loop += 1
    gdict.update({'funcname': 'while_loop', 'funcname_cum': 'while_loop_cum'})
    while vpi < max_vpi and loop < max_loop:
        p_updates = 0
        s_updates = 0
        t0 = time.perf_counter()
        pressure[:], edges_flux = update_pressure_only_ms(
            relative_perm,
            biphasic_mobility,
            saturation,
            fp,
            bc,
            lsds,
            initial_fine_vols,
            alpha_lim_finescale,
            beta_lim,
            OP,
            OR,
            coarse_struct,
            iterative_ms,
            tol_iterative,
            pressure,
            **gdict
        )
        t1 = time.perf_counter()
        TimeProfile.dt_total_pressure_update = t1 - t0
        TimeProfile.update_data()
        TimeProfile.show_data()
        if loop > 1:
            TimeProfile.export_to_data(**gdict)
            import pdb; pdb.set_trace()
        TimeProfile.reset_data()
        
        p_updates += 1
        
        fp.insert_or_update_data({
            'edges_flux1': edges_flux
        })
        
        dtnew = calculate_time_step_mimpes(fp, lsds, dvtol, Rdtmax, Rdtmin)
        dt_total = 0
        while dt_total < dtnew:
            dtmax = dtnew - dt_total   
            saturation_plot[:] = saturation         
            saturation[:], dt, plot_vpi, vpi, cumulative_oil, cumulative_water, water_flux, oil_flux, faces_flux, water_faces_flux = update_saturation_only(
                edges_flux,
                relative_perm,
                biphasic_mobility,
                saturation,
                fp,
                bc,
                lsds,
                porosity,
                total_area_reservoir,
                vpi,
                cumulative_oil,
                cumulative_water,
                dtmax,
                cfl,
                vpis_to_plot=vpis_to_plot
            )
            s_updates += 1
            dt_total += dt
            print(f"Loop: {loop}")
            print(f"Dt: {dt_total} / Dt_total: {dtnew}")
            print(f"VPI: {vpi}")
            loop+=1
            
            if plot_vpi == True:
                plot_results_ms(
                    loop,
                    pressure,
                    saturation,
                    simulation_data,
                    vpi,
                    cumulative_oil,
                    cumulative_water,
                    water_flux,
                    oil_flux,
                    fp,
                    fine_mesh_path,
                    mesh_data,
                    coarse_struct,
                    faces_flux,
                    water_faces_flux,
                    saturation_plot
                )
        

        fp.insert_or_update_data({
                'dt1': dtnew,
                'edges_flux0': fp['edges_flux1'].copy()
            })
        
        pressures_update.append(p_updates)
        sat_update.append(s_updates)
        vpi_update.append(vpi[0])
    
    vpi_update = np.array(vpi_update)
    pressures_update = np.array(pressures_update)
    sat_update = np.array(sat_update)
    ratio_updates = sat_update/pressures_update
    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot()
    ax1.plot(
        vpi_update,
        ratio_updates,
        marker='o',
        color='blue'
    )
    
    ax1.set_yscale('log')
    ax1.set_xlabel('PVI')
    ax1.set_ylabel('Saturation Update/Pressure Update')
    
    fig.savefig('pressure_saturation_updates_ameba_coarse4_ms.png', dpi=500)
    
    
        