from packs import defpaths
import numpy as np
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess
import os
from packs.manager.meshmanager import MeshProperty
from packs.manager.mesh_data import MeshData
from packs.manager.boundary_conditions import BoundaryConditions
from packs.mpfa_methods.weight_interpolation.lpew import get_lpew2_weights
from packs.mpfa_methods.flux_calculation.diamond_method import DiamondFluxCalculation, get_xi_params_ds_flux
from scipy.sparse.linalg import spsolve
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation

"""
teste de convergencia problema 4.5 da dissertacao do professor Fernando
Escoamento Monofasico em um Dominio Heterogeneo e Suavemente Anisotropico
"""

def preprocess_mesh(mesh_name, mesh_properties_name) -> MeshProperty:
    mpfa_preprocess = MpfaPreprocess()
    mesh_properties = mpfa_preprocess.create_properties_if_not_exists(mesh_name, mesh_properties_name)
    mpfa_preprocess.preprocess_data_lsds(mesh_properties)
    mpfa_preprocess.calculate_areas(mesh_properties)
    mpfa_preprocess.calculate_h_dist(mesh_properties)

    return mesh_properties

def exact_solution(centroids, alpha):
    x = centroids[:, 0]
    y = centroids[:, 1]

    resp = np.zeros(len(x))

    test = x <= 0
    resp[test] = alpha*(2*np.sin(y[test]) + np.cos(y[test]))*x[test] + np.sin(y[test])
    test2 = ~test
    resp[test2] = np.exp(x[test2])*np.sin(y[test2])

    return resp 

def get_permeability(centroids, alpha):
    x = centroids[:, 0]
    permeability = np.zeros((len(x), 2, 2))
    test = x < 0

    k1 = np.array([
        [1, 0],
        [0, 1]
    ])

    k2 = alpha*np.array([
        [2, 1],
        [1, 2]
    ])
    permeability[test] = k1
    test2 = ~test
    permeability[test2] = k2

    return {'permeability': permeability}

def get_permeability_2(centroids, alpha):

    k = np.array([
        [2, 1],
        [1, 2]
    ])

    permeability = np.zeros((len(centroids), 2, 2))
    permeability[:] = k
    return {'permeability': permeability}

def exact_solution_2(centroids, alpha):
    x = centroids[:,0]
    y = centroids[:,1]
    return np.exp(x*y)

def get_source(centroids, alpha):
    x = centroids[:, 0]
    y = centroids[:, 1]
    source = np.zeros(len(x))

    test = x <= 0
    source[test] = alpha*(2*np.sin(x[test]) + np.cos(y[test]))*x[test] + np.sin(y[test])

    test = ~test
    source[test] = -2*alpha*np.exp(x[test])*np.cos(y[test])
    return source

def get_source_2(centroids, alpha):
    x = centroids[:, 0]
    y = centroids[:, 1]

    source = -2*np.exp(x*y)*(1 + np.power(y, 2) + np.power(x, 2) + x*y)
    return source

def define_boundary_conditions(mesh_properties: MeshProperty, alpha):
    bc = BoundaryConditions()

    xmin, ymin = mesh_properties.faces_centroids[:, 0:2].min(axis=0)
    xmax, ymax = mesh_properties.faces_centroids[:, 0:2].max(axis=0)
    
    mesh_delta = 1e-7
    nodes_xmin = mesh_properties.nodes[
        mesh_properties.nodes_centroids[:, 0] < xmin + mesh_delta
    ]
    nodes_xmax = mesh_properties.nodes[
        mesh_properties.nodes_centroids[:, 0] > xmax - mesh_delta
    ]
    nodes_ymin = mesh_properties.nodes[
        mesh_properties.nodes_centroids[:, 1] < ymin + mesh_delta
    ]
    nodes_ymax = mesh_properties.nodes[
        mesh_properties.nodes_centroids[:, 1] > ymax - mesh_delta
    ]
    
    nodes_bc = np.unique(np.concatenate([
        nodes_xmin, nodes_xmax, nodes_ymin, nodes_ymax
    ])).astype(np.uint64)
    
    pressures_bc = exact_solution(mesh_properties.nodes_centroids[nodes_bc], alpha)
    bc.set_boundary('dirichlet_nodes', nodes_bc, pressures_bc)
    bc.set_boundary('neumann_edges', np.array([]), np.array([]))

    return bc

def define_boundary_conditions_2(mesh_properties: MeshProperty, alpha):
    bc = BoundaryConditions()

    xmin, ymin = mesh_properties.faces_centroids[:, 0:2].min(axis=0)
    xmax, ymax = mesh_properties.faces_centroids[:, 0:2].max(axis=0)
    
    mesh_delta = 1e-7
    nodes_xmin = mesh_properties.nodes[
        mesh_properties.nodes_centroids[:, 0] < xmin + mesh_delta
    ]
    nodes_xmax = mesh_properties.nodes[
        mesh_properties.nodes_centroids[:, 0] > xmax - mesh_delta
    ]
    nodes_ymin = mesh_properties.nodes[
        mesh_properties.nodes_centroids[:, 1] < ymin + mesh_delta
    ]
    nodes_ymax = mesh_properties.nodes[
        mesh_properties.nodes_centroids[:, 1] > ymax - mesh_delta
    ]
    
    nodes_bc = np.unique(np.concatenate([
        nodes_xmin, nodes_xmax, nodes_ymin, nodes_ymax
    ])).astype(np.uint64)
    
    pressures_bc = exact_solution_2(mesh_properties.nodes_centroids[nodes_bc], alpha)
    bc.set_boundary('dirichlet_nodes', nodes_bc, pressures_bc)
    bc.set_boundary('neumann_edges', np.array([]), np.array([]))

    return bc

def get_Emax(pressure, exact_sol):
    return np.absolute(pressure - exact_sol).max()

def get_Erms(pressure, exact_sol):
    dif = pressure - exact_sol
    n = len(pressure)
    return np.linalg.norm(dif/n)

def get_L2(pressure, exact_sol):
    dif = pressure - exact_sol
    return np.linalg.norm(dif)

def get_Eu(pressure, exact_sol, areas):
    error = pressure - exact_sol
    error_2 = np.power(error, 2)
    resp = areas*error_2
    resp = np.sqrt(resp.sum())
    return resp

def show_results(mesh_properties: MeshProperty, alpha, mesh_path, tags_to_plot):
    export_name = 'result_' + mesh_properties['mesh_name'][0] + '_' + str(alpha)
    mesh_data = MeshData(mesh_path=mesh_path)
    for tag in tags_to_plot:
        mesh_data.create_tag(tag)
        mesh_data.insert_tag_data(tag, mesh_properties[tag], elements_type='faces', elements_array=mesh_properties['faces'])
    
    mesh_data.export_all_elements_type_to_vtk(export_name, element_type='faces')





def run_problem(mesh_name, mesh_properties_name, alpha):

    dtype_E = [('Emax', np.float64), ('Erms', np.float64), ('L2', np.float64),
               ('Eu', np.float64)]

    diamond_flux = DiamondFluxCalculation()

    mesh_properties = preprocess_mesh(mesh_name, mesh_properties_name)
    mesh_properties.insert_or_update_data(
        get_permeability(mesh_properties.faces_centroids, alpha)
    )
    bc = define_boundary_conditions(mesh_properties, alpha)
    mesh_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })
    get_lpew2_weights(mesh_properties)
    get_xi_params_ds_flux(mesh_properties)

    error_tag = 'error_alpha_' + str(alpha)
    pressure_tag = 'pressure_alpha_' + str(alpha)
    error_plot_tag = 'error_plot_' + str(alpha)
    exact_solution_tag = 'exact_solution_' + str(alpha)

    if not mesh_properties.verify_name_in_data_names(error_tag):

        resp = diamond_flux.mount_problem(
            bc,
            **mesh_properties.get_all_data()
        )

        pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])

        # edges_flux = diamond_flux.get_edges_flux(
        #     bc,
        #     pressure,
        #     **mesh_properties.get_all_data()        
        # )
        
        # faces_flux = diamond_flux.get_faces_flux(
        #     edges_flux,
        #     mesh_properties.adjacencies,
        #     mesh_properties.bool_boundary_edges
        # )

        exact_faces_solution = exact_solution(mesh_properties.faces_centroids, alpha)

        array = np.zeros(1, dtype=dtype_E)
        array['Emax'] = get_Emax(pressure, exact_faces_solution)
        array['Erms'] = get_Erms(pressure, exact_faces_solution)
        array['L2'] = get_L2(pressure, exact_faces_solution)
        array['Eu'] = get_Eu(pressure, exact_faces_solution, mesh_properties['areas'])

        mesh_properties.insert_or_update_data({
            error_tag: array,
            pressure_tag: pressure,
            error_plot_tag: np.absolute(pressure - exact_faces_solution),
            exact_solution_tag: exact_faces_solution
        })
        mesh_properties.export_data()
        
    print(f'MESH: {mesh_name}')
    print(f'Emax: {mesh_properties[error_tag]["Emax"]}')
    print(f'Erms: {mesh_properties[error_tag]["Erms"]}')
    print(f'L2: {mesh_properties[error_tag]["L2"]}')
    print(f'Eu: {mesh_properties[error_tag]["Eu"]}')
    print()

    show_results(mesh_properties, alpha, mesh_name, [pressure_tag, error_plot_tag, exact_solution_tag])

def run_problem_lsds(mesh_name, mesh_properties_name, alpha):
    dtype_E = [('Emax', np.float64), ('Erms', np.float64), ('L2', np.float64),
               ('Eu', np.float64)]

    lsds = LsdsFluxCalculation()

    mesh_properties = preprocess_mesh(mesh_name, mesh_properties_name)
    mesh_properties.insert_or_update_data(
        get_permeability(mesh_properties.faces_centroids, alpha)
    )
    bc = define_boundary_conditions(mesh_properties, alpha)
    mesh_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })
    mesh_properties.insert_or_update_data(
        {'nodes_to_calculate': mesh_properties['nodes']}
    )
    mesh_properties.insert_or_update_data(
        get_gls_nodes_weights(**mesh_properties.get_all_data())
    )

    mesh_properties.insert_or_update_data(
        lsds.get_all_edges_flux_params(**mesh_properties.get_all_data())
    )

    lsds_prefix = 'lsds_'
    error_tag = lsds_prefix + 'error_alpha_' + str(alpha)
    pressure_tag = lsds_prefix + 'pressure_alpha_' + str(alpha)
    error_plot_tag = lsds_prefix + 'error_plot_' + str(alpha)
    exact_solution_tag = lsds_prefix + 'exact_solution_' + str(alpha)

    if not mesh_properties.verify_name_in_data_names(error_tag):

        resp = lsds.mount_problem_v6(
            bc,
            **mesh_properties.get_all_data()
        )

        pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])

        # edges_flux = diamond_flux.get_edges_flux(
        #     bc,
        #     pressure,
        #     **mesh_properties.get_all_data()        
        # )
        
        # faces_flux = diamond_flux.get_faces_flux(
        #     edges_flux,
        #     mesh_properties.adjacencies,
        #     mesh_properties.bool_boundary_edges
        # )

        exact_faces_solution = exact_solution(mesh_properties.faces_centroids, alpha)

        array = np.zeros(1, dtype=dtype_E)
        array['Emax'] = get_Emax(pressure, exact_faces_solution)
        array['Erms'] = get_Erms(pressure, exact_faces_solution)
        array['L2'] = get_L2(pressure, exact_faces_solution)
        array['Eu'] = get_Eu(pressure, exact_faces_solution, mesh_properties['areas'])

        mesh_properties.insert_or_update_data({
            error_tag: array,
            pressure_tag: pressure,
            error_plot_tag: np.absolute(pressure - exact_faces_solution),
            exact_solution_tag: exact_faces_solution
        })
        mesh_properties.export_data()
        
    print(f'MESH: {mesh_name}')
    print(f'Emax: {mesh_properties[error_tag]["Emax"]}')
    print(f'Erms: {mesh_properties[error_tag]["Erms"]}')
    print(f'L2: {mesh_properties[error_tag]["L2"]}')
    print(f'Eu: {mesh_properties[error_tag]["Eu"]}')
    print()

    show_results(mesh_properties, alpha, mesh_name, [pressure_tag, error_plot_tag, exact_solution_tag])


def run_problem_2(mesh_name, mesh_properties_name, alpha):

    dtype_E = [('Emax', np.float64), ('Erms', np.float64), ('L2', np.float64),
               ('Eu', np.float64)]

    diamond_flux = DiamondFluxCalculation()

    mesh_properties = preprocess_mesh(mesh_name, mesh_properties_name)
    mesh_properties.insert_or_update_data(
        get_permeability_2(mesh_properties.faces_centroids, alpha)
    )
    bc = define_boundary_conditions_2(mesh_properties, alpha)
    mesh_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })
    get_lpew2_weights(mesh_properties)
    get_xi_params_ds_flux(mesh_properties)

    error_tag = 'error_alpha_' + str(alpha)
    pressure_tag = 'pressure_alpha_' + str(alpha)
    error_plot_tag = 'error_plot_' + str(alpha)
    exact_solution_tag = 'exact_solution_' + str(alpha)

    if not mesh_properties.verify_name_in_data_names(error_tag):

        resp = diamond_flux.mount_problem(
            bc,
            **mesh_properties.get_all_data()
        )
        resp['source'] += get_source_2(mesh_properties['faces_centroids'], 1)

        # import pdb; pdb.set_trace()

        pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])

        # edges_flux = diamond_flux.get_edges_flux(
        #     bc,
        #     pressure,
        #     **mesh_properties.get_all_data()        
        # )
        
        # faces_flux = diamond_flux.get_faces_flux(
        #     edges_flux,
        #     mesh_properties.adjacencies,
        #     mesh_properties.bool_boundary_edges
        # )

        exact_faces_solution = exact_solution_2(mesh_properties.faces_centroids, alpha)

        array = np.zeros(1, dtype=dtype_E)
        array['Emax'] = get_Emax(pressure, exact_faces_solution)
        array['Erms'] = get_Erms(pressure, exact_faces_solution)
        array['L2'] = get_L2(pressure, exact_faces_solution)
        array['Eu'] = get_Eu(pressure, exact_faces_solution, mesh_properties['areas'])

        mesh_properties.insert_or_update_data({
            error_tag: array,
            pressure_tag: pressure,
            error_plot_tag: np.absolute(pressure - exact_faces_solution),
            exact_solution_tag: exact_faces_solution
        })
        mesh_properties.export_data()
        
    print(f'MESH: {mesh_name}')
    print(f'Emax: {mesh_properties[error_tag]["Emax"]}')
    print(f'Erms: {mesh_properties[error_tag]["Erms"]}')
    print(f'L2: {mesh_properties[error_tag]["L2"]}')
    print(f'Eu: {mesh_properties[error_tag]["Eu"]}')
    print()
    import pdb; pdb.set_trace()

    show_results(mesh_properties, alpha, mesh_name, [pressure_tag, error_plot_tag, exact_solution_tag])


def test_problem1():
    mesh_prefix = 'str_trimesh_'
    mesh_sufix = '.msh'
    list_of_meshs = ['16x16', '32x32', '64x64', '128x128']
    alpha = 1
    for n_mesh in list_of_meshs:
        mesh_prop_name = mesh_prefix + n_mesh
        mesh_name = mesh_prop_name + mesh_sufix
        mesh_name = os.path.join(defpaths.lpew2_mesh_folder, mesh_name)
        run_problem(mesh_name, mesh_prop_name, alpha)
        # run_problem_lsds(mesh_name, mesh_prop_name, alpha)

def test_problem2():
    mesh_prefix = 'uns_trimesh_'
    mesh_sufix = '.msh'
    list_of_meshs = ['16x16', '32x32', '64x64', '128x128']
    alpha = 1
    for n_mesh in list_of_meshs:
        mesh_prop_name = mesh_prefix + n_mesh
        mesh_name = mesh_prop_name + mesh_sufix
        mesh_name = os.path.join(defpaths.lpew2_mesh_folder, mesh_name)
        run_problem_2(mesh_name, mesh_prop_name, alpha)
        # run_problem_lsds(mesh_name, mesh_prop_name, alpha)



