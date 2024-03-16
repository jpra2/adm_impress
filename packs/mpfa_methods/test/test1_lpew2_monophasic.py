
from packs import defpaths
import numpy as np
from packs.mpfa_methods.mesh_preprocess import preprocess_mesh
import os
from packs.manager.meshmanager import MeshProperty
from packs.manager.mesh_data import MeshData
from packs.manager.boundary_conditions import BoundaryConditions
from packs.mpfa_methods.weight_interpolation.lpew import get_lpew2_weights
from packs.mpfa_methods.flux_calculation.diamond_method import DiamondFluxCalculation, get_xi_params_ds_flux
from scipy.sparse.linalg import spsolve
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
import sympy

"""
teste de convergencia problema 4.5 da dissertacao do professor Fernando
Escoamento Monofasico em um Dominio Heterogeneo e Suavemente Anisotropico
"""

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

def get_permeability_3(centroids):
    epsilon = 5e-2
    x = centroids[:, 0]
    y = centroids[:, 1]
    permeability = np.zeros((len(centroids), 2, 2))

    permeability[:, 0, 0] = np.power(y, 2) + epsilon*np.power(x, 2)
    permeability[:, 0, 1] = -(1-epsilon)*x*y
    permeability[:, 1, 0] = -(1-epsilon)*x*y
    permeability[:, 1, 1] = epsilon*np.power(y, 2) + np.power(x, 2)

    return {'permeability': permeability}

def get_permeability_4(centroids):
    x, y, kxx, kxy, kyx, kyy, uxy, fonte = problem_4_sympy()

    x_ = centroids[:, 0]
    y_ = centroids[:, 1]

    permeability = np.zeros((len(centroids), 2, 2))

    f = sympy.lambdify((x, y), kxx, 'numpy')
    
    permeability[:, 0, 0] = f(x_, y_)

    f = sympy.lambdify((x, y), kxy, 'numpy')

    permeability[:, 0, 1] = f(x_, y_)

    f = sympy.lambdify((x, y), kyx, 'numpy')

    permeability[:, 1, 0] = f(x_, y_)

    f = sympy.lambdify((x, y), kyy, 'numpy')

    permeability[:, 1, 1] = f(x_, y_)

    return {'permeability': permeability}

def problem_4_sympy():
    """
        Test 5: heterogeneous rotating anisotropy
        paper: Benchmark on Discretization Schemes
               for Anisotropic Diffusion Problems
               on General Grids

    """
    alpha = 1e-3
    x, y = sympy.symbols('x y')
    v1 = 1/(x**2 + y**2)
    kxx = v1*(alpha*(x**2) + y**2)
    kxy = v1*((alpha-1)*x*y)
    kyx = kxy
    kyy = v1*(x**2 + alpha*(y**2))


    k = [
        [kxx, kxy], 
        [kyx, kyy]
    ]

    uxy = sympy.sin(sympy.pi*x)*sympy.sin(sympy.pi*y)

    grad_uxy = [uxy.diff(x), uxy.diff(y)]

    k_grad_uxy = [
        k[0][0]*grad_uxy[0] + k[0][1]*grad_uxy[1],
        k[1][0]*grad_uxy[0] + k[1][1]*grad_uxy[1]
    ]

    # k_grad_u_x = kxx*grad_uxy_x + kxy*grad_uxy_y
    # k_grad_u_y = kyx*grad_uxy_x + kyy*grad_uxy_y

    fonte = -1*(k_grad_uxy[0].diff(x) + k_grad_uxy[1].diff(y))

    # fonte = -1*(k_grad_u_x.diff(x) + k_grad_u_y.diff(y))

    return x, y, kxx, kxy, kyx, kyy, uxy, fonte

def exact_solution_2(centroids, alpha):
    x = centroids[:,0]
    y = centroids[:,1]
    return np.exp(x*y)

def exact_solution_4(centroids):
    x, y, kxx, kxy, kyx, kyy, uxy, fonte = problem_4_sympy()

    x_ = centroids[:, 0]
    y_ = centroids[:, 1]
    f = sympy.lambdify((x, y), uxy, 'numpy')
    return f(x_, y_)


def get_source(centroids, alpha):
    x = centroids[:, 0]
    y = centroids[:, 1]
    source = np.zeros(len(x))

    test = x <= 0
    source[test] = alpha*(2*np.sin(y[test]) + np.cos(y[test]))*x[test] + np.sin(y[test])

    test = ~test
    source[test] = -2*alpha*np.exp(x[test])*np.cos(y[test])
    return source

def get_source_2(centroids, alpha):
    x = centroids[:, 0]
    y = centroids[:, 1]

    source = -2*np.exp(x*y)*(1 + np.power(y, 2) + np.power(x, 2) + x*y)
    return source

def get_source_3(centroids):
    x = centroids[:, 0]
    y = centroids[:, 1]

    mesh_delta = 1e-7

    test1 = x > 3/8 - mesh_delta
    test2 = x < 5/8 + mesh_delta
    test3 = y > 3/8 - mesh_delta
    test4 = y < 5/8 + mesh_delta

    test = test1 & test2 & test3 & test4
    source = np.zeros(len(x))
    source[test] = 1.0
    return source

def get_source_4(centroids):
    x, y, kxx, kxy, kyx, kyy, uxy, fonte = problem_4_sympy()

    f = sympy.lambdify((x, y), fonte, 'numpy')
    x_ = centroids[:, 0]
    y_ = centroids[:, 1]
    return f(x_, y_)






 
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

def define_boundary_conditions_3(mesh_properties: MeshProperty):
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
    
    pressures_bc = np.zeros(len(nodes_bc))
    bc.set_boundary('dirichlet_nodes', nodes_bc, pressures_bc)
    bc.set_boundary('neumann_edges', np.array([]), np.array([]))

    return bc

def define_boundary_conditions_4(mesh_properties: MeshProperty):
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
    
    pressures_bc = exact_solution_4(mesh_properties.nodes_centroids[nodes_bc])
    bc.set_boundary('dirichlet_nodes', nodes_bc, pressures_bc)
    bc.set_boundary('neumann_edges', np.array([]), np.array([]))

    return bc

def get_Emax(pressure, exact_sol):
    return np.absolute(pressure - exact_sol).max()

def get_Erms(pressure, exact_sol):
    dif = pressure - exact_sol
    n = len(pressure)
    return np.sqrt(np.sum(np.power(dif, 2))/n)

def get_L2(pressure, exact_sol):
    dif = pressure - exact_sol
    return np.linalg.norm(dif)

def get_Eu(pressure, exact_sol, areas):
    error = pressure - exact_sol
    error_2 = np.power(error, 2)
    resp = areas*error_2
    resp = np.sqrt(resp.sum())
    return resp

def get_erl2(pressure, exact_sol, areas):
    v1 = np.power(pressure - exact_sol, 2)*areas
    v2 = np.power(pressure, 2)*areas
    erl2 = np.sqrt(v1.sum()/v2.sum())
    return erl2


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
        resp['source'] += get_source(mesh_properties['faces_centroids'], alpha)*mesh_properties['areas']
        # resp['source'] += get_source(mesh_properties['faces_centroids'], alpha)


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
        resp['source'] += get_source(mesh_properties['faces_centroids'], alpha)*mesh_properties['areas']


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
        resp['source'] += get_source_2(mesh_properties['faces_centroids'], 1)*mesh_properties['areas']

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

    show_results(mesh_properties, alpha, mesh_name, [pressure_tag, error_plot_tag, exact_solution_tag])

def run_problem_3(mesh_name, mesh_properties_name):

    diamond_flux = DiamondFluxCalculation()

    mesh_properties = preprocess_mesh(mesh_name, mesh_properties_name)
    mesh_properties.insert_or_update_data(
        get_permeability_3(mesh_properties.faces_centroids)
    )
    bc = define_boundary_conditions_3(mesh_properties)
    mesh_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })
    get_lpew2_weights(mesh_properties)
    get_xi_params_ds_flux(mesh_properties)
    


    pressure_tag = 'pressure'

    if not mesh_properties.verify_name_in_data_names(pressure_tag):

        resp = diamond_flux.mount_problem(
            bc,
            **mesh_properties.get_all_data()
        )
        resp['source'] += get_source_3(mesh_properties['faces_centroids'])*mesh_properties['areas']

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

        mesh_properties.insert_or_update_data({
            pressure_tag: pressure
        })
        mesh_properties.export_data()

    print(f'MESH: {mesh_name}')
    print(f'N elements: {len(mesh_properties.faces)}')
    print(f'Pmax: {mesh_properties.pressure.max()}')
    print(f'Pmin: {mesh_properties.pressure.min()}')
    print()

    show_results(mesh_properties, 1, mesh_name, [pressure_tag])

    
def run_problem_3_lslds(mesh_name, mesh_properties_name):

    lsds = LsdsFluxCalculation()

    mesh_properties = preprocess_mesh(mesh_name, mesh_properties_name)
    mesh_properties.insert_or_update_data(
        get_permeability_3(mesh_properties.faces_centroids)
    )
    bc = define_boundary_conditions_3(mesh_properties)
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
    pressure_tag = lsds_prefix + 'pressure'

    if not mesh_properties.verify_name_in_data_names(pressure_tag):

        resp = lsds.mount_problem_v6(
            bc,
            **mesh_properties.get_all_data()
        )
        resp['source'] += get_source_3(mesh_properties['faces_centroids'])*mesh_properties['areas']


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
        mesh_properties.insert_or_update_data({
            pressure_tag: pressure
        })
        mesh_properties.export_data()
        
    print(f'MESH: {mesh_name}')
    print(f'N elements: {len(mesh_properties.faces)}')
    print(f'Pmax: {mesh_properties[pressure_tag].max()}')
    print(f'Pmin: {mesh_properties[pressure_tag].min()}')
    print()

def run_problem_4(mesh_name, mesh_properties_name):
    dtype_E = [('Emax', np.float64), ('Erms', np.float64), ('L2', np.float64),
               ('Eu', np.float64)]

    diamond_flux = DiamondFluxCalculation()

    mesh_properties = preprocess_mesh(mesh_name, mesh_properties_name)
    mesh_properties.insert_or_update_data(
        get_permeability_4(mesh_properties.faces_centroids)
    )
    bc = define_boundary_conditions_4(mesh_properties)
    mesh_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })
    get_lpew2_weights(mesh_properties, update=False)
    get_xi_params_ds_flux(mesh_properties, update=False)

    error_tag = 'error'
    pressure_tag = 'pressure'
    error_plot_tag = 'error_plot'
    exact_solution_tag = 'exact_solution'

    if not mesh_properties.verify_name_in_data_names(error_tag):

        resp = diamond_flux.mount_problem(
            bc,
            **mesh_properties.get_all_data()
        )
        resp['source'] += get_source_4(mesh_properties['faces_centroids'])*mesh_properties['areas']

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

        exact_faces_solution = exact_solution_4(mesh_properties.faces_centroids)

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
    
    erl2 = get_erl2(mesh_properties[pressure_tag], mesh_properties[exact_solution_tag], mesh_properties['areas'])
    print(f'MESH: {mesh_name}')
    print(f'Emax: {mesh_properties[error_tag]["Emax"]}')
    print(f'Erms: {mesh_properties[error_tag]["Erms"]}')
    print(f'L2: {mesh_properties[error_tag]["L2"]}')
    print(f'Eu: {mesh_properties[error_tag]["Eu"]}')
    print(f"N elements: {len(mesh_properties['faces'])}")
    print(f"Erl2: {erl2}")
    print()

def run_problem_4_lsds(mesh_name, mesh_properties_name):
    dtype_E = [('Emax', np.float64), ('Erms', np.float64), ('L2', np.float64),
               ('Eu', np.float64)]

    lsds = LsdsFluxCalculation()

    mesh_properties = preprocess_mesh(mesh_name, mesh_properties_name)
    mesh_properties.insert_or_update_data(
        get_permeability_4(mesh_properties.faces_centroids)
    )
    bc = define_boundary_conditions_4(mesh_properties)
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

    error_tag = 'error'
    pressure_tag = 'pressure'
    error_plot_tag = 'error_plot'
    exact_solution_tag = 'exact_solution'

    if not mesh_properties.verify_name_in_data_names(pressure_tag):

        resp = lsds.mount_problem_v6(
            bc,
            **mesh_properties.get_all_data()
        )
        resp['source'] += get_source_4(mesh_properties['faces_centroids'])*mesh_properties['areas']


        pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])

        exact_faces_solution = exact_solution_4(mesh_properties.faces_centroids)

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
        
    erl2 = get_erl2(mesh_properties[pressure_tag], mesh_properties[exact_solution_tag], mesh_properties['areas'])
    print(f'MESH: {mesh_name}')
    print(f'Emax: {mesh_properties[error_tag]["Emax"]}')
    print(f'Erms: {mesh_properties[error_tag]["Erms"]}')
    print(f'L2: {mesh_properties[error_tag]["L2"]}')
    print(f'Eu: {mesh_properties[error_tag]["Eu"]}')
    print(f"N elements: {len(mesh_properties['faces'])}")
    print(f"Erl2: {erl2}")
    print()




def test_problem1():
    mesh_prefix = 'str_trimesh_'
    mesh_sufix = '.msh'
    # list_of_meshs = ['16x16', '32x32', '64x64', '128x128']
    list_of_meshs = ['8x8', '16x16', '32x32', '64x64']
    alpha = 1000
    for n_mesh in list_of_meshs:
        mesh_prop_name = mesh_prefix + n_mesh
        mesh_name = mesh_prop_name + mesh_sufix
        mesh_name = os.path.join(defpaths.lpew2_mesh_folder, mesh_name)
        # run_problem(mesh_name, mesh_prop_name, alpha)
        # import pdb; pdb.set_trace()
        run_problem_lsds(mesh_name, mesh_prop_name, alpha)

def test_problem2():
    mesh_prefix = 'uns_trimesh_'
    mesh_sufix = '.msh'
    # list_of_meshs = ['16x16', '32x32', '64x64', '128x128']
    list_of_meshs = ['16x16', '32x32', '64x64']
    alpha = 1
    for n_mesh in list_of_meshs:
        mesh_prop_name = mesh_prefix + n_mesh
        mesh_name = mesh_prop_name + mesh_sufix
        mesh_name = os.path.join(defpaths.lpew2_mesh_folder, mesh_name)
        run_problem_2(mesh_name, mesh_prop_name, alpha)
        # run_problem_lsds(mesh_name, mesh_prop_name, alpha)



def test_problem3():
    mesh_prefix = 'uns_trimesh_'
    mesh_sufix = '.msh'
    # list_of_meshs = ['16x16', '32x32', '64x64', '128x128']
    list_of_meshs = ['8x8', '16x16', '32x32', '64x64']
    for n_mesh in list_of_meshs:
        mesh_prop_name = mesh_prefix + n_mesh
        mesh_name = mesh_prop_name + mesh_sufix
        mesh_name = os.path.join(defpaths.lpew2_mesh_folder, mesh_name)
        run_problem_3(mesh_name, mesh_prop_name)
        # run_problem_3_lslds(mesh_name, mesh_prop_name)
    
def test_problem4():
    # mesh_prefix = 'mesh1_'
    mesh_prefix = 'str_quadmesh_'
    # mesh_prefix = 'str2_trimesh_'
    # mesh_prefix = 'uns_trimesh_'
    mesh_sufix = '.msh'
    # list_of_meshs = ['16x16', '32x32', '64x64', '128x128']
    list_of_meshs = ['8x8', '16x16', '32x32', '64x64', '128x128']
    # list_of_meshs = ['8_8', '16_16', '32_32', '64_64', '128_128']
    for n_mesh in list_of_meshs:
        mesh_prop_name = mesh_prefix + n_mesh
        mesh_name = mesh_prop_name + mesh_sufix
        # mesh_name = mesh_prefix + n_mesh + mesh_sufix
        # mesh_name = os.path.join(defpaths.mpfad_mesh_folder, mesh_name)
        mesh_name = os.path.join(defpaths.lpew2_mesh_folder, mesh_name)
        run_problem_4(mesh_name, mesh_prop_name)
        # import pdb; pdb.set_trace()
        # mesh_prop_name = 'lsds_' + mesh_prop_name
        # run_problem_4_lsds(mesh_name, mesh_prop_name)
        # import pdb; pdb.set_trace()