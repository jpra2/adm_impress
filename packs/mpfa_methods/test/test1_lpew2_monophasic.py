from packs import defpaths
import numpy as np
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess
import os
from packs.manager.meshmanager import MeshProperty
from packs.manager.boundary_conditions import BoundaryConditions
from packs.mpfa_methods.weight_interpolation.lpew import get_lpew2_weights
from packs.mpfa_methods.flux_calculation.diamond_method import DiamondFluxCalculation, get_xi_params_ds_flux
from scipy.sparse.linalg import spsolve

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
    resp[test] = (2*np.sin(y[test]) + np.cos(y[test]))*alpha*x[test] + np.sin(y[test])
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

def define_boundary_conditions(mesh_properties: MeshProperty, alpha):
    bc = BoundaryConditions()

    xmin, ymin = mesh_properties.faces_centroids[:, 0:2].min(axis=0)
    xmax, ymax = mesh_properties.faces_centroids[:, 0:2].max(axis=0)
    
    mesh_delta = 1e-5
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

def get_Emax(pressure, exact_sol):
    return np.absolute(pressure - exact_sol).max()

def get_Erms(pressure, exact_sol):
    dif = pressure - exact_sol
    n = len(pressure)
    return np.linalg.norm(dif/n)

def run_problem(mesh_name, mesh_properties_name, alpha):

    dtype_E = [('Emax', np.float64), ('Erms', np.float64)]

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

    if not mesh_properties.verify_name_in_data_names('error'):

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

        Emax = get_Emax(pressure, exact_faces_solution)
        Erms = get_Erms(pressure, exact_faces_solution)

        array = np.zeros(1, dtype=dtype_E)
        array['Emax'] = Emax
        array['Erms'] = Erms

        print(f'Emax: {Emax}')
        print(f'Erms: {Erms}')

        mesh_properties.insert_or_update_data({
            'error': array
        })
        mesh_properties.export_data()

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



