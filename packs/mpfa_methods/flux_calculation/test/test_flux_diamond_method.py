from files_to_test.lpew2_test_weights.info_test1 import get_filenames
from packs.mpfa_methods.weight_interpolation.test.test_lpew_weights import initialize_mesh
from packs.manager.meshmanager import MeshProperty, load_mesh_properties
from packs.mpfa_methods.flux_calculation.diamond_method import get_xi_params_ds_flux, DiamondFluxCalculation
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess
import numpy as np
from packs.manager.boundary_conditions import BoundaryConditions
from packs.mpfa_methods.weight_interpolation.lpew import LpewWeight, get_lpew2_weights
from scipy.sparse.linalg import spsolve
from packs import defnames, defpaths
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists, mesh_verify, create_properties
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.manager.mesh_data import MeshData
import os
import matplotlib.pyplot as plt
import pandas as pd

def get_permeability(faces):
    permeability = np.zeros((len(faces), 2, 2))
    permeability[:, 1, 1] = 1
    permeability[:, 0, 0] = 1
    
    return {'permeability': permeability}

def get_edges_centroids(nodes_centroids, nodes_of_edges):
    resp = (nodes_centroids[nodes_of_edges[: 1]] + nodes_centroids[nodes_of_edges[: 0]])/2
    return resp

def get_dirichlet_nodes(nodes, nodes_centroids):
    
    k1 = 100
    k2 = 1
    
    delta = 1e-10
    nodes1 = nodes[nodes_centroids[:, 0] <= 0 + delta] 
    nodes2 = nodes[nodes_centroids[:, 0] >= nodes_centroids[:, 0].max() - delta] 
    
    values1 = np.repeat(k1, len(nodes1))
    values2 = np.repeat(k2, len(nodes2))
    
    dr_nodes = np.concatenate([nodes1, nodes2])
    values_dirichlet = np.concatenate([values1, values2])
    
    return dr_nodes, values_dirichlet

def get_neumann_edges(edges, edges_centroids):
    k1 = 0
    k2 = 0
    
    delta = 1e-10
    edges1 = edges[edges_centroids[:, 1] <= 0 + delta] 
    edges2 = edges[edges_centroids[:, 1] >= edges_centroids[:, 1].max() - delta] 
    
    values1 = np.repeat(k1, len(edges1))
    values2 = np.repeat(k2, len(edges2))
    
    neumann_edges = np.concatenate([edges1, edges2])
    values_neumann = np.concatenate([values1, values2])
    
    return neumann_edges, values_neumann


def setup5():
    mesh_test_name = defpaths.linear_2k_test
    mesh_properties_name = defpaths.mesh_prop_linear_2k
    mesh_verify(mesh_test_name)
    mesh_properties: MeshProperty = create_properties_if_not_exists(mesh_test_name, mesh_properties_name)
    
    mpfapreprocess = MpfaPreprocess()
    mpfapreprocess.calculate_h_dist(mesh_properties)
    mpfapreprocess.calculate_areas(mesh_properties)

    def get_permeability_linear(centroids):
        x = centroids[:, 0]
        n = len(x)

        k1 = 10
        k2 = 1

        test = x <= 0.5
        test2 = ~test

        permeability = np.zeros((n, 2, 2))
        permeability[:, 0, 0] = 1
        permeability[:, 1, 1] = 1

        permeability[test] = permeability[test]*k1
        permeability[test2] = permeability[test2]*k2

        return {'permeability': permeability}

    mesh_properties.insert_or_update_data(
        get_permeability_linear(
            mesh_properties.faces_centroids[:, 0:2]
        )
    )

    #define boundary conditions    
    problem_name = 'linear_2k_test'
    bc = BoundaryConditions()
    bc.insert_name(problem_name)
    
    xmin, ymin = mesh_properties.faces_centroids[:, 0:2].min(axis=0)
    xmax, ymax = mesh_properties.faces_centroids[:, 0:2].max(axis=0)
    
    mesh_delta = 1e-7
    nodes_xmin = mesh_properties.nodes[
        mesh_properties.nodes_centroids[:, 0] < xmin + mesh_delta
    ]
    nodes_xmax = mesh_properties.nodes[
        mesh_properties.nodes_centroids[:, 0] > xmax - mesh_delta
    ]
    # nodes_ymin = mesh_properties.nodes[
    #     mesh_properties.nodes_centroids[:, 1] < ymin + mesh_delta
    # ]
    # nodes_ymax = mesh_properties.nodes[
    #     mesh_properties.nodes_centroids[:, 1] > ymax - mesh_delta
    # ]
    
    nodes_bc = np.concatenate([
        nodes_xmin, nodes_xmax
    ]).astype(np.uint64)

    pressures_xmin = np.repeat(100.0, len(nodes_xmin))
    pressures_xmax = np.repeat(1.0, len(nodes_xmax))
    
    pressures_bc = np.concatenate([pressures_xmin, pressures_xmax])
    
    bc.set_boundary('dirichlet_nodes', nodes_bc, pressures_bc)

    edges_centroids = mesh_properties.edges_centroids
    edges_ymax = mesh_properties.edges[edges_centroids[:, 1] >= edges_centroids[:, 1].max() - mesh_delta]
    edges_ymin = mesh_properties.edges[edges_centroids[:, 1] <= 0 + mesh_delta]

    edges_bc = np.concatenate([edges_ymax, edges_ymin])
    neumann_values = np.repeat(0.0, len(edges_bc))

    bc.set_boundary('neumann_edges', edges_bc, neumann_values)

    mesh_properties.insert_or_update_data({
        'neumann_edges': edges_bc,
        'neumann_edges_value': neumann_values
    })
    
    ## create weights and xi params ds for flux calculation
    lsds = LsdsFluxCalculation()
    mesh_properties.insert_or_update_data(
        lsds.preprocess(mesh_properties)
    )

    if not mesh_properties.verify_name_in_data_names('nodes_weights'):
        # define nodes to calculate_weights
        mesh_properties.insert_data({'nodes_to_calculate': mesh_properties.nodes.copy()})
        get_lpew2_weights(mesh_properties)
        mesh_properties.remove_data(['nodes_to_calculate'])
        mesh_properties.export_data()

    if not mesh_properties.verify_name_in_data_names('xi_params_dsflux'):
        get_xi_params_ds_flux(mesh_properties)

    diamond_flux = DiamondFluxCalculation()
    
    # resp = lsds.mount_problem(
    #     mesh_properties.nodes_weights,
    #     mesh_properties.xi_params,
    #     mesh_properties.faces,
    #     mesh_properties.edges,
    #     mesh_properties.bool_boundary_edges,
    #     mesh_properties.adjacencies,
    #     bc,
    #     mesh_properties.nodes_of_edges,
    #     mesh_properties.neumann_weights
    # )

    # resp = lsds.mount_problem_v2(
    #     mesh_properties.nodes_weights,
    #     mesh_properties.xi_params,
    #     mesh_properties.faces,
    #     mesh_properties.nodes,
    #     mesh_properties.bool_boundary_edges,
    #     mesh_properties.adjacencies,
    #     bc,
    #     mesh_properties.nodes_of_edges,
    #     mesh_properties.neumann_weights,
    #     mesh_properties.edges_of_nodes
    # )

    # resp = lsds.mount_problem_v3(bc, **mesh_properties.get_all_data())
    resp = diamond_flux.mount_problem(bc, **mesh_properties.get_all_data())

    pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])
    edges_flux = lsds.get_edges_flux(
        mesh_properties.xi_params_dsflux,
        mesh_properties.nodes_weights,
        mesh_properties.nodes_of_edges,
        pressure,
        mesh_properties.adjacencies,
        bc,
        mesh_properties.neumann_weights        
    )
    
    faces_flux = lsds.get_faces_flux(
        edges_flux,
        mesh_properties.adjacencies,
        mesh_properties.bool_boundary_edges
    )
    
    mesh_path = os.path.join(defpaths.mesh, mesh_test_name)
    mesh_data = MeshData(mesh_path=mesh_path)
    mesh_data.create_tag('pressure')
    mesh_data.insert_tag_data('pressure', pressure, 'faces', mesh_properties.faces)
    mesh_data.create_tag('permeability')
    perm = mesh_properties.permeability.reshape((len(mesh_properties.faces), 4))
    mesh_data.insert_tag_data('permeability', perm[:, 0], 'faces', mesh_properties.faces)
    mesh_data.create_tag('faces_flux')
    mesh_data.insert_tag_data('faces_flux', faces_flux, 'faces', mesh_properties.faces)
    mesh_data.create_tag('nodes_pressure_presc')
    mesh_data.insert_tag_data('nodes_pressure_presc', pressures_bc, 'nodes', nodes_bc)
    # mesh_data.create_tag('pressure_error')
    # mesh_data.insert_tag_data('pressure_error', error, 'faces', mesh_properties.faces)
    mesh_data.export_all_elements_type_to_vtk('linear_test_nodes', 'nodes')
    mesh_data.export_all_elements_type_to_vtk('linear_test_faces', 'faces')
    mesh_data.export_only_the_elements('linear_test_nodes_pressure_boundary', 'nodes', nodes_bc)
    
    plt.clf()
    fig1, ax1 = plt.subplots(1)
    facesxcentroids = mesh_properties.faces_centroids[:, 0]
    pd_series = pd.Series(facesxcentroids, index=np.arange(len(facesxcentroids)))
    pd_series_ordenated = pd_series.sort_values(ascending=True)
    
    ax1.plot(facesxcentroids[pd_series_ordenated.index], pressure[pd_series_ordenated.index])
    ax1.set_xlabel('X position')
    ax1.set_ylabel('Presssure')
    ax1.set_title('Pressure x position')
    # ax1.legend()

    ext = 'svg'
    fig1.savefig(os.path.join('results', 'LinearTest.' + ext), format=ext)

    
    import pdb; pdb.set_trace()


def test_xi_params_ds():
    setup5()
    # filenames = get_filenames()
    # initialize_mesh(**filenames)
    # mesh_properties: MeshProperty = load_mesh_properties(filenames['mesh_properties_name'])    
    # get_lpew2_weights(mesh_properties)
    # get_xi_params_ds_flux(mesh_properties)
    
    
    # bc = BoundaryConditions()
    # dr_nodes, values_dirichlet = get_dirichlet_nodes(mesh_properties.edges, mesh_properties.edges_centroids)
    # neumann_edges, values_neumann = get_neumann_edges(mesh_properties.edges, mesh_properties.edges_centroids)
    
    # bc.set_boundary('dirichlet_nodes', dr_nodes, values_dirichlet)
    # bc.set_boundary('neumann_edges', neumann_edges, values_neumann)
    
    # lpew2 = LpewWeight()
    # neumann_weights = lpew2.create_lpew2_neumann_weights(
    #     edges_neumann=neumann_edges,
    #     neumann_values=values_neumann,
    #     **mesh_properties.get_all_data()        
    # )
    
    # mesh_properties.insert_or_update_data(neumann_weights)
    
    # diamondflux = DiamondFluxCalculation()
    # resp = diamondflux.mount_problem(bc, **mesh_properties.get_all_data())

    # pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])
    
    
    
    
    import pdb; pdb.set_trace()
    
    
    
    

    



