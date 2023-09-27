from files_to_test.lpew2_test_weights.info_test1 import get_filenames
from packs.mpfa_methods.weight_interpolation.test.test_lpew_weights import initialize_mesh
from packs.manager.meshmanager import MeshProperty, load_mesh_properties
from packs.mpfa_methods.flux_calculation.diamond_method import get_xi_params_ds_flux, DiamondFluxCalculation
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess
import numpy as np
from packs.manager.boundary_conditions import BoundaryConditions
from packs.mpfa_methods.weight_interpolation.lpew import LpewWeight, get_lpew2_weights

def get_permeability(faces):
    permeability = np.zeros((len(faces), 2, 2))
    permeability[:, 1, 1] = 1
    permeability[:, 0, 0] = 1
    
    return {'permeability': permeability}

def get_edges_centroids(nodes_centroids, nodes_of_edges):
    resp = (nodes_centroids[nodes_of_edges[: 1]] + nodes_centroids[nodes_of_edges[: 0]])/2
    return resp

def get_dirichlet_edges(edges, edges_centroids):
    
    k1 = 100
    k2 = 1
    
    delta = 1e-10
    edges1 = edges[edges_centroids[:, 0] <= 0 + delta] 
    edges2 = edges[edges_centroids[:, 0] >= edges_centroids[:, 0].max() - delta] 
    
    values1 = np.repeat(k1, len(edges1))
    values2 = np.repeat(k2, len(edges2))
    
    dr_edges = np.concatenate([edges1, edges2])
    values_dirichlet = np.concatenate([values1, values2])
    
    return dr_edges, values_dirichlet

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



def test_xi_params_ds():
    filenames = get_filenames()
    initialize_mesh(**filenames)
    mesh_properties: MeshProperty = load_mesh_properties(filenames['mesh_properties_name'])    
    get_lpew2_weights(mesh_properties)
    get_xi_params_ds_flux(mesh_properties)
    
    # edges_centroids = get_edges_centroids(mesh_properties.nodes_centroids, mesh_properties.nodes_of_edges)
    bc = BoundaryConditions()
    dr_edges, values_dirichlet = get_dirichlet_edges(mesh_properties.edges, mesh_properties.edges_centroids)
    neumann_edges, values_neumann = get_neumann_edges(mesh_properties.edges, mesh_properties.edges_centroids)
    
    bc.set_boundary('dirichlet_edges', dr_edges, values_dirichlet)
    bc.set_boundary('neumann_edges', neumann_edges, values_neumann)
    
    lpew2 = LpewWeight()
    neumann_weights = lpew2.create_lpew2_neumann_weights(
        edges_neumann=neumann_edges,
        neumann_values=values_neumann,
        **mesh_properties.get_all_data()        
    )
    
    mesh_properties.insert_or_update_data(neumann_weights)
    
    diamondflux = DiamondFluxCalculation()
    diamondflux.mount_problem(mesh_properties.edges_dim, bc, **mesh_properties.get_all_data())
    
    
    
    import pdb; pdb.set_trace()
    
    
    
    

    



