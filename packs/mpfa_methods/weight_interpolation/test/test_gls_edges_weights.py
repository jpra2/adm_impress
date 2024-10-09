from packs.multiscale.unstructured.test.create_primal_test import get_fine_mesh_path_and_mesh_properties_name_for_test
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs.mpfa_methods.weight_interpolation.edge_gls_weight import EdgeGlsWeight
from packs.manager import MeshProperty, MeshData, BoundaryConditions
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights

import numpy as np
from scipy.sparse.linalg import spsolve
import scipy.sparse as sp
import matplotlib.pyplot as plt

def define_boundary_conditions1(fine_mesh_properties: MeshProperty):
    bc_name = 'test_edges_weights'
    bc = BoundaryConditions(name=bc_name)

    dirichlet_edges_1 = fine_mesh_properties['Inflow']
    dirichlet_edges_0 = fine_mesh_properties['Outflow']
    walls_edges = fine_mesh_properties['Walls']

    p0 = 100.0
    pL = 1.0

    all_dirichlet_edges = np.concatenate([dirichlet_edges_1, dirichlet_edges_0])
    values_edges_dirichlet = np.concatenate([
        np.repeat(p0, dirichlet_edges_1.shape[0]),
        np.repeat(pL, dirichlet_edges_0.shape[0])
    ])
    bc.set_boundary('dirichlet_edges', all_dirichlet_edges, values_edges_dirichlet)

    dirichlet_nodes1 = np.unique(
        fine_mesh_properties['nodes_of_edges'][
            dirichlet_edges_1
        ])
    
    dirichlet_nodes0 = np.unique(
        fine_mesh_properties['nodes_of_edges'][
            dirichlet_edges_0
        ])
    
    bc_nodes = np.concatenate([dirichlet_nodes1, dirichlet_nodes0])
    nodes_values = np.concatenate([
        np.repeat(p0, dirichlet_nodes1.shape[0]),
        np.repeat(pL, dirichlet_nodes0.shape[0])
    ])

    bc.set_boundary('dirichlet_nodes', bc_nodes, nodes_values)

    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    return bc

def define_boundary_conditions2(fine_mesh_properties: MeshProperty):
    bc_name = 'test_edges_weights'
    bc = BoundaryConditions(name=bc_name)

    neumann_edges_1 = fine_mesh_properties['Inflow']
    dirichlet_edges_0 = fine_mesh_properties['Outflow']
    walls_edges = fine_mesh_properties['Walls']

    q0 = 1.0
    pL = 1.0

    all_dirichlet_edges = dirichlet_edges_0
    values_edges_dirichlet = np.repeat(pL, dirichlet_edges_0.shape[0])

    bc.set_boundary('dirichlet_edges', all_dirichlet_edges, values_edges_dirichlet)
    
    dirichlet_nodes0 = np.unique(
        fine_mesh_properties['nodes_of_edges'][
            dirichlet_edges_0
        ])
    
    bc_nodes = dirichlet_nodes0
    nodes_values = np.repeat(pL, dirichlet_nodes0.shape[0])

    bc.set_boundary('dirichlet_nodes', bc_nodes, nodes_values)

    neumann_edges = np.concatenate([
        neumann_edges_1,
        walls_edges
    ])

    neumann_edges_values = np.concatenate([
        np.repeat(q0, neumann_edges_1.shape[0]),
        np.repeat(0.0, walls_edges.shape[0])
    ])
    
    bc.set_boundary('neumann_edges', neumann_edges, neumann_edges_values)

    return bc

def update_permeability1(fine_mesh_properties: MeshProperty):
    n_faces = fine_mesh_properties['faces'].shape[0]
    permeability = np.zeros((n_faces, 2, 2))
    permeability[:, 0, 0] = 3
    permeability[:, 1, 1] = 3
    fine_mesh_properties.insert_or_update_data({
        'permeability': permeability
    })


    


def run():
    fine_mesh_path, fine_mesh_properties_name, fine_mesh_path_v4 = get_fine_mesh_path_and_mesh_properties_name_for_test()
    
    fine_mesh_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name, mesh_name_v4=fine_mesh_path_v4)
    update_permeability1(fine_mesh_properties)

    bc = define_boundary_conditions2(fine_mesh_properties)
    fine_mesh_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    lsds = LsdsFluxCalculation()

    if not fine_mesh_properties.verify_name_in_data_names('nodes_weights'):
        fine_mesh_properties.insert_or_update_data(
            {'nodes_to_calculate': fine_mesh_properties['nodes']}
        )
        fine_mesh_properties.insert_or_update_data(
            get_gls_nodes_weights(**fine_mesh_properties.get_all_data())
        )
        fine_mesh_properties.export_data()
    
    if not fine_mesh_properties.verify_name_in_data_names('xi_params'):
        fine_mesh_properties.insert_or_update_data(
            lsds.get_all_edges_flux_params(**fine_mesh_properties.get_all_data())
        )
        fine_mesh_properties.export_data()

    edge_weights = EdgeGlsWeight()
    edges_weight, faces_weight, weight, all_neumann_weights = edge_weights.create_weights(**fine_mesh_properties.get_all_data())

    

    resp = lsds.mount_transmissibility_matrix(
        bc,
        **fine_mesh_properties.get_all_data()
    )

    pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])
    x_faces_centroids = fine_mesh_properties['faces_centroids'][:, 0]

    n_edges = fine_mesh_properties['edges'].shape[0]
    n_faces = fine_mesh_properties['faces'].shape[0]
    edges_matrix_weight = sp.csc_matrix((weight, (edges_weight, faces_weight)), shape=(n_edges, n_faces))

    edges_pressure = edges_matrix_weight.dot(pressure)
    edges_pressure[fine_mesh_properties['neumann_edges']] += all_neumann_weights
    edges_pressure[bc['dirichlet_edges']['id']] = bc['dirichlet_edges']['value']
    x_edges_centroids = fine_mesh_properties.edges_centroids[:, 0]


    # import pdb; pdb.set_trace()

    # test = np.nonzero(edges_pressure)[0]
    # edges_pressure = edges_pressure[test]
    # x_edges_centroids = x_edges_centroids[test]

    # import pdb; pdb.set_trace()

    nodes_pressures = lsds.get_nodes_pressures(
        boundary_conditions=bc,
        faces_pressures=pressure,
        **fine_mesh_properties.get_all_data()
    )
    x_nodes_centroids = fine_mesh_properties['nodes_centroids'][:, 0]

    plt.plot(x_faces_centroids, pressure, 'bo', label='faces', markersize=5)
    plt.plot(x_nodes_centroids, nodes_pressures, 'go', label='nodes', markersize=3)
    plt.plot(x_edges_centroids, edges_pressure, 'ro', label='edges', markersize=1)
    plt.legend()
    plt.savefig('test_edges_all.svg')








    import pdb; pdb.set_trace()
    