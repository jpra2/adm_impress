import os
from packs import defpaths
from packs.utils import calculate_face_properties
import numpy as np
from packs.manager.meshmanager import MeshProperty, create_initial_mesh_properties, load_mesh_properties
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights

def mesh_verify(mesh_name):
    mesh_path = os.path.join(defpaths.mesh, mesh_name)

    if os.path.exists(mesh_path):
        pass
    else:
        raise FileExistsError

def constant_function(x):
    return np.repeat(1, x.shape[0])

def create_properties(mesh_name, mesh_properties_name):
    mesh_path = os.path.join(defpaths.mesh, mesh_name)
    mesh_properties: MeshProperty = create_initial_mesh_properties(mesh_path, mesh_properties_name)
    mesh_properties.update_data(
        {
            'nodes_centroids': mesh_properties.nodes_centroids[:, 0:2],
            'faces_centroids': mesh_properties.faces_centroids[:, 0:2]
        }
    )

    datas_to_rename = {
        'faces_adj_by_nodes': 'faces_of_nodes',
        'faces_adj_by_edges': 'adjacencies',
        'nodes_adj_by_nodes': 'nodes_of_nodes',
        'edges_adj_by_nodes': 'edges_of_nodes'
    }

    mesh_properties.rename_data(datas_to_rename)

    norma, unitary_normal_edges = calculate_face_properties.create_unitary_normal_edges_xy_plane(
        mesh_properties.nodes_of_edges,
        mesh_properties.nodes_centroids,
        mesh_properties.adjacencies,
        mesh_properties.faces_centroids,
        mesh_properties.bool_boundary_edges
    )

    mesh_properties.insert_data({
        'unitary_normal_edges': unitary_normal_edges
    })

    permeability = np.zeros((len(mesh_properties.faces), 2, 2))
    permeability[:, [0, 1], [0, 1]] = 1
    mesh_properties.insert_data({'permeability': permeability})

    mesh_properties.export_data()

def mount_weight_matrix(nodes_weights):
    n_faces = nodes_weights['face_id'].max() + 1
    n_nodes = nodes_weights['node_id'].max() + 1
    mweight = np.zeros((n_nodes, n_faces))

    for node in np.unique(nodes_weights['node_id']):
        faces = nodes_weights['face_id'][nodes_weights['node_id'] == node]
        weights = nodes_weights['weight'][nodes_weights['node_id'] == node]
        mweight[node, faces] = weights
    
    return mweight

def constant_verify(weights_matrix, mesh_properties):
    const_faces_solution = constant_function(mesh_properties.faces_centroids[:, 0])
    nodes_const_solution_interpolated = np.dot(weights_matrix, const_faces_solution)
    n_nodes = weights_matrix.shape[0]
    assert nodes_const_solution_interpolated.sum() == n_nodes

def linear_function(x):
    return x

def linear_verify(weights_matrix, mesh_properties):
    tolerance = 0.000001
    linear_faces_solution = linear_function(mesh_properties.faces_centroids[:, 0])
    nodes_linear_solution_interpolated = np.dot(weights_matrix, linear_faces_solution)
    nodes_linear_solution = linear_function(mesh_properties.nodes_centroids[:, 0])
    erro = np.absolute(nodes_linear_solution - nodes_linear_solution_interpolated) > tolerance
    # TODO terminar

def quadratic_function(x):
    return np.power(x, 2)

def quadratic_verify(weights_matrix, mesh_properties):
    tolerance = 0.000001
    q_faces_solution = quadratic_function(mesh_properties.faces_centroids[:, 0])
    nodes_q_solution_interpolated = np.dot(weights_matrix, q_faces_solution)
    nodes_q_solution = quadratic_function(mesh_properties.nodes_centroids[:, 0])
    erro = np.absolute(nodes_q_solution - nodes_q_solution_interpolated) > tolerance
    import pdb; pdb.set_trace()
    # TODO terminar


    

def test_weights():
    ## verify if exists the mesh test
    mesh_test_name_export = '2d_unstructured.msh'
    mesh_properties_name = 'gls_test_weights'
    mesh_verify(mesh_test_name_export)
    create_properties(mesh_test_name_export, mesh_properties_name)
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    data = mesh_properties.get_all_data()
    data.update({
        'nodes_to_calculate': data['nodes']
    })

    nodes_weights = get_gls_nodes_weights(**data)

    weights_matrix = mount_weight_matrix(nodes_weights)

    constant_verify(weights_matrix, mesh_properties)
    linear_verify(weights_matrix, mesh_properties)
    quadratic_verify(weights_matrix, mesh_properties)
    
    








    import pdb; pdb.set_trace()