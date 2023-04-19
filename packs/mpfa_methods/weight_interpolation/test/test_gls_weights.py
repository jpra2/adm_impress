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
        'unitary_normal_edges': unitary_normal_edges,
        'edges_dim': norma
    })

    permeability = np.zeros((len(mesh_properties.faces), 2, 2))
    permeability[:, [0, 1], [0, 1]] = 1
    mesh_properties.insert_data({'permeability': permeability})

    mesh_properties.export_data()

def mount_weight_matrix(nodes_weights):
    n_faces = nodes_weights['face_id'].max() + 1
    n_nodes = nodes_weights['node_id'].max() + 1
    mweight = np.zeros((n_nodes, n_faces))

    lines = np.array([], dtype=int)
    cols = lines.copy()
    data = np.array([], dtype=np.float64)

    for node in np.unique(nodes_weights['node_id']):
        faces = nodes_weights['face_id'][nodes_weights['node_id'] == node]
        weights = nodes_weights['weight'][nodes_weights['node_id'] == node]
        lines = np.append(lines, np.repeat(node, faces.shape[0]))
        cols = np.append(cols, faces)
        data = np.append(data, weights)
    
    mweight[lines, cols] = data
    
    return mweight

def constant_verify(weights_matrix, mesh_properties):
    const_faces_solution = constant_function(mesh_properties.faces_centroids[:, 0])
    nodes_const_solution_interpolated = np.dot(weights_matrix, const_faces_solution)
    n_nodes = weights_matrix.shape[0]
    assert nodes_const_solution_interpolated.sum() == n_nodes

def linear_function(x):
    return x

def linear_verify(weights_matrix, mesh_properties):
    tolerance = 1e-12
    linear_faces_solution = linear_function(mesh_properties.faces_centroids[:, 0])
    nodes_linear_solution_interpolated = np.dot(weights_matrix, linear_faces_solution)
    x_min_faces = mesh_properties.faces_centroids[:,0].min()
    nodes_x_min = mesh_properties.nodes[mesh_properties.nodes_centroids[:,0] < x_min_faces]
    result_nodes_x_min = linear_function(mesh_properties.nodes_centroids[:,0][nodes_x_min])
    nodes_linear_solution_interpolated[nodes_x_min] = result_nodes_x_min
    x_max_faces = mesh_properties.faces_centroids[:,0].max()
    nodes_x_max = mesh_properties.nodes[mesh_properties.nodes_centroids[:,0] > x_max_faces]
    result_nodes_x_max = linear_function(mesh_properties.nodes_centroids[:,0][nodes_x_max])
    nodes_linear_solution_interpolated[nodes_x_max] = result_nodes_x_max
    nodes_linear_solution = linear_function(mesh_properties.nodes_centroids[:, 0])
    erro = np.absolute(nodes_linear_solution - nodes_linear_solution_interpolated) > tolerance
    
    assert erro.sum() == 0

def quadratic_function(x):
    return np.power(x, 2)

def quadratic_verify(weights_matrix, mesh_properties):
    tolerance = 0.000001
    q_faces_solution = quadratic_function(mesh_properties.faces_centroids[:, 0])
    nodes_q_solution_interpolated = np.dot(weights_matrix, q_faces_solution)
    nodes_q_solution = quadratic_function(mesh_properties.nodes_centroids[:, 0])
    erro = np.absolute(nodes_q_solution - nodes_q_solution_interpolated) > tolerance
    # TODO terminar

def create_properties_if_not_exists(mesh_name, mesh_properties_name):
    mesh_properties = MeshProperty()
    mesh_properties.insert_mesh_name(mesh_properties_name)
    if mesh_properties.exists():
        mesh_properties.load_data()
        return mesh_properties
    else:
        create_properties(mesh_name, mesh_properties_name)
        return load_mesh_properties(mesh_properties_name)
    
def test_weights():
    ## verify if exists the mesh test
    mesh_test_name = defpaths.mpfad_test_mesh
    mesh_properties_name = defpaths.mpfad_mesh_properties_name
    mesh_verify(mesh_test_name)
    mesh_properties = create_properties_if_not_exists(mesh_test_name, mesh_properties_name)
    data = mesh_properties.get_all_data()
    data.update({
        'nodes_to_calculate': data['nodes']
    })

    nodes_weights = get_gls_nodes_weights(**data)
    weights_matrix = mount_weight_matrix(nodes_weights)

    constant_verify(weights_matrix, mesh_properties)
    linear_verify(weights_matrix, mesh_properties)
    quadratic_verify(weights_matrix, mesh_properties)