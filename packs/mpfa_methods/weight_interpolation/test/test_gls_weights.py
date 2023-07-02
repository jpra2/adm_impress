import os
from packs import defpaths
from packs.utils import calculate_face_properties
import numpy as np
from packs.manager.meshmanager import MeshProperty, create_initial_mesh_properties, load_mesh_properties
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights, mount_weight_matrix
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.manager.boundary_conditions import BoundaryConditions
from scipy.sparse.linalg import spsolve
from packs.manager.mesh_data import MeshData

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

def test_neumann_weights():
    problem_name = 'mpfad_test_neumman_weights'
    mesh_test_name = defpaths.mpfad_test_mesh
    mesh_properties_name = defpaths.mpfad_mesh_properties_neumann_name
    mesh_verify(mesh_test_name)
    mesh_properties = create_properties_if_not_exists(mesh_test_name, mesh_properties_name)
    data = mesh_properties.get_all_data()
    data.update({
        'nodes_to_calculate': data['nodes']
    })
    
    lsds = LsdsFluxCalculation()
    
    mesh_properties.update_data(
            lsds.preprocess(**mesh_properties.get_all_data())
        )
    
    mesh_properties.remove_data(['nodes_to_calculate'])
    
    bc = BoundaryConditions()
    bc.insert_name(problem_name)
    
    xmin, ymin = mesh_properties.faces_centroids[:, 0:2].min(axis=0)
    xmax, ymax = mesh_properties.faces_centroids[:, 0:2].max(axis=0)
    
    mesh_delta = 1e-13
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
    
    nodes_bc = nodes_xmax.astype(np.uint64)
    # nodes_bc = np.concatenate([nodes_xmin, nodes_xmax]).astype(np.uint64)
    centroids_nodes_bc = mesh_properties.nodes_centroids[nodes_bc, 0:2]
    pressures_bc = linear_function(centroids_nodes_bc[:, 0])
    bc.set_boundary('nodes_pressures', nodes_bc, pressures_bc)
    
    nodes_of_edges = mesh_properties.nodes_of_edges
    edges_centroids = np.mean(mesh_properties.nodes_centroids[nodes_of_edges], axis=1)
    
    edges_x_min = mesh_properties.edges[
        edges_centroids[:,0] < xmin + mesh_delta
    ]
    
    values_neumann_edges_x_min = np.repeat(10.0, len(edges_x_min))
    
    bc.set_boundary('neumann_edges', edges_x_min, values_neumann_edges_x_min)
    # bc.set_boundary('neumann_edges', np.array([]), np.array([]))
    
    mesh_properties.insert_data({
        'neumann_edges': bc.neumann_edges['id'],
        'neumann_edges_value': bc.neumann_edges['value'],
        'nodes_to_calculate': mesh_properties.nodes[:]
    })
    # import pdb; pdb.set_trace()
    
    
    mesh_properties.insert_data(
        get_gls_nodes_weights(**mesh_properties.get_all_data())
    )
    mesh_properties.insert_data(
        lsds.get_all_edges_flux_params(**mesh_properties.get_all_data())
    )
    
    
    resp = lsds.mount_problem(
        mesh_properties.nodes_weights,
        mesh_properties.xi_params,
        mesh_properties.faces,
        mesh_properties.edges,
        mesh_properties.bool_boundary_edges,
        mesh_properties.adjacencies,
        bc,
        mesh_properties.nodes_of_edges,
        mesh_properties.neumann_nodes_weights
    )
    
    pressure = spsolve(resp['transmissibility'], resp['source'])
    edges_flux = lsds.get_edges_flux(
        mesh_properties.xi_params,
        mesh_properties.nodes_weights,
        mesh_properties.nodes_of_edges,
        pressure,
        mesh_properties.adjacencies,
        bc,
        mesh_properties.neumann_nodes_weights    
    )
    
    faces_flux = lsds.get_faces_flux(
        edges_flux,
        mesh_properties.adjacencies,
        mesh_properties.bool_boundary_edges
    )
    
    _mesh_path = os.path.join(defpaths.mesh, defpaths.mpfad_test_mesh)
    mesh_data = MeshData(mesh_path=_mesh_path)
    mesh_data.create_tag('pressure')
    mesh_data.insert_tag_data('pressure', pressure, 'faces', mesh_properties.faces)
    mesh_data.create_tag('faces_flux')
    mesh_data.insert_tag_data('faces_flux', faces_flux, 'faces', mesh_properties.faces)
    mesh_data.export_all_elements_type_to_vtk('testglsneumann', 'faces')
    
    # print(mesh_properties.keys())
    # print(bc.keys())
    
    
    
    import pdb; pdb.set_trace()
    
    
    
    

    
    