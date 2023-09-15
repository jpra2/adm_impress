from packs import defpaths, defnames
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists, mesh_verify, create_properties
from packs.manager.meshmanager import MeshProperty
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights, mount_sparse_weight_matrix
from packs.mpfa_methods.weight_interpolation.lpew import get_lpew2_weights
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.manager.boundary_conditions import BoundaryConditions
import numpy as np
from scipy.sparse.linalg import spsolve
import scipy.sparse as sp
import os
from packs.manager.mesh_data import MeshData
import matplotlib.pyplot as plt
from packs.utils import calculate_face_properties

all_pr_names = ['problem1', 'problem2', 'problem3', 'problem4', 'problem5']
all_weight_interpolation_names = ['gls', 'inverse_distance']

def exact_solution_p1(centroids):
    """example 5.2 - Mild anisotropy
        DOI: 10.1002/fld.5031
    """
    x = centroids[:, 0]
    y = centroids[:, 1]
    term1 = np.sin((1 - x)*(1 - y))/np.sin(1)
    term2 = np.power(1-x, 3)*np.power(1-y, 2)
    return 0.5*(term1 + term2)

def exact_solution_p2(centroids):
    """Test 1.1 Mild anisotropy
    paper: Benchmark on Discretization Schemes
           for Anisotropic Diffusion Problems
           on General Grids    
    """
    x = centroids[:, 0]
    y = centroids[:, 1]
    resp = 16*x*(1-y)*y*(1-y)
    return resp
    
def exact_solution_p3(centroids):
    """example 5.3 - Discontinuous diffusion tensor with one discontinuity line
        DOI: 10.1002/fld.5031
    """
    x = centroids[:, 0]
    y = centroids[:, 1]
    
    solution = np.zeros(len(x))
    
    t1 = x <= 0.5
    t2 = x > 0.5
    
    x1 = x[t1]
    y1 = y[t1]
    
    x2 = x[t2]
    y2 = y[t2]
    
    solution[t1] = 1 - 2*np.power(y1, 2) + 4*x1*y1 + 6*x1 + 2*y1
    solution[t2] = -2*np.power(y2, 2) + 1.6*x2*y2 -0.6*x2 + 3.2*y2 + 4.3
    
    return solution
    
def exact_solution_p4(centroids):
    """5.5 Heterogeneous rotating anisotropy
        DOI: 10.1002/fld.5031
    """
    x = centroids[:,0]
    y = centroids[:,1]
    x2 = np.pi*x
    y2 = np.pi*y
    
    resp = np.sin(x2)*np.sin(y2)
    return resp
    
def exact_solution_p5(centroids):
    """5.6 Anisotropy of Righi-Leduc type
        DOI: 10.1002/fld.5031
    """
    
    a = 1
    b = 2
    
    x = centroids[:,0]
    y = centroids[:,1]
    
    term1 = np.exp(b*x/a) - 1
    term2 = np.exp(b/a) - 1
    
    solution = term1/term2
    return solution


def get_permeability_p5(n_elements, centroids):
    a = 1
    b = 2
    permeability = np.zeros((n_elements, 2, 2))
    x = centroids[:,0]
    y = centroids[:,1]
    
    permeability[:, 0, 0] = np.repeat(a, n_elements)
    permeability[:, 1, 1] = np.repeat(a, n_elements)
    permeability[:, 0, 1] = b*y
    permeability[:, 1, 0] = -b*y
    
    return {
        'permeability': permeability
    }
    
def get_permeability_p4(n_elements, centroids):
    x = centroids[:,0]
    y = centroids[:,1]
    alpha = 1e-3
    permeability = np.zeros((n_elements, 2, 2))
    permeability[:, 0, 0] = alpha*np.power(x, 2) + np.power(y, 2)
    permeability[:, 0, 1] = (alpha - 1)*x*y
    permeability[:, 1, 0] = (alpha - 1)*x*y
    permeability[:, 1, 1] = np.power(x, 2) + alpha*np.power(y, 2)
    
    return {
        'permeability': permeability
    }
    
def get_permeability_p3(n_elements, centroids):
    x = centroids[:, 0]
    permeability = np.zeros((n_elements, 2, 2))
    
    k1 = np.array([
        [1, 0],
        [0, 1]
    ])
    
    k2 = np.array([
        [10, 3],
        [3, 1]
    ])
    
    t1 = x <= 0.5
    t2 = x > 0.5
    
    permeability[t1,:] = k1
    permeability[t2,:] = k2
    return {'permeability': permeability}
    
def get_permeability_p2(n_elements, centroids=0):
    return get_permeability_p1(n_elements)
    
def get_permeability_p1(n_elements, centroids=0):
    permeability = np.zeros((n_elements, 2, 2))
    K = np.array([
        [1.5, 0.5],
        [0.5, 1.5]
    ])
    permeability[:,:] = K
    return {'permeability': permeability}

def get_permeability_and_exact_solution_func(pr_name):
    global all_pr_names
    
    if pr_name == all_pr_names[0]:
        get_permeability = get_permeability_p1
        exact_solution = exact_solution_p1
    elif pr_name == all_pr_names[1]:
        get_permeability = get_permeability_p2
        exact_solution = exact_solution_p2
    elif pr_name == all_pr_names[2]:
        get_permeability = get_permeability_p3
        exact_solution = exact_solution_p3
    elif pr_name == all_pr_names[3]:
        get_permeability = get_permeability_p4
        exact_solution = exact_solution_p4
    elif pr_name == all_pr_names[4]:
        get_permeability = get_permeability_p5
        exact_solution = exact_solution_p5
    else:
        raise NameError
    
    return get_permeability, exact_solution

def get_Eu(absolute_error, areas):
    error_2 = np.power(absolute_error, 2)
    resp = areas*error_2
    resp = np.sqrt(resp.sum())
    return resp

def get_Eq(edges_dim, absolute_error, adjacencies, nodes_of_edges, nodes_centroids, faces_centroids, bool_boundary_edges, edges_of_faces):
    ### get d_sigma

    d_k_sigma = np.zeros(adjacencies.shape)
    bool_internal_edges = ~bool_boundary_edges

    centroids_nodes_of_edges = nodes_centroids[nodes_of_edges]
    edges_centroid = np.mean(centroids_nodes_of_edges, axis=1)

    d_k_sigma[:, 0] = np.linalg.norm(
        faces_centroids[adjacencies[:, 0]] - edges_centroid,
        axis=1
    )

    d_k_sigma[:, 1] = np.linalg.norm(
        faces_centroids[adjacencies[:, 1]] - edges_centroid,
        axis=1
    )

    d_k_sigma[bool_boundary_edges, 1] = 0

    d_sigma = d_k_sigma.sum()

    term1 = absolute_error[adjacencies]
    term1[bool_boundary_edges, 1] = 0
    term1 = np.power(term1[:, 0] - term1[:, 1], 2)
    term1 = term1*edges_dim

    faces_term = term1[edges_of_faces]
    faces_term = faces_term.sum(axis=1).sum()
    resp = np.sqrt(faces_term/d_sigma)

    return resp

def nodes_weights_test(mesh_properties: MeshProperty, exact_solution_func, pr_name):
    tag_test = pr_name + defnames.tag_node_weight_test_sufix
    
    exact_solution_faces = exact_solution_func(mesh_properties.faces_centroids[:, 0:2])
    exact_solution_nodes = exact_solution_func(mesh_properties.nodes_centroids[:, 0:2])
    bnodes = mesh_properties.bool_boundary_nodes
    
    nodes_weight_matrix = mount_sparse_weight_matrix(mesh_properties.nodes_weights)
    weighted_nodes_pressures = nodes_weight_matrix.dot(exact_solution_faces)
    
    weighted_nodes_pressures[bnodes] = exact_solution_nodes[bnodes]
    
    error = np.absolute(weighted_nodes_pressures - exact_solution_nodes)
    
    mesh_properties.insert_or_update_data(
        {tag_test: error}
    )

def get_tags(pr_name, sufix=''):
    tags = ['pressure_' + pr_name + sufix,
    'edges_flux_' + pr_name + sufix,
    'faces_flux_' + pr_name + sufix,
    'p_exact_' + pr_name + sufix,
    'error_' + pr_name + sufix,
    'error2_' + pr_name + sufix]

    return tags

def backup_tags_process(mesh_properties: MeshProperty, pr_name: str, sufix: str):
    tags = get_tags(pr_name)
    tags_backup = get_tags(pr_name, sufix=sufix)
    dict_tags = dict(zip(tags, tags_backup))
    mesh_properties.backup_datas(dict_tags)
    mesh_properties.export_data()


def remove_process_tags(mesh_properties: MeshProperty, pr_name: str):
    tags = get_tags(pr_name)
    mesh_properties.remove_data(tags)

def backup_weights(mesh_properties: MeshProperty, sufix_name: str):
    tags = ['nodes_weights', 'neumann_nodes_weights']
    for tag in tags:
        mesh_properties.backup_data(tag, tag+sufix_name)
    
    mesh_properties.export_data()

def calculate_areas(mesh_properties: MeshProperty):
    centroids_nodes = mesh_properties.nodes_centroids
    z_centroids = np.zeros((len(centroids_nodes), 1))
    centroids_nodes = np.hstack([centroids_nodes, z_centroids])
    nodes_of_faces = mesh_properties.nodes_of_faces
    cnodes_faces = centroids_nodes[nodes_of_faces]
    n_faces = len(mesh_properties.faces)

    if not mesh_properties.verify_name_in_data_names('areas'):
        areas = np.zeros(n_faces)
        for i in range(n_faces):
            areas[i] = calculate_face_properties.polygon_area(cnodes_faces[i])
        mesh_properties.insert_data({'areas': areas})
        mesh_properties.export_data()

def calculate_h_dist(mesh_properties: MeshProperty):
    if not mesh_properties.verify_name_in_data_names('h_dist'):
        
        h_dist = calculate_face_properties.create_face_to_edge_distances(
            mesh_properties.faces_centroids,
            mesh_properties.adjacencies,
            mesh_properties.nodes_of_edges,
            mesh_properties.edges,
            mesh_properties.nodes_centroids,
            mesh_properties.bool_boundary_edges
        )
        mesh_properties.insert_data({'h_dist': h_dist})

        bedges = mesh_properties.bool_boundary_edges
        iedges = ~bedges

        m1 = np.mean(h_dist[iedges, 0])
        m2 = np.mean(h_dist[iedges, 1])
        m3 = np.mean(h_dist[bedges, 0])
        m_hdist = np.mean([m1, m2, m3])
        mesh_properties.insert_data({'m_hdist': np.array([m_hdist])})
        mesh_properties.export_data()




   
def run(pr_name, mesh_type, ns, n):    
    
    mesh_test_name = defpaths.load_mpfad_meshtest_by_type_and_number(mesh_type, ns[n])
    # mesh_properties_name = mesh_type + '_' + str(ns[n])
    prop_name = os.path.split(mesh_test_name)[1]
    prop_name = os.path.splitext(prop_name)[0]
    mesh_properties_name = prop_name
    mesh_properties: MeshProperty = create_properties_if_not_exists(mesh_test_name, mesh_properties_name)
    pressure_tag = 'pressure_' + pr_name

    # mesh_properties.remove_data(['nodes_weights', 'neumann_nodes_weights'])
    # mesh_properties.remove_data(get_tags(pr_name))

    lsds = LsdsFluxCalculation()
    
    get_permeability, exact_solution = get_permeability_and_exact_solution_func(pr_name)
    
    calculate_areas(mesh_properties)
    calculate_h_dist(mesh_properties)
    
    if  not mesh_properties.verify_name_in_data_names('nodes_weights'):
        
        # define nodes to calculate_weights
        mesh_properties.insert_or_update_data({'nodes_to_calculate': mesh_properties.nodes.copy()})
        
        ## create weights and xi params for flux calculation
        mesh_properties.insert_or_update_data(
            lsds.preprocess(mesh_properties)
        )
    
        mesh_properties.insert_data(
            get_gls_nodes_weights(**mesh_properties.get_all_data())
        )
        
        # dtype_neumann = [('node_id', np.int), ('nweight', np.float)]
        # zero_neumann = np.zeros(len([]), dtype=dtype_neumann)
        # get_lpew2_weights(mesh_properties)
        # mesh_properties.insert_data({
        #     'neumann_nodes_weights': zero_neumann
        # })

        if not mesh_properties.verify_name_in_data_names('xi_params'):
           
            mesh_properties.insert_data(
                lsds.get_all_edges_flux_params(**mesh_properties.get_all_data())
            )
    
        mesh_properties.remove_data(['nodes_to_calculate'])

        backup_weights(mesh_properties, '_lpew2')
        # backup_weights(mesh_properties, '_gls')
        
        mesh_properties.export_data()  
    
    
    if not mesh_properties.verify_name_in_data_names(pressure_tag):
           
        mesh_properties.update_data(
            get_permeability(len(mesh_properties.faces), mesh_properties.faces_centroids[:, 0:2])
        )
        
        #define boundary conditions    
        problem_name = pr_name
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
        
        nodes_bc = np.unique(np.concatenate([
            nodes_xmin, nodes_xmax, nodes_ymin, nodes_ymax
        ])).astype(np.uint64)
        
        centroids_nodes_bc = mesh_properties.nodes_centroids[nodes_bc, 0:2]
        
        pressures_bc = exact_solution(centroids_nodes_bc)
        
        bc.set_boundary('nodes_pressures', nodes_bc, pressures_bc)

        resp = lsds.mount_problem_v2(
            mesh_properties.nodes_weights,
            mesh_properties.xi_params,
            mesh_properties.faces,
            mesh_properties.nodes,
            mesh_properties.bool_boundary_edges,
            mesh_properties.adjacencies,
            bc,
            mesh_properties.nodes_of_edges,
            mesh_properties.neumann_nodes_weights,
            mesh_properties.edges_of_nodes
        )
        
        pressure = spsolve(resp['transmissibility'].tocsc(), resp['source'])
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
        
        p_exact = exact_solution(mesh_properties.faces_centroids[:, 0:2])
        
        error = np.absolute((pressure - p_exact))
        
        error2 = np.power(error, 2)
        
        mesh_properties.insert_data({
            pressure_tag: pressure,
            'edges_flux_' + pr_name: edges_flux,
            'faces_flux_' + pr_name: faces_flux,
            'p_exact_' + pr_name: p_exact,
            'error_' + pr_name: error,
            'error2_' + pr_name: error2
        })
        backup_tags_process(mesh_properties, pr_name, '_lpew2')
        mesh_properties.export_data()
    
    nodes_weights_test(mesh_properties, exact_solution, pr_name)
    tag_weight_test = pr_name + defnames.tag_node_weight_test_sufix
    
    error_weighted_node_pressure_abs = np.absolute(mesh_properties[tag_weight_test])
    
    l1_weighted_error = error_weighted_node_pressure_abs.max()
    l2_weighted_error = np.linalg.norm(error_weighted_node_pressure_abs)
    
    eq = get_Eq(
        mesh_properties['edges_dim'],
        mesh_properties['error_' + pr_name],
        mesh_properties['adjacencies'],
        mesh_properties['nodes_of_edges'],
        mesh_properties['nodes_centroids'],
        mesh_properties['faces_centroids'],
        mesh_properties.bool_boundary_edges,
        mesh_properties['edges_of_faces']  
    )

    eu = get_Eu(
        mesh_properties['error_' + pr_name],
        mesh_properties['areas']
    )

    l1_norm = mesh_properties['error_' + pr_name].max()
    
    # Eu = np.sqrt(np.dot(mesh_properties.areas, mesh_properties.error2))
    # l1_norm = Eu
    l2_norm = np.linalg.norm(mesh_properties['error_' + pr_name])
    
    
    
    # mesh_path = os.path.join(defpaths.mesh, mesh_test_name)
    # # mesh_path = mesh_test_name
    # mesh_data = MeshData(mesh_path=mesh_path)
    
    # mesh_data.create_tag('pressure')
    # mesh_data.insert_tag_data('pressure', pressure, 'faces', mesh_properties.faces)
    # mesh_data.create_tag('permeability', data_size=4)
    # perm = mesh_properties.permeability.reshape((len(mesh_properties.faces), 4))
    # mesh_data.insert_tag_data('permeability', perm, 'faces', mesh_properties.faces)
    # mesh_data.create_tag('faces_flux')
    # mesh_data.insert_tag_data('faces_flux', faces_flux, 'faces', mesh_properties.faces)
    # mesh_data.create_tag('nodes_pressure_presc')
    # mesh_data.insert_tag_data('nodes_pressure_presc', pressures_bc, 'nodes', nodes_bc)
    # mesh_data.create_tag('pressure_abs_error')
    # mesh_data.insert_tag_data('pressure_abs_error', error, 'faces', mesh_properties.faces)
    # to_export_name = pr_name + mesh_type + str(ns[n])
    # mesh_data.export_all_elements_type_to_vtk(to_export_name + 'nodes', 'nodes')
    # mesh_data.export_all_elements_type_to_vtk(to_export_name + 'faces', 'faces')
    # mesh_data.export_only_the_elements('test_7_nodes_pressure_boundary', 'nodes', nodes_bc)
    
    return l1_norm, l2_norm, len(mesh_properties.faces), mesh_properties.m_hdist[0], eu, eq, l1_weighted_error, l2_weighted_error

def get_tag_prefix(pr_name, weight_interpolation_name):
    return pr_name + '_' + weight_interpolation_name + '_'
   
def run_by_mesh_name_problem_and_weight_interpolation(mesh_properties, pr_name, weight_interpolation_name='gls'):
    global all_pr_names, all_weight_interpolation_names
    assert weight_interpolation_name in all_weight_interpolation_names
        
    tag_prefix = get_tag_prefix(pr_name, weight_interpolation_name)
    # define nodes to calculate_weights
    mesh_properties.insert_data({'nodes_to_calculate': mesh_properties.nodes.copy()})
    ## create weights and xi params for flux calculation
    centroids_nodes = mesh_properties.nodes_centroids
    z_centroids = np.zeros((len(centroids_nodes), 1))
    centroids_nodes = np.hstack([centroids_nodes, z_centroids])
    nodes_of_faces = mesh_properties.nodes_of_faces
    cnodes_faces = centroids_nodes[nodes_of_faces]
    n_faces = len(mesh_properties.faces)
    
    keys_prop = list(mesh_properties.keys())
    if 'areas' not in keys_prop:
        areas = np.zeros(n_faces)
        for i in range(n_faces):
            areas[i] = calculate_face_properties.polygon_area(cnodes_faces[i])
        mesh_properties.insert_data({'areas': areas})
    
    lsds = LsdsFluxCalculation()
    mesh_properties.update_data(
        lsds.preprocess(**mesh_properties.get_all_data())
    )        
    
    if weight_interpolation_name == 'gls':
        mesh_properties.insert_data(
            get_gls_nodes_weights(**mesh_properties.get_all_data())
        )
    elif weight_interpolation_name == 'inverse_distance':
        from packs.mpfa_methods.weight_interpolation.inverse_distance_method import InverseDistanceWeight2D
        mesh_properties.insert_data(
            InverseDistanceWeight2D.get_weights(
                mesh_properties.nodes_centroids,
                mesh_properties.faces_of_nodes,
                mesh_properties.faces_centroids,
                mesh_properties.nodes
            )
        )
    else:
        raise NameError
    
    if 'xi_params' not in keys_prop:    
        mesh_properties.insert_data({
            'xi_params': lsds.get_all_edges_flux_params(**mesh_properties.get_all_data())
        })
    
    mesh_properties.remove_data(['nodes_to_calculate'])
    
    #define boundary conditions    
    problem_name = pr_name
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
    
    nodes_bc = np.unique(np.concatenate([
        nodes_xmin, nodes_xmax, nodes_ymin, nodes_ymax
    ])).astype(np.uint64)
    
    centroids_nodes_bc = mesh_properties.nodes_centroids[nodes_bc, 0:2]
    
    if pr_name == all_pr_names[0]:
        mesh_properties.update_data(
            get_permeability_p1(len(mesh_properties.faces), mesh_properties.faces_centroids[:, 0:2])
        )
        exact_solution = exact_solution_p1
    elif pr_name == all_pr_names[1]:
        mesh_properties.update_data(
            get_permeability_p2(len(mesh_properties.faces), mesh_properties.faces_centroids[:, 0:2])
        )
        exact_solution = exact_solution_p2
    elif pr_name == all_pr_names[2]:
        mesh_properties.update_data(
            get_permeability_p3(len(mesh_properties.faces), mesh_properties.faces_centroids[:, 0:2])
        )
        exact_solution = exact_solution_p3
    elif pr_name == all_pr_names[3]:
        mesh_properties.update_data(
            get_permeability_p4(len(mesh_properties.faces), mesh_properties.faces_centroids[:, 0:2])
        )
        exact_solution = exact_solution_p4
    elif pr_name == all_pr_names[4]:
        mesh_properties.update_data(
            get_permeability_p5(len(mesh_properties.faces), mesh_properties.faces_centroids[:, 0:2])
        )
        exact_solution = exact_solution_p5
    else:
        raise NameError
    
    pressures_bc = exact_solution(centroids_nodes_bc)
    
    bc.set_boundary('nodes_pressures', nodes_bc, pressures_bc)
    
    lsds = LsdsFluxCalculation()
    resp = lsds.mount_problem(
        mesh_properties.nodes_weights,
        mesh_properties.xi_params,
        mesh_properties.faces,
        mesh_properties.edges,
        mesh_properties.bool_boundary_edges,
        mesh_properties.adjacencies,
        bc,
        mesh_properties.nodes_of_edges
    )
    
    pressure = spsolve(resp['transmissibility'], resp['source'])
    edges_flux = lsds.get_edges_flux(
        mesh_properties.xi_params,
        mesh_properties.nodes_weights,
        mesh_properties.nodes_of_edges,
        pressure,
        mesh_properties.adjacencies,
        bc        
    )
    
    faces_flux = lsds.get_faces_flux(
        edges_flux,
        mesh_properties.adjacencies,
        mesh_properties.bool_boundary_edges
    )
    
    p_exact = exact_solution(mesh_properties.faces_centroids[:, 0:2])
    
    error = np.absolute((pressure - p_exact))
    error2 = np.power(error, 2)

    mesh_properties.insert_data({
        tag_prefix + pressure_tag: pressure,
        tag_prefix + 'edges_flux': edges_flux,
        tag_prefix + 'faces_flux': faces_flux,
        tag_prefix + 'p_exact': p_exact,
        tag_prefix + 'error': error,
        tag_prefix + 'error2': error2
    })
    mesh_properties.export_data()
    
    # l1_norm = mesh_properties.error.max()
    # # Eu = np.sqrt(np.dot(mesh_properties.areas, mesh_properties.error2))
    # # l1_norm = Eu
    # l2_norm = np.linalg.norm(mesh_properties.error)
    
    
    # # mesh_path = os.path.join(defpaths.mesh, mesh_test_name)
    # # # mesh_path = mesh_test_name
    # # mesh_data = MeshData(mesh_path=mesh_path)
    
    # # mesh_data.create_tag('pressure')
    # # mesh_data.insert_tag_data('pressure', pressure, 'faces', mesh_properties.faces)
    # # mesh_data.create_tag('permeability', data_size=4)
    # # perm = mesh_properties.permeability.reshape((len(mesh_properties.faces), 4))
    # # mesh_data.insert_tag_data('permeability', perm, 'faces', mesh_properties.faces)
    # # mesh_data.create_tag('faces_flux')
    # # mesh_data.insert_tag_data('faces_flux', faces_flux, 'faces', mesh_properties.faces)
    # # mesh_data.create_tag('nodes_pressure_presc')
    # # mesh_data.insert_tag_data('nodes_pressure_presc', pressures_bc, 'nodes', nodes_bc)
    # # mesh_data.create_tag('pressure_abs_error')
    # # mesh_data.insert_tag_data('pressure_abs_error', error, 'faces', mesh_properties.faces)
    # # to_export_name = pr_name + mesh_type + str(ns[n])
    # # mesh_data.export_all_elements_type_to_vtk(to_export_name + 'nodes', 'nodes')
    # # mesh_data.export_all_elements_type_to_vtk(to_export_name + 'faces', 'faces')
    # # mesh_data.export_only_the_elements('test_7_nodes_pressure_boundary', 'nodes', nodes_bc)
    
    # return l1_norm, l2_norm, len(mesh_properties.faces)
    
    

def testp1_by_meshtype(mesh_type, ns, pr_name):
    
    all_l1_error = []
    all_l2_error = []
    all_n_faces = []
    all_m_hdist = []
    all_eu = []
    all_eq = []
    all_l1_weighted_error = []
    all_l2_weighted_error = []
    
    
    for n in range(len(ns)):
        l1_norm, l2_norm, n_faces, m_hdist, eu, eq, l1_weighted_error, l2_weighted_error = run(pr_name, mesh_type, ns, n)
        # import pdb; pdb.set_trace()

        all_l1_error.append(l1_norm)
        all_l2_error.append(l2_norm)
        all_n_faces.append(n_faces)
        all_m_hdist.append(m_hdist)
        all_eu.append(eu)
        all_eq.append(eq)
        all_l1_weighted_error.append(l1_weighted_error)
        all_l2_weighted_error.append(l2_weighted_error)
    
    return {
        'l1_norm': all_l1_error,
        'l2_norm': all_l2_error,
        'n_faces': all_n_faces, 
        'm_hdist': all_m_hdist,
        'eu': all_eu,
        'eq': all_eq,
        'l1_weighted': all_l1_weighted_error,
        'l2_weighted': all_l2_weighted_error
    }     
        
    
def plot_errors():
    # 'mesh1': [8, 32, 64, 128]
    global all_pr_names
    pr_name = all_pr_names[0]
    
    mesh_types_dict = {
        # 'mesh1': [8, 32, 64, 128],
        # # 'mesh1': [8, 32, 64],
        # 'mesh2': [0, 1, 2, 3, 4, 5, 6, 7],
        # 'mesh5':  [12, 24, 48, 96, 192, 384],
        'mesh6': [1, 2, 3, 4]   
    }
    
    fig1, ax1 = plt.subplots(1)
    fig2, ax2 = plt.subplots(1)
    fig3, ax3 = plt.subplots(1)
    fig4, ax4 = plt.subplots(1)
    fig5, ax5 = plt.subplots(1)
    fig6, ax6 = plt.subplots(1)
    
    
    
    all_resps = []
    
    mesh_types = list(mesh_types_dict.keys())
    for mesh_type in mesh_types:
        resp = testp1_by_meshtype(mesh_type, mesh_types_dict[mesh_type], pr_name)
        all_resps.append(resp.update({'mesh_type': mesh_type}))
        all_resps.append(resp.update({'mesh_type': mesh_type}))
        ax1.plot(np.log10(resp['m_hdist']), np.log10(resp['l1_norm']), label=mesh_type)
        ax2.plot(np.log10(resp['m_hdist']), np.log10(resp['l2_norm']), label=mesh_type)
        ax3.plot(np.log10(resp['m_hdist']), np.log10(resp['eu']), label=mesh_type)
        ax4.plot(np.log10(resp['m_hdist']), np.log10(resp['eq']), label=mesh_type)
        ax5.plot(np.log10(resp['m_hdist']), np.log10(resp['l1_weighted']), label=mesh_type)
        ax6.plot(np.log10(resp['m_hdist']), np.log10(resp['l2_weighted']), label=mesh_type)
        
    
    ax1.set_xlabel('Log10 H medio')
    ax1.set_ylabel('Log10 Linf norm')
    ax1.set_title(pr_name)
    ax1.legend()
    
    ax2.set_xlabel('Log10 H medio')
    ax2.set_ylabel('Log10 L2 norm')
    ax2.set_title(pr_name)
    ax2.legend()

    ax3.set_xlabel('Log10 H medio')
    ax3.set_ylabel('Log10 Eu')
    ax3.set_title(pr_name)
    ax3.legend()

    ax4.set_xlabel('Log10 H medio')
    ax4.set_ylabel('Log10 Eq')
    ax4.set_title(pr_name)
    ax4.legend()
    
    ax5.set_xlabel('Log10 H medio')
    ax5.set_ylabel('Log10 L1 Weight test')
    ax5.set_title(pr_name)
    ax5.legend()
    
    ax6.set_xlabel('Log10 H medio')
    ax6.set_ylabel('Log10 L2 Weight test')
    ax6.set_title(pr_name)
    ax6.legend()
    
    ext = 'svg'
    
    fig1.savefig(os.path.join('results', pr_name + '_' + 'L1_Norm.' + ext), format=ext)
    fig2.savefig(os.path.join('results', pr_name + '_' + 'L2_Norm.' + ext), format=ext)
    fig3.savefig(os.path.join('results', pr_name + '_' + 'Eu.' + ext), format=ext)
    fig4.savefig(os.path.join('results', pr_name + '_' + 'Eq.' + ext), format=ext)
    fig5.savefig(os.path.join('results', pr_name + '_' + 'L1_weight_test.' + ext), format=ext)
    fig6.savefig(os.path.join('results', pr_name + '_' + 'L2_weight_test.' + ext), format=ext)

    





