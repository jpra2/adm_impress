from packs import defpaths
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists, mesh_verify, create_properties
from packs.manager.meshmanager import MeshProperty
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.manager.boundary_conditions import BoundaryConditions
import numpy as np
from scipy.sparse.linalg import spsolve
import os
from packs.manager.mesh_data import MeshData
import matplotlib.pyplot as plt
from packs.utils import calculate_face_properties

all_pr_names = ['problem1', 'problem2', 'problem3']

def exact_solution_p1(centroids):
    x = centroids[:, 0]
    y = centroids[:, 1]
    term1 = np.sin((1 - x)*(1 - y))/np.sin(1)
    term2 = np.power(1-x, 3)*np.power(1-y, 2)
    return 0.5*(term1 + term2)

def exact_solution_p2(centroids):
    x = centroids[:, 0]
    y = centroids[:, 1]
    resp = 16*x*(1-y)*y*(1-y)
    return resp
    
def exact_solution_p3(centroids):
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

def run(pr_name, mesh_type, ns, n):
    global all_pr_names
    
    
    mesh_test_name = defpaths.load_mpfad_meshtest_by_type_and_number(mesh_type, ns[n])
    mesh_properties_name = pr_name + '_' + mesh_type + '_' + str(ns[n])
    mesh_properties: MeshProperty = create_properties_if_not_exists(mesh_test_name, mesh_properties_name)
    pressure_tag = 'pressure'
    keys_prop = list(mesh_properties.keys())
    
    
    if pressure_tag not in keys_prop:
        # define nodes to calculate_weights
        mesh_properties.insert_data({'nodes_to_calculate': mesh_properties.nodes.copy()})
        ## create weights and xi params for flux calculation
        centroids_nodes = mesh_properties.nodes_centroids
        z_centroids = np.zeros((len(centroids_nodes), 1))
        centroids_nodes = np.hstack([centroids_nodes, z_centroids])
        nodes_of_faces = mesh_properties.nodes_of_faces
        cnodes_faces = centroids_nodes[nodes_of_faces]
        n_faces = len(mesh_properties.faces)        
        areas = np.zeros(n_faces)
        for i in range(n_faces):
            areas[i] = calculate_face_properties.polygon_area(cnodes_faces[i])
        mesh_properties.insert_data({'areas': areas})
        
        lsds = LsdsFluxCalculation()
        mesh_properties.update_data(
            lsds.preprocess(**mesh_properties.get_all_data())
        )
        
        mesh_properties.insert_data({
            'nodes_weights': get_gls_nodes_weights(**mesh_properties.get_all_data()),
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
        
        p_exact = exact_solution_p1(mesh_properties.faces_centroids[:, 0:2])
        
        error = np.absolute((pressure - p_exact))
        
        error2 = np.power(error, 2)
        
        mesh_properties.insert_data({
            pressure_tag: pressure,
            'edges_flux': edges_flux,
            'faces_flux': faces_flux,
            'p_exact': p_exact,
            'error': error,
            'error2': error2
        })
        mesh_properties.export_data()
    
    l1_norm = mesh_properties.error.max()
    # Eu = np.sqrt(np.dot(mesh_properties.areas, mesh_properties.error2))
    # l1_norm = Eu
    l2_norm = np.linalg.norm(mesh_properties.error)
    
    
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
    
    return l1_norm, l2_norm, len(mesh_properties.faces)
   
        

def testp1_by_meshtype(mesh_type, ns, pr_name):
    
    all_l1_error = []
    all_l2_error = []
    all_n_faces = []
    
    for n in range(len(ns)):
        l1_norm, l2_norm, n_faces = run(pr_name, mesh_type, ns, n)
        # import pdb; pdb.set_trace()

        all_l1_error.append(l1_norm)
        all_l2_error.append(l2_norm)
        all_n_faces.append(n_faces)
    
    return {
        'l1_norm': all_l1_error,
        'l2_norm': all_l2_error,
        'n_faces': all_n_faces
    }

        
        
        
    
def plot_errors():
    # 'mesh1': [8, 32, 64, 128]
    global all_pr_names
    pr_name = all_pr_names[2]
    
    mesh_types_dict = {
        # 'mesh1': [8, 32, 64, 128],
        # 'mesh2': [0, 1, 2, 3, 4, 5, 6, 7],
        'mesh5': [12, 24, 48, 96, 192, 384],
        # 'mesh6': [1, 2, 3, 4]   
    }
    
    fig1, ax1 = plt.subplots(1)
    fig2, ax2 = plt.subplots(1)
    
    all_resps = []
    
    mesh_types = list(mesh_types_dict.keys())
    for mesh_type in mesh_types:
        resp = testp1_by_meshtype(mesh_type, mesh_types_dict[mesh_type], pr_name)
        all_resps.append(resp.update({'mesh_type': mesh_type}))
        all_resps.append(resp.update({'mesh_type': mesh_type}))
        ax1.plot(np.log10(resp['n_faces']), np.log10(resp['l1_norm']), label=mesh_type)
        ax2.plot(np.log10(resp['n_faces']), np.log10(resp['l2_norm']), label=mesh_type)
    
    ax1.set_xlabel('N faces')
    ax1.set_ylabel('Log10 Linf norm')
    ax2.set_xlabel('Log 10 N faces')
    ax2.set_ylabel('Log10 L2 norm')
    ax1.set_title(pr_name)
    ax2.set_title(pr_name)
    ax1.legend()
    ax2.legend()
    
    fig1.savefig(os.path.join('results', 'L1_Norm.svg'), format='svg')
    fig2.savefig(os.path.join('results', 'L2_Norm.svg'), format='svg')

    





