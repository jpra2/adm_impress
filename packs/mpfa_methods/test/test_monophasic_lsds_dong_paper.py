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

def exact_solution_p1(centroids):
    x = centroids[:, 0]
    y = centroids[:, 1]
    term1 = np.sin((1 - x)*(1 - y))/np.sin(1)
    term2 = np.power(1-x, 3)*np.power(1-y, 2)
    return 0.5*(term1 + term2)

def get_permeability_p1(n_elements):
    permeability = np.zeros((n_elements, 2, 2))
    K = np.array([
        [1.5, 0.5],
        [0.5, 1.5]
    ])
    permeability[:,:] = K
    return {'permeability': permeability}

def run(pr_name, mesh_type, ns, n):
    
    mesh_test_name = defpaths.load_mpfad_meshtest_by_type_and_number(mesh_type, ns[n])
    mesh_properties_name = pr_name + '_' + mesh_type + '_' + str(ns[n])
    mesh_properties: MeshProperty = create_properties_if_not_exists(mesh_test_name, mesh_properties_name)
    pressure_tag = 'pressure'
    keys_prop = list(mesh_properties.keys())
    
    if pressure_tag not in keys_prop:
        # define nodes to calculate_weights
        mesh_properties.insert_data({'nodes_to_calculate': mesh_properties.nodes.copy()})
        ## create weights and xi params for flux calculation
        lsds = LsdsFluxCalculation()
        mesh_properties.update_data(
            lsds.preprocess(**mesh_properties.get_all_data())
        )
        
        mesh_properties.insert_data({
            'nodes_weights': get_gls_nodes_weights(**mesh_properties.get_all_data()),
            'xi_params': lsds.get_all_edges_flux_params(**mesh_properties.get_all_data())
        })
        
        mesh_properties.remove_data(['nodes_to_calculate'])
        mesh_properties.update_data(
            get_permeability_p1(len(mesh_properties.faces))
        )
        
        #define boundary conditions    
        problem_name = pr_name
        bc = BoundaryConditions()
        bc.insert_name(problem_name)
        
        xmin, ymin = mesh_properties.faces_centroids[:, 0:2].min(axis=0)
        xmax, ymax = mesh_properties.faces_centroids[:, 0:2].max(axis=0)
        
        mesh_delta = 1e-10
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
        pressures_bc = exact_solution_p1(centroids_nodes_bc)
        
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
        
        mesh_properties.insert_data({
            pressure_tag: pressure,
            'edges_flux': edges_flux,
            'faces_flux': faces_flux,
            'p_exact': p_exact,
            'error': error
        })
        mesh_properties.export_data()
    
    l1_norm = mesh_properties.error.max()
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

    
        

def testp1_by_meshtype(mesh_type, ns):
    pr_name = 'problem1'
    
    all_l1_error = []
    all_l2_error = []
    all_n_faces = []
    
    for n in range(len(ns)):
        l1_norm, l2_norm, n_faces = run(pr_name, mesh_type, ns, n)
        import pdb; pdb.set_trace()

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
    
    mesh_types_dict = {
        'mesh1': [8, 32, 64, 128],
        'mesh2': [0, 1, 2, 3, 4, 5, 6, 7],
        'mesh5': [12, 24, 48, 96, 192, 384],
        'mesh6': [1, 2, 3, 4]  
    }   
    
    fig1, ax1 = plt.subplots(1)
    fig2, ax2 = plt.subplots(1)
    
    mesh_types = list(mesh_types_dict.keys())
    for mesh_type in mesh_types:
        resp = testp1_by_meshtype(mesh_type, mesh_types_dict[mesh_type])
        ax1.plot(resp['n_faces'], np.log10(resp['l1_norm']), label=mesh_type)
        ax2.plot(resp['n_faces'], np.log10(resp['l2_norm']), label=mesh_type)
    
    ax1.set_xlabel('N faces')
    ax1.set_ylabel('Log10 L1 norm')
    ax2.set_xlabel('N faces')
    ax2.set_ylabel('Log10 L2 norm')
    ax1.set_title('L1 Norm')
    ax2.set_title('L2 Norm')
    ax1.legend()
    ax2.legend()
    
    fig1.savefig('L1 Norm.svg', format='svg')
    fig2.savefig('L1 Norm.svg', format='svg')

    





