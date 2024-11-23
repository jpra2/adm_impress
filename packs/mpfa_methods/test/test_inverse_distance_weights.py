from packs import defpaths
import os
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists, mesh_verify, create_properties
from packs.manager.meshmanager import MeshProperty
import numpy as np



def test_1():
       
    su_meshs = defpaths.load_su_mesh_paths()
    weights_interpolation = [gls]
    pr_names = ['problem1']
    
    for mesh in su_meshs:
        only_mesh = os.path.split(mesh)[1]
        only_mesh, ext = os.path.splitext(only_mesh)
        mesh_properties_name = only_mesh
        mesh_properties: MeshProperty = create_properties_if_not_exists(mesh, mesh_properties_name)
        
        keys_prop = list(mesh_properties.keys())
        pr_output_tag = pr_name + '_' + 'output'
        dtype_output = [('log10l1norm', np.float64), ('log10l2norm', np.float64), ('len_faces', int)]
        
        if pr_output_tag not in keys_prop:
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