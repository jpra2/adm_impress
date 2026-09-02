from packs.mpfa_methods.mesh_preprocess import preprocess_mesh
from packs import defpaths
from packs.manager import BoundaryConditions
from packs.examples.benchmarks_monophasic.cross.test_cross_2 import set_weights_nodes, set_fine_transmissibility_v2, set_fine_transmissibility_without_bc_v2
from packs.manager.meshmanager2 import MeshProperty

import os
import numpy as np
import anndata as ad
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve


def get_R(theta):
    
    R = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])
    
    return R
 
    
def set_permeability(fp:MeshProperty):
    
    k1 = np.array([
        [1000, 0],
        [0, 1]
    ])
    
    k2 = np.array([
            [100, 0],
            [0, 10]
    ])*1
    
    theta = np.deg2rad(45)
    R = get_R(theta)
    
    k_rot = R @ k1 @ R.T
    # k_rot = k1*100
    # k2_rot = R @ k2 @ R.T
    
    nfaces = len(fp.faces)
    perm = np.zeros((nfaces, 2, 2))
    perm[:] = k1
    
    xcentroids = fp.faces_centroids[:,0]
    
    test1 = xcentroids>3
    
    perm[test1] = k_rot
    # perm[test1] = k2
    
    
    
    
    
    fp.insert_or_update_data({'permeability': perm})
    fp.export_data()



def set_boundary_conditions(fine_properties: MeshProperty):
    bc = BoundaryConditions()

    nodes_centroids = fine_properties['nodes_centroids']
    faces = fine_properties['faces']
    faces_centroids = fine_properties['faces_centroids']

    xmin, ymin = nodes_centroids.min(axis=0)
    xmax, ymax = nodes_centroids.max(axis=0)

    # c_p1 = np.array([xmin, ymin])
    # c_p0 = np.array([xmax, ymax])

    c_p1 = np.array([xmin, ymax])
    c_p0 = np.array([xmax, ymin])

    dists = np.linalg.norm(faces_centroids - c_p1, axis=1)
    face_p1 = faces[dists <= dists.min()][0]
    dists[:] = np.linalg.norm(faces_centroids - c_p0, axis=1)
    face_p0 = faces[dists <= dists.min()][0]
    
    faces_pressure = np.array([face_p0, face_p1])
    pressure_presc = np.array([101, 1.0])

    # faces_neumann = np.array([face_p1])
    # neummann_presc_faces = np.array([2.0])

    # areas = fine_properties['areas']
    # edges_dim = fine_properties.edges_dim
    # area_face_p1 = areas[face_p1]

    # import pdb; pdb.set_trace()

    bc.set_boundary('dirichlet_volumes', faces_pressure, pressure_presc)

    # bc.set_boundary('neumann_volumes', faces_neumann, neummann_presc_faces)

    walls_edges = fine_properties['edges'][fine_properties['bool_boundary_edges']]

    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    fine_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    # bc.set_boundary('water_saturation_volumes', np.array([face_p1]), np.array([1.0]))
    # bc.set_boundary('water_saturation_edges', np.array([]), np.array([]))

    bc.set_boundary('injectors', np.array([face_p1]), np.array([True]))
    bc.set_boundary('producers', np.array([face_p0]), np.array([True]))

    bc.update_zero_bcs()

    return bc


def get_properties():
    
    fine_mesh_path = os.path.join('mesh', 'quadrado_estruturado.msh')
    fine_mesh_properties_name = 'quadrado_estruturado' 
    fine_mesh_path_v4 = fine_mesh_path

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name, mesh_name_v4=fine_mesh_path_v4)

    return fine_properties, fine_mesh_path


def test_edge(x, y, fp: MeshProperty, T_withoutbc, debug=False):
    edges_centroids = fp.edges_centroids
    selected_edge = fp.edges[
        ((np.abs(edges_centroids[:, 0] - x)) < 1e-10) &
        ((np.abs(edges_centroids[:, 1] - y)) < 1e-10) 
    ]
    
    nodes_selected_edge = fp.nodes_of_edges[selected_edge[0]]
    faces_of_nodes_edge = fp.faces_of_nodes[nodes_selected_edge]
    kl = fp.adjacencies[selected_edge[0]]
    
    nodes_weight = fp.nodes_weights
    xi_params = fp.xi_params
    
    xi_params_edge = xi_params[selected_edge]
    
    node_weight_A = nodes_weight[nodes_weight['node_id'] == nodes_selected_edge[0]]
    node_weight_B = nodes_weight[nodes_weight['node_id'] == nodes_selected_edge[1]]
    
    fA = node_weight_A['face_id']
    fB = node_weight_B['face_id']
    
    transm_K = T_withoutbc[kl[0]]
    transm_L = T_withoutbc[kl[1]]
    
    if debug == True:
        print(transm_K)
        print('#'*30)
        print(transm_L)
        print('#'*30)
        print(node_weight_A)
        print('#'*30)
        print(node_weight_B)
        print('#'*30)
        print(fA)
        print('#'*30)
        print(fB)
        print('#'*30)
        print(kl)
        import ipdb; ipdb.set_trace()
    

def run5():

    fp, fine_mesh_path = get_properties()
    set_permeability(fp)
    bc = set_boundary_conditions(fp)
    set_weights_nodes(fp, update=True)
    resp = set_fine_transmissibility_v2(fp, bc)
    T_bc = resp['transmissibility']
    b_bc = resp['source']
    solution = spsolve(T_bc, b_bc)
    resp2 = set_fine_transmissibility_without_bc_v2(fp)
    T_complete = resp2['transmissibility_without_bc']
    
    
    
    ## edge na interface
    test_edge(3, 2.5, fp, T_complete, debug=True)
    ## edge na regiao x > media
    test_edge(4, 2.5, fp, T_complete, debug=True)
    ## edge na regiao x < media
    test_edge(2, 2.5, fp, T_complete, debug=True)

    
    
    import ipdb; ipdb.set_trace()
    
    
    
    
    