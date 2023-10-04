from packs.mpfa_methods.weight_interpolation.lpew import LpewWeight, get_lpew2_weights
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists
from packs import defpaths, defnames
from packs.manager.meshmanager import MeshProperty, load_mesh_properties
from packs.utils import calculate_face_properties
import os
from files_to_test.lpew2_test_weights.info_test1 import get_filenames
import scipy.io as sio
import numpy as np

def initialize_mesh(mesh_properties_name, mesh_name, **kwargs):
        
    mesh_properties: MeshProperty = create_properties_if_not_exists(mesh_name, mesh_properties_name)
    
    h_dist = calculate_face_properties.create_face_to_edge_distances(
            mesh_properties.faces_centroids,
            mesh_properties.adjacencies,
            mesh_properties.nodes_of_edges,
            mesh_properties.edges,
            mesh_properties.nodes_centroids,
            mesh_properties.bool_boundary_edges
        )
    mesh_properties.insert_or_update_data({'h_dist': h_dist})

    mesh_properties.export_data()

def preprocess_lpew2(mesh_properties_name):

    lpew = LpewWeight()
    mesh_properties = load_mesh_properties(mesh_properties_name)
    resp = lpew.preprocess(**mesh_properties.get_all_data())
    mesh_properties.insert_or_update_data(resp)

    h_dist = calculate_face_properties.create_face_to_edge_distances(
            mesh_properties.faces_centroids,
            mesh_properties.adjacencies,
            mesh_properties.nodes_of_edges,
            mesh_properties.edges,
            mesh_properties.nodes_centroids,
            mesh_properties.bool_boundary_edges
        )
    mesh_properties.insert_or_update_data({'h_dist': h_dist})
    mesh_properties.export_data()

def create_Tk(mesh_properties_name):
    """criar k normal e tangente
    """

    lpew = LpewWeight()
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    resp = lpew.create_Tk_points(
        mesh_properties['nodes_centroids'],
        mesh_properties['nodes_of_edges']
    )
    mesh_properties.insert_or_update_data(resp)
    mesh_properties.export_data()

def create_neta(mesh_properties_name):
    
    lpew = LpewWeight()
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    resp = lpew.create_neta(**mesh_properties.get_all_data())
    mesh_properties.insert_or_update_data(resp)
    mesh_properties.export_data()

def create_knt_vef(mesh_properties_name):
    
    lpew = LpewWeight()
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    resp = lpew.create_knt_barra_vef(**mesh_properties.get_all_data())
    mesh_properties.insert_or_update_data(resp)
    mesh_properties.export_data()

def create_zeta(mesh_properties_name):
    
    lpew = LpewWeight()
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    resp = lpew.create_zeta(**mesh_properties.get_all_data())
    mesh_properties.insert_or_update_data(resp)
    mesh_properties.export_data()

def create_lambda_barra(mesh_properties_name):
    
    lpew = LpewWeight()
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    resp = lpew.create_lambda_barra(**mesh_properties.get_all_data())
    mesh_properties.insert_or_update_data(resp)
    mesh_properties.export_data()

def create_lpew2_weights(mesh_properties_name):
    
    lpew = LpewWeight()
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    resp = lpew.create_lpew2_weights(**mesh_properties.get_all_data())
    mesh_properties.insert_or_update_data(resp)
    mesh_properties.export_data()
    print(mesh_properties.keys())

def sequence(mesh_properties_name, mesh_name, **kwargs):
    # mesh_properties_name = defnames.lpew2_test_mesh_prop_name
    # mesh_name = os.path.join(defpaths.mpfad_mesh_folder, 'mesh1_8_8.msh')
    
    initialize_mesh(mesh_properties_name, mesh_name)
    preprocess_lpew2(mesh_properties_name)
    create_Tk(mesh_properties_name)
    create_neta(mesh_properties_name)
    create_knt_vef(mesh_properties_name)
    create_zeta(mesh_properties_name)
    create_lambda_barra(mesh_properties_name)
    create_lpew2_weights(mesh_properties_name)

def verify_values(mesh_properties: MeshProperty):

    edges = mesh_properties['edges']
    edges_of_nodes = mesh_properties['edges_of_nodes']
    faces_of_nodes = mesh_properties['faces_of_nodes']
    adjacencies = mesh_properties['adjacencies']
    nodes_centroids = mesh_properties['nodes_centroids']
    faces_centroids = mesh_properties['faces_centroids']
    zeta = mesh_properties['zeta']
    tk_points = mesh_properties['tk_points']
    nodes_of_edges = mesh_properties['nodes_of_edges']
    kn_kt_theta_phi_vangle = mesh_properties['kn_kt_theta_phi_vangle']
    neta = mesh_properties['neta']
    kn_kt_barra = mesh_properties['kn_kt_barra']
    lambda_barra = mesh_properties['lambda_barra']
    lpew = LpewWeight()

    for node in mesh_properties.nodes:
        node_centroid = nodes_centroids[node]
        faces_node = faces_of_nodes[node]
        edges_node = edges_of_nodes[node]
        for face in faces_node:
            face_centroid = faces_centroids[face]
            edges_face = edges[(adjacencies[:, 0] == face) | (adjacencies[:, 1] == face)]
            edges_face = np.intersect1d(edges_face, edges_node)
            tkpoints_edges_face = tk_points[edges_face]
            tkpoints_edges_face = tkpoints_edges_face[nodes_of_edges[edges_face]==node]

            zetas = zeta['zeta'][
                (zeta['node_id'] == node)
            ]

            thetas = kn_kt_theta_phi_vangle['theta'][
                (kn_kt_theta_phi_vangle['node_id']==node) &
                (kn_kt_theta_phi_vangle['face_id']==face)
            ]

            phis = kn_kt_theta_phi_vangle['phi'][
                (kn_kt_theta_phi_vangle['node_id']==node) &
                (kn_kt_theta_phi_vangle['face_id']==face)
            ]

            netas = neta['neta'][
                (neta['node_id']==node) &
                (neta['face_id']==face)
            ] # OK

            vangles = kn_kt_theta_phi_vangle['v_angle'][
                (kn_kt_theta_phi_vangle['node_id']==node) &
                (kn_kt_theta_phi_vangle['face_id']==face)
            ] # ok

            ks = kn_kt_theta_phi_vangle[
                (kn_kt_theta_phi_vangle['node_id'] == node) &
                (kn_kt_theta_phi_vangle['face_id'] == face)
            ]

            ksbarra = kn_kt_barra[
                (kn_kt_barra['node_id'] == node) &
                (kn_kt_barra['face_id'] == face)
            ]

            lambda_node = lambda_barra['lambda_barra'][
                lambda_barra['node_id'] == node
            ]



            import pdb; pdb.set_trace()

def get_data_from_matfile(matfile, **kwargs):
    mdata = sio.loadmat(matfile)
    return mdata

def preprocess_mdata(mdata, mesh_properties: MeshProperty):

    nodes_centroids = mesh_properties['nodes_centroids']
    nodes = mesh_properties['nodes']
    faces = mesh_properties['faces']
    faces_centroids = mesh_properties['faces_centroids']

    mdata_coord_nos = mdata['coord_nos'][:, 0:2]
    mdata_coord_faces = mdata['coord_faces'][:, 0:2]

    delta = 1e-10
    all_nodes_id = np.repeat(-1, mdata_coord_nos.shape[0])
    all_faces_id = all_nodes_id.copy()

    for i, coord_no in enumerate(mdata_coord_nos):
        d1 = np.linalg.norm(nodes_centroids - coord_no, axis=1)
        test1 = d1 <= delta
        no_id = nodes[test1][0]
        all_nodes_id[i] = no_id

        d2 = np.linalg.norm(faces_centroids - mdata_coord_faces[i], axis=1)
        test2 = d2 <= delta
        face_id = faces[test2]
        all_faces_id[i] = face_id
    
    nodes_id_dirichlet = []
    values_dirichlet = []
    for data in mdata['dirichlet_nodes']:
        coord_no = data[0:2]
        d1 = np.linalg.norm(nodes_centroids - coord_no, axis=1)
        test1 = d1 <= delta
        no_id = nodes[test1][0]
        if no_id in nodes_id_dirichlet:
            continue
        nodes_id_dirichlet.append(no_id)
        values_dirichlet.append(data[2])
    
    nodes_id_dirichlet = np.array(nodes_id_dirichlet)
    values_dirichlet = np.array(values_dirichlet)
    mdata['nodes_id_dirichlet'] = nodes_id_dirichlet
    mdata['values_dirichlet_nodes'] = values_dirichlet

    mdata['node_id'] = all_nodes_id
    mdata['face_id'] = all_faces_id

    global_id_nodes = []
    global_id_faces = []

    for coord_no in mdata['coord']:
        d1 = np.linalg.norm(nodes_centroids - coord_no[0:2], axis=1)
        test1 = d1 <= delta
        no_id = nodes[test1][0]
        global_id_nodes.append(no_id)
    
    for coord_face in mdata['centelem']:
        d2 = np.linalg.norm(faces_centroids - coord_face[0:2], axis=1)
        test2 = d2 <= delta
        face_id = faces[test2]
        global_id_faces.append(face_id)
    
    global_id_nodes = np.array(global_id_nodes)
    global_id_faces = np.array(global_id_faces).flatten()

    mdata['global_id_nodes'] = global_id_nodes
    mdata['global_id_faces'] = global_id_faces

    mdata['rowM'] = mdata['rowM'].astype(np.int).flatten()
    mdata['colM'] = mdata['colM'].astype(np.int).flatten()
    mdata['dataM'] = mdata['dataM'].flatten()

    edges_centroids = mesh_properties.edges_centroids
    mdata_edge_ids = np.repeat(-1, len(edges_centroids))
    edges = mesh_properties.edges
    for i, edge_centroid in enumerate(mdata['edges_centroids']):
        d1 = np.linalg.norm(edge_centroid - edges_centroids, axis=1)
        test1 = d1 <= delta
        edge_id = edges[test1][0]
        mdata_edge_ids[i] = edge_id
    
    mdata['edges'] = mdata_edge_ids




def test_all_weights(mdata, mesh_properties: MeshProperty):

    mdata_weights = mdata['vecw']
    mdata_node_id = mdata['node_id']
    mdata_face_id = mdata['face_id']

    mesh_properties_weights = mesh_properties['nodes_weights']
    wtest = np.zeros(mesh_properties_weights.shape[0])

    for i, weight in enumerate(mdata_weights):
        node_mdata = mdata_node_id[i]
        face_mdata = mdata_face_id[i]

        test = (mesh_properties_weights['node_id']==node_mdata) & \
        (mesh_properties_weights['face_id']==face_mdata)

        wtest[test] = weight
    
    
    weights_error = np.abs(wtest - mesh_properties_weights['weight'])

    delta = 1e-10

    v1 = weights_error <= delta

    assert v1.sum() - len(mdata_weights) == 0 

    
def test_weights_prof_fernando():

    filenames = get_filenames()
    mdata = get_data_from_matfile(**filenames)
    initialize_mesh(**filenames)
    mesh_properties: MeshProperty = load_mesh_properties(filenames['mesh_properties_name'])
    get_lpew2_weights(mesh_properties)



    # sequence(**filenames)
    # mesh_properties: MeshProperty = load_mesh_properties(filenames['mesh_properties_name'])
    preprocess_mdata(mdata, mesh_properties)
    test_all_weights(mdata, mesh_properties)






