from packs.mpfa_methods.weight_interpolation.lpew import LpewWeight
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists
from packs import defpaths, defnames
from packs.manager.meshmanager import MeshProperty, load_mesh_properties
from packs.utils import calculate_face_properties
import os

mesh_properties_name = defnames.lpew2_test_mesh_prop_name

def initialize_mesh():
    global mesh_properties_name
    mesh_name = os.path.join(defpaths.mpfad_mesh_folder, 'mesh1_8_8.msh')
    mesh_properties = create_properties_if_not_exists(mesh_name, mesh_properties_name)

def preprocess_lpew2():
    global mesh_properties_name
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

def create_Tk():
    """criar k normal e tangente
    """
    global mesh_properties_name
    lpew = LpewWeight()
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    resp = lpew.create_Tk_points(
        mesh_properties['nodes_centroids'],
        mesh_properties['nodes_of_edges']
    )
    mesh_properties.insert_data(resp)
    mesh_properties.export_data()

def create_neta():
    global mesh_properties_name
    lpew = LpewWeight()
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    resp = lpew.create_neta(**mesh_properties.get_all_data())
    mesh_properties.insert_data(resp)
    mesh_properties.export_data()

def create_knt_vef():
    global mesh_properties_name
    lpew = LpewWeight()
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    resp = lpew.create_neta(**mesh_properties.get_all_data())


def sequence():
    # initialize_mesh()
    # preprocess_lpew2()
    # create_Tk()
    # create_neta()
    # create_knt_vef()






