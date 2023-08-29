from packs.mpfa_methods.weight_interpolation.lpew2 import Lpew2Weight
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists
from packs import defpaths, defnames
from packs.manager.meshmanager import MeshProperty, load_mesh_properties
import os

mesh_properties_name = defnames.lpew2_test_mesh_prop_name

def initialize_mesh():
    global mesh_properties_name
    mesh_name = os.path.join(defpaths.mpfad_mesh_folder, 'mesh1_8_8.msh')
    mesh_properties = create_properties_if_not_exists(mesh_name, mesh_properties_name)

def preprocess_lpew2():
    global mesh_properties_name
    lpew2 = Lpew2Weight()
    mesh_properties = load_mesh_properties(mesh_properties_name)
    resp = lpew2.preprocess(**mesh_properties.get_all_data())
    mesh_properties.update_data(resp)
    mesh_properties.export_data()

def sequence():
    initialize_mesh()
    preprocess_lpew2()






