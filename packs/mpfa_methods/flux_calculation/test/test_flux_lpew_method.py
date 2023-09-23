from files_to_test.lpew2_test_weights.info_test1 import get_filenames
from packs.mpfa_methods.weight_interpolation.test.test_lpew_weights import initialize_mesh
from packs.manager.meshmanager import MeshProperty, load_mesh_properties
from packs.mpfa_methods.flux_calculation.lpew import get_xi_params_lpew_flux

def test_xi_params_lpew():
    filenames = get_filenames()
    initialize_mesh(**filenames)
    mesh_properties: MeshProperty = load_mesh_properties(filenames['mesh_properties_name'])
    get_xi_params_lpew_flux(mesh_properties)

    



