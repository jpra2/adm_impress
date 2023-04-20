from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import mesh_verify, linear_function, create_properties_if_not_exists
from packs.manager.meshmanager import load_mesh_properties, MeshProperty
from packs import defpaths
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation



def Skl_func(mesh_properties: MeshProperty):

    lsds_flux = LsdsFluxCalculation()
    Skl = lsds_flux.get_Skl(**mesh_properties.get_all_data())

def x_and_y_k_sigma_func(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    x_and_y_k_sigma = lsds.get_x_and_y_k_sigma(**mesh_properties.get_all_data())

def D_and_mi_func(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    D_and_mi = lsds.get_D_and_mi(**mesh_properties.get_all_data())

def epsilon_alpha_func(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    epsilon_alpha = lsds.get_epsilon_alpha(**mesh_properties.get_all_data())


def test_lsds_flux():
    mesh_properties_name = defpaths.mpfad_mesh_properties_name
    mesh_name = defpaths.mpfad_test_mesh
    mesh_verify(mesh_name)
    mesh_properties = create_properties_if_not_exists(mesh_name, mesh_properties_name)

    Skl_func(mesh_properties)
    x_and_y_k_sigma_func(mesh_properties)
    D_and_mi_func(mesh_properties)
    epsilon_alpha_func(mesh_properties)

