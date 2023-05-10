from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import mesh_verify, linear_function, create_properties_if_not_exists
from packs.manager.meshmanager import load_mesh_properties, MeshProperty
from packs import defpaths
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
import numpy as np


def xi_verify(xi_alpha, mesh_properties: MeshProperty):
    bool_internal_edges = ~mesh_properties.bool_boundary_edges
    test = np.array(xi_alpha.tolist())[bool_internal_edges]
    line_sum = test.sum(axis=1)

def Skl_func(mesh_properties: MeshProperty):

    lsds_flux = LsdsFluxCalculation()
    Skl = lsds_flux.get_Skl(**mesh_properties.get_all_data())

def x_and_y_k_sigma_func(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    x_and_y_k_sigma = lsds.get_x_and_y_k_sigma(**mesh_properties.get_all_data())

def D_and_mi_func(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    D_and_mi, x_alpha, y_alpha = lsds.get_D_and_mi(**mesh_properties.get_all_data())

def xi_alpha_func(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    xi_alpha = lsds.get_xi_alpha_internal_edges(**mesh_properties.get_all_data())
    xi_verify(xi_alpha, mesh_properties)

def M_matrix_func(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    M_matrix = lsds.get_M_matrix(**mesh_properties.get_all_data())

def get_internal_edges_flux_params(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    internal_flux_params = lsds.get_internal_edges_flux_params(**mesh_properties.get_all_data())

def xi_params_test_func(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    xi_alpha = lsds.get_xi_alpha_internal_edges(**mesh_properties.get_all_data())
    internal_flux_params = lsds.get_internal_edges_flux_params(**mesh_properties.get_all_data())

    bool_internal_edges = ~mesh_properties.bool_boundary_edges
    test1 = np.array(xi_alpha.tolist())[bool_internal_edges]
    test2 = internal_flux_params[bool_internal_edges]
    line1_sum = test1.sum(axis=1)
    line2_sum = test2.sum(axis=1)

def get_boundary_weights_func(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    boundary_weights = lsds.get_boundary_edges_flux_params(**mesh_properties.get_all_data())

def get_all_edges_flux_params_func(mesh_properties: MeshProperty):
    lsds = LsdsFluxCalculation()
    all_edges_params = lsds.get_all_edges_flux_params(**mesh_properties.get_all_data())

def test_lsds_flux():
    mesh_properties_name = defpaths.mpfad_mesh_properties_name
    mesh_name = defpaths.mpfad_test_mesh
    mesh_verify(mesh_name)
    mesh_properties = create_properties_if_not_exists(mesh_name, mesh_properties_name)
    
    lsds = LsdsFluxCalculation()
    properties_updated = lsds.preprocess(**mesh_properties.get_all_data())
    mesh_properties.update_data(properties_updated)

    # Skl_func(mesh_properties)
    # x_and_y_k_sigma_func(mesh_properties)
    # D_and_mi_func(mesh_properties)
    # xi_alpha_func(mesh_properties)
    # M_matrix_func(mesh_properties)
    # get_internal_edges_flux_params(mesh_properties)
    # xi_params_test_func(mesh_properties)
    # get_boundary_weights_func(mesh_properties)
    get_all_edges_flux_params_func(mesh_properties)



