from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import mesh_verify, linear_function, create_properties_if_not_exists
from packs.manager.meshmanager import load_mesh_properties, MeshProperty
from packs import defpaths
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation



def mount_S_alpha_sigma(mesh_properties: MeshProperty):

    lsds_flux = LsdsFluxCalculation()

    lsds_flux.S_alpha_sigma(**mesh_properties.get_all_data())

    pass

def test_lsds_flux():
    mesh_properties_name = defpaths.mpfad_mesh_properties_name
    mesh_name = defpaths.mpfad_test_mesh
    mesh_verify(mesh_name)
    mesh_properties = create_properties_if_not_exists(mesh_name, mesh_properties_name)



    mount_S_alpha_sigma(mesh_properties)









    pass