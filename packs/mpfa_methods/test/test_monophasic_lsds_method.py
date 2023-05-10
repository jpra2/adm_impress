
import numpy as np
from packs.manager.arraydatamanager import SuperArrayManager
from packs.manager.generic_data import BoundaryConditions
from packs import defpaths
from packs.manager.meshmanager import MeshProperty
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists, mesh_verify
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation


def setup1():
    """Define monophasic problem with pressure presciption only
    """

    # create mesh
    mesh_test_name = defpaths.mpfad_test_mesh
    mesh_properties_name = defpaths.mpfad_mesh_properties_name
    mesh_verify(mesh_test_name)
    mesh_properties: MeshProperty = create_properties_if_not_exists(mesh_test_name, mesh_properties_name)
    
    # define nodes to calculate_weights
    mesh_properties.insert_data({'nodes_to_calculate': mesh_properties.nodes.copy()})

    ## create weights and xi params for flux calculation
    lsds = LsdsFluxCalculation()
    mesh_properties.update_data(
        lsds.preprocess(**mesh_properties.get_all_data())
    )
    nodes_weights = get_gls_nodes_weights(**mesh_properties.get_all_data())
    xi_params = lsds.get_all_edges_flux_params(**mesh_properties.get_all_data())

    mesh_properties.insert_data({
        'nodes_weights': nodes_weights,
        'xi_params': xi_params
    })

    mesh_properties.remove_data(['nodes_to_calculate'])
    





    

    import pdb; pdb.set_trace()





    















    import pdb; pdb.set_trace()










    pass





def test_monophasic_problem_with_pressure_prescription():
    """test the lsds flux calculation with gls vertex interpolation
       on a problem with only pressure prescription in 1D physics
    """

    setup1()




    pass
