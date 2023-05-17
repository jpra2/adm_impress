
import numpy as np
from packs.manager.boundary_conditions import BoundaryConditions
from packs import defpaths
from packs.manager.meshmanager import MeshProperty
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists, mesh_verify
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.utils.utils_old import get_box
from scipy.sparse.linalg import spsolve
from packs.manager.mesh_data import MeshData
import os

def setup1():
    """Define monophasic problem with pressure presciption only
    """

    # create mesh
    mesh_test_name = defpaths.mpfad_test_mesh
    mesh_properties_name = defpaths.mpfad_mesh_properties_name
    mesh_verify(mesh_test_name)
    mesh_properties: MeshProperty = create_properties_if_not_exists(mesh_test_name, mesh_properties_name)
    
    # # define nodes to calculate_weights
    # mesh_properties.insert_data({'nodes_to_calculate': mesh_properties.nodes.copy()})

    # ## create weights and xi params for flux calculation
    # lsds = LsdsFluxCalculation()
    # mesh_properties.update_data(
    #     lsds.preprocess(**mesh_properties.get_all_data())
    # )

    # mesh_properties.insert_data({
    #     'nodes_weights': get_gls_nodes_weights(**mesh_properties.get_all_data()),
    #     'xi_params': lsds.get_all_edges_flux_params(**mesh_properties.get_all_data())
    # })

    # mesh_properties.remove_data(['nodes_to_calculate'])
    # mesh_properties.export_data()

    ## define boundary conditions
    problem_name = 'monophasic_1'
    bc = BoundaryConditions()
    bc.insert_name(problem_name)

    # left_value = 100
    # right_value = 0

    # vcentroids = mesh_properties.nodes_centroids

    # delta = 1e-5
    # x_left = vcentroids[:,0].min()
    # x_right = vcentroids[:,0].max()

    # nodes_left = mesh_properties.nodes[vcentroids[:,0] < x_left + delta]
    # nodes_right = mesh_properties.nodes[vcentroids[:,0] > x_right - delta]

    # pressure_left = np.repeat(left_value, len(nodes_left))
    # pressure_right = np.repeat(right_value, len(nodes_right))

    # nodes_bc = np.concatenate([nodes_left, nodes_right])
    # pressures_bc = np.concatenate([pressure_left, pressure_right])

    # dtype_bc_array = [('id', np.uint64), ('value', np.float64)]
    # bc_array = np.zeros(len(nodes_bc), dtype=dtype_bc_array)
    # bc_array['id'][:] = nodes_bc
    # bc_array['value'][:] = pressures_bc
    
    # bc.insert_data({'nodes_pressures': bc_array})

    # bc.export_data()
    bc.load_data()
    
    
    lsds = LsdsFluxCalculation()
    resp = lsds.mount_problem(
        mesh_properties.nodes_weights,
        mesh_properties.xi_params,
        mesh_properties.faces,
        mesh_properties.edges,
        mesh_properties.bool_boundary_edges,
        mesh_properties.adjacencies,
        bc,
        mesh_properties.nodes,
        mesh_properties.bool_boundary_nodes,
        mesh_properties.nodes_of_edges
    )
    
    pressure = spsolve(resp['transmissibility'], resp['source'])
    mesh_path = os.path.join(defpaths.mesh, defpaths.mpfad_test_mesh)
    mesh_data = MeshData(mesh_path=mesh_path)
    mesh_data.create_tag('pressure')
    mesh_data.insert_tag_data('pressure', pressure, 'faces', mesh_properties.faces)
    mesh_data.export_all_elements_type_to_vtk('test_pressure', 'faces')


def test_monophasic_problem_with_pressure_prescription():
    """test the lsds flux calculation with gls vertex interpolation
       on a problem with only pressure prescription in 1D physics
    """

    setup1()
    import pdb; pdb.set_trace()



    pass
