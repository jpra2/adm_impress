from files_to_test.lpew2_test_weights.info_test1 import get_filenames
from packs.mpfa_methods.weight_interpolation.test.test_lpew_weights import get_data_from_matfile, preprocess_mdata, initialize_mesh
from packs.manager.meshmanager import MeshProperty, load_mesh_properties
from scipy import sparse as sp
from packs.manager.boundary_conditions import BoundaryConditions
import numpy as np
from packs.mpfa_methods.flux_calculation.diamond_method import DiamondFluxCalculation, get_xi_params_ds_flux
from packs.mpfa_methods.weight_interpolation.lpew import get_lpew2_weights
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess


def run():
    diamond_flux = DiamondFluxCalculation()

    filenames = get_filenames()
    mdata = get_data_from_matfile(**filenames)
    initialize_mesh(**filenames)
    mesh_properties: MeshProperty = load_mesh_properties(filenames['mesh_properties_name'])
    preprocess_mdata(mdata, mesh_properties)
    n_faces = len(mesh_properties.faces)

    mpfa_preprocess = MpfaPreprocess()
    mpfa_preprocess.preprocess_data_lsds(mesh_properties)
    mpfa_preprocess.calculate_areas(mesh_properties)
    mpfa_preprocess.calculate_h_dist(mesh_properties)

    M = sp.csr_matrix((mdata['dataM'], (mdata['rowM']-1, mdata['colM']-1)), shape=(n_faces, n_faces))

    source2 = mdata['source'].flatten()
    fonte = mdata['fonte'].flatten()

    bc = BoundaryConditions()
    bc.set_boundary('dirichlet_nodes', mdata['nodes_id_dirichlet'], mdata['values_dirichlet_nodes'])
    bc.set_boundary('neumann_edges', np.array([]), np.array([]))

    permeability = mdata['kmap'][:, 1:5].reshape((n_faces, 2, 2))
    mesh_properties.insert_or_update_data({'permeability': permeability})
    mesh_properties.insert_or_update_data(
        {
            'neumann_edges': np.array([]),
            'neumann_edges_value': np.array([])
        }
    )

    get_lpew2_weights(mesh_properties)
    get_xi_params_ds_flux(mesh_properties)

    resp = diamond_flux.mount_problem(
            bc,
            **mesh_properties.get_all_data()
        )
    resp['source'] += fonte

    erro_source = resp['source'] - source2
    
    v1 = sp.find(resp['transmissibility'])
    v2 = sp.find(M)

    err_data = v1[2] - v2[2]
    print(np.linalg.norm(err_data))

    import pdb; pdb.set_trace()
    











    import pdb; pdb.set_trace()
    
    