from files_to_test.lpew2_test_weights.info_test1 import get_filenames
from packs.mpfa_methods.weight_interpolation.test.test_lpew_weights import get_data_from_matfile, preprocess_mdata, initialize_mesh
from packs.manager.meshmanager import MeshProperty, load_mesh_properties
from scipy import sparse as sp
from packs.manager.boundary_conditions import BoundaryConditions
import numpy as np
from packs.mpfa_methods.flux_calculation.diamond_method import DiamondFluxCalculation, get_xi_params_ds_flux
from packs.mpfa_methods.weight_interpolation.lpew import get_lpew2_weights
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess
from packs.mpfa_methods.weight_interpolation.test.test_lpew_weights import test_all_weights

def test_internal_faces(data_mdata, data_mesh, nodes, bool_boundary_nodes, nodes_of_nodes, edges_of_nodes, adjacencies, **kwargs):
    
    bnodes = nodes[bool_boundary_nodes]
    nodes_adj_bnodes = np.unique(np.concatenate(nodes_of_nodes[bnodes]))
    test = np.isin(nodes, nodes_adj_bnodes)
    test2 = ~test

    internal_edges = np.unique(np.concatenate(edges_of_nodes[test2]))

    for edge in internal_edges:
        faces_adj = adjacencies[edge]
        for face in faces_adj:

            test1 = data_mdata[0] == face
            test2 = data_mesh[0] == face

            lines_mdata = data_mdata[0][test1]
            lines_mesh = data_mesh[0][test2]

            cols_mdata = data_mdata[1][test1]
            cols_mesh = data_mesh[1][test2]

            dataa_mdata = data_mdata[2][test1]
            dataa_mesh = data_mesh[2][test2]

def test_kn_bedges(mesh_properties: MeshProperty, mdata):

    delta = 1e-8
    kn_kt_mesh = mesh_properties['kn_kt_dsflux']
    kn_data = mdata['kn_data']
    edge_ids_data = mdata['edges']
    edge_ids_mesh = mesh_properties.edges
    adjacencies = mesh_properties.adjacencies
    biedges = ~mesh_properties.bool_boundary_edges
    edges_dim = mesh_properties.edges_dim

    for edge in mesh_properties.edges[mesh_properties.bool_boundary_edges]:
        kn_edge_mesh = kn_kt_mesh['kn'][(kn_kt_mesh['edge_id']==edge) & (kn_kt_mesh['face_id']==adjacencies[edge, 0])]
        kn_edge_mdata = kn_data[edge_ids_data==edge]
        kt_edge_mesh = kn_kt_mesh['kt'][(kn_kt_mesh['edge_id']==edge) & (kn_kt_mesh['face_id']==adjacencies[edge, 0])]
        kt_edge_mdata = mdata['kt_data'][edge_ids_data==edge]

        assert abs(kn_edge_mdata.flatten()[0] - kn_edge_mesh.flatten()[0]) <= delta
        assert abs(kt_edge_mdata.flatten()[0] - kt_edge_mesh.flatten()[0]) <= delta

    
    for edge in mesh_properties.edges[biedges]:
        kappa_mesh = abs(mesh_properties['kappa_D_dsflux']['kappa'][edge]*edges_dim[edge])
        kappa_mdata = np.abs(mdata['ked_data'][edge_ids_data==edge])
        D_mesh = abs(mesh_properties['kappa_D_dsflux']['D'][edge])
        D_mdata = np.abs(mdata['ded_data'][edge_ids_data==edge])
        
        try:
            assert abs(kappa_mesh - kappa_mdata.flatten()[0]) <= delta
            assert abs(D_mesh - D_mdata.flatten()[0]) <= delta
        except:
            import pdb; pdb.set_trace()


    import pdb; pdb.set_trace()
    pass







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

    # source2 = mdata['source'].flatten()
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

    test_all_weights(mdata, mesh_properties)

    resp = diamond_flux.mount_problem(
            bc,
            **mesh_properties.get_all_data()
        )
    resp['source'] += fonte

    test_kn_bedges(mesh_properties, mdata)

    # erro_source = resp['source'] - source2
    
    v1 = sp.find(resp['transmissibility'])
    v2 = sp.find(M)



    err_data = v1[2] - v2[2]
    print(np.linalg.norm(err_data))
    import pdb; pdb.set_trace()

    import pdb; pdb.set_trace()
    











    import pdb; pdb.set_trace()
    
    