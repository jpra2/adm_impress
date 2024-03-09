## boundary mpfa
nodes_pressure_prescription_name = 'dirichlet_nodes'
neumann_edges = 'neumann_edges'

mpfa_boundary_names = [neumann_edges, 
                       'dirichlet_nodes', 'neumann_nodes']

tag_node_weight_test_sufix = '_nodes_weights_test_error'
lpew2_test_mesh_prop_name = 'lpew2_test'

nodes_weights_matrix_structure = 'nodes_weights_matrix_structure'

fine_primal_id = 'primal_id'
fine_dual_id = 'dual_id'

def dual_ids(name:str):
    my_dict = {
        'vertice_id': 3,
        'edge_id': 2,
        'face_id': 1,
        'inernal_id': 0 
    }

    assert name in my_dict.keys()
    return my_dict[name]
