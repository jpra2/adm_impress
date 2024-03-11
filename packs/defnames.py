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
dual_volumes_str = 'dual_volumes'
dual_interation_region = 'dual_interation_region'

def get_primal_id_name_by_level(level:int):
    return '_'.join(fine_primal_id, 'level' + str(level))

def get_dual_id_name_by_level(level:int):
    return '_'.join([fine_dual_id, 'level' + str(level)])

def get_dual_volumes_name_by_level(level):
    return '_'.join([dual_volumes_str, 'level' + str(level)])

def get_dual_interation_region_name_by_level(level):
    return '_'.join([dual_interation_region, 'level' + str(level)])

def dual_ids(name:str):
    my_dict = {
        'vertice_id': 3,
        'edge_id': 2,
        'face_id': 1,
        'inernal_id': 0 
    }

    assert name in my_dict.keys()
    return my_dict[name]
