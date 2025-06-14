## boundary mpfa
import os

nodes_pressure_prescription_name = 'dirichlet_nodes'
neumann_edges = 'neumann_edges'

mpfa_boundary_names = [neumann_edges, 
                       'dirichlet_nodes', 'neumann_nodes', 'dirichlet_edges',
                       'dirichlet_volumes', 'neumann_volumes',
                       'water_saturation_volumes', 'water_saturation_edges',
                       'injectors', 'producers', 'dirichlet_faces', 'initial_saturation',
                       'edges_injector', 'edges_producer']

tag_node_weight_test_sufix = '_nodes_weights_test_error'
lpew2_test_mesh_prop_name = 'lpew2_test'

nodes_weights_matrix_structure = 'nodes_weights_matrix_structure'

fine_primal_id = 'primal_id'
fine_dual_id = 'dual_id'
dual_volumes_str = 'dual_volumes'
dual_interation_region = 'dual_interation_region'
vertices_selected = 'vertices_selected'
edges_selected = 'edges_selected'
boundary_dual_interaction = 'boundary_of_interaction'
internal_dual_path = 'internal_dual_path'
dual_initial_ccs = 'initial_ccs'


list_op_toget = ['AMS', 'AMS-U', 'MsRSB']


def get_primal_id_name_by_level(level:int):
    return '_'.join([fine_primal_id, 'level' + str(level)])

def get_dual_id_name_by_level(level:int):
    return '_'.join([fine_dual_id, 'level' + str(level)])

def get_dual_volumes_name_by_level(level):
    return '_'.join([dual_volumes_str, 'level' + str(level)])

def get_dual_interation_region_name_by_level(level):
    return '_'.join([dual_interation_region, 'level' + str(level)])

def dual_ids(name:str):
    
    """
        my_dict = {
            'vertice_id': 3,
            'edge_id': 2,
            'face_id': 1,
            'internal_id': 0 
        }
    """
    
    my_dict = {
        'vertice_id': 3,
        'edge_id': 2,
        'face_id': 1,
        'internal_id': 0 
    }

    assert name in my_dict.keys()
    return my_dict[name]

def level_str(level):
    return '_level' + str(level)



