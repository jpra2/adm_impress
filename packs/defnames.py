## boundary mpfa
nodes_pressure_prescription_name = 'nodes_pressures'
neumann_edges = 'neumann_edges'
dirichlet_edges = 'dirichlet_edges'

mpfa_boundary_names = [nodes_pressure_prescription_name, neumann_edges, 
                       dirichlet_edges, 'dirichlet_nodes', 'neumann_nodes']

tag_node_weight_test_sufix = '_nodes_weights_test_error'
lpew2_test_mesh_prop_name = 'lpew2_test'

nodes_weights_matrix_structure = 'nodes_weights_matrix_structure'