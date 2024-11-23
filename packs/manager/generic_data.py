from packs.manager.arraydatamanager import SuperArrayManager

class SimulationData(SuperArrayManager):
    my_data_names = ['all_loops', 'all_vpi', 'all_cumulative_oil', 'all_cumulative_water', 'pressure_', 'saturation_']

class PrimalCoarseData(SuperArrayManager):
    my_data_names = [
        'adjacencies', 'nodes', 'faces', 'edges',
        'bool_boundary_edges', 'bool_intersect_edges',
        'map_nodes', 'map_faces', 'map_edges', 'nodes_weights_internal',
        'bool_boundary_nodes', 'bool_intersect_nodes', 'coarse_id',
        'nodes_to_calculate', 'nodes_of_nodes', 'edges_of_nodes',
        'faces_of_nodes', 'nodes_centroids', 'faces_centroids',
        'permeability', 'unitary_normal_edges', 'neumann_edges', 
        'neumann_edges_value', 'dirichlet_nodes', 'edges_multiplier',
        'dual_id', 'edges_flux', 'pressure', 'edges_dim', 'xi_params', 
        'xi_params_backup', 'nodes_weight_select', 'nodes_of_edges'
    ]