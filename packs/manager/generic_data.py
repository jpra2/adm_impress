from packs.manager.arraydatamanager import SuperArrayManager
import numpy as np
from packs.utils.utils_old import remap_values
import os
from packs import defpaths

class SimulationData(SuperArrayManager):
    my_data_names = ['all_loops', 'all_vpi', 'all_cumulative_oil', 'all_cumulative_water', 'pressure_', 'saturation_', 'water_flux', 'oil_flux']

    @property
    def class_path(self):
        return os.path.join(defpaths.data_simulation, self.class_name() + '_' +  self.name + '.npz')

    @property
    def label(self):
        try: 
            return self['label'][0]
        except KeyError:
            return ''
    
    @property
    def name(self):
        return self['name'].copy()[0]


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
        'xi_params_backup', 'nodes_weight_select', 'nodes_of_edges', 
        'other_side_flux', 'local_boundary_edge_flux_for_test', 'local_saturation_test',
        'max_delta_sat', 'percent_var_flux', 'boundary_nodes_weights', 'neumann_weights'
    ]

    def get_local_nodes_weights(self, global_nodes_weight, map_from_nodes, map_to_nodes, map_from_faces, map_to_faces, local_bool_boundary_nodes, global_nodes_to_calculate):
        
        if np.any(np.isin(map_from_nodes, global_nodes_to_calculate)):
            pass
        else:
            return 
        nodes_weight = global_nodes_weight
        test1 = np.isin(nodes_weight['node_id'], map_from_nodes)
        test2 =  np.isin(nodes_weight['face_id'], map_from_faces)
        test3 = test1 & test2

        local_nodes_weight = nodes_weight[test3]

        local_nodes_weight['node_id'][:] = remap_values(map_from_nodes, map_to_nodes, local_nodes_weight['node_id'])
        local_nodes_weight['face_id'][:] = remap_values(map_from_faces, map_to_faces, local_nodes_weight['face_id'])

        test4 = np.isin(local_nodes_weight['node_id'], map_to_nodes[local_bool_boundary_nodes])
        test4 = ~test4
        local_nodes_weight = local_nodes_weight[test4]
        
        self.insert_or_update_data({
            self.my_data_names[9]: local_nodes_weight
        })
    
    def get_local_nodes_weights_explicit(self, global_nodes_weight, map_from_nodes, map_to_nodes, map_from_faces, map_to_faces, local_bool_boundary_nodes):
        nodes_weight = global_nodes_weight
        test1 = np.isin(nodes_weight['node_id'], map_from_nodes)
        test2 =  np.isin(nodes_weight['face_id'], map_from_faces)
        test3 = test1 & test2

        local_nodes_weight = nodes_weight[test3]

        local_nodes_weight['node_id'][:] = remap_values(map_from_nodes, map_to_nodes, local_nodes_weight['node_id'])
        local_nodes_weight['face_id'][:] = remap_values(map_from_faces, map_to_faces, local_nodes_weight['face_id'])

        test4 = np.isin(local_nodes_weight['node_id'], map_to_nodes[local_bool_boundary_nodes])
        test4 = ~test4
        local_nodes_weight = local_nodes_weight[test4]
        
        return local_nodes_weight, test3
    
    def update_all_nodes_weights(self, nodes_to_calculate, boundary_nodes_weights_all):
        
        local_nodes_weights = self['nodes_weights_internal'].copy()
    
        if np.all(np.isin(self['nodes'][self['bool_boundary_nodes']], nodes_to_calculate)):
            ## se foi o primeiro looop ou se todos os nos do contorno do primal precisaram ser atualizados
            self.insert_or_update_data({
                'boundary_nodes_weights': boundary_nodes_weights_all['nodes_weights'],
                'neumann_weights': boundary_nodes_weights_all['neumann_weights']
            })
        else:
            ## se apenas alguns precisaram ser atualizados, atualiza-se apenas os modificados
            bnodes_weight = self['boundary_nodes_weights']
            nodes_to_modify_bool = np.isin(bnodes_weight['node_id'], nodes_to_calculate)
            nodes_to_maintain_bool = ~nodes_to_modify_bool
            new_bnodes_weight = bnodes_weight[nodes_to_maintain_bool]
            
            local_bnodes_weights = np.hstack([new_bnodes_weight, boundary_nodes_weights_all['nodes_weights']])
            
            neumann_weights = self['neumann_weights']
            nodes_to_modify_bool_neumann = np.isin(neumann_weights['node_id'], nodes_to_calculate)
            nodes_to_maintain_bool_neumann = ~nodes_to_modify_bool_neumann
            new_neumann_weights = neumann_weights[nodes_to_maintain_bool_neumann]
            
            local_neumann_weights = np.hstack([new_neumann_weights, boundary_nodes_weights_all['neumann_weights']])
            
            self.insert_or_update_data({
                'boundary_nodes_weights': local_bnodes_weights,
                'neumann_weights': local_neumann_weights
            })

        # local_nodes_weights = np.hstack([local_nodes_weights, boundary_nodes_weights['nodes_weights']])
        local_nodes_weights = np.hstack([local_nodes_weights, self['boundary_nodes_weights']])
        
        self.insert_or_update_data({
            'nodes_weights': local_nodes_weights
        })