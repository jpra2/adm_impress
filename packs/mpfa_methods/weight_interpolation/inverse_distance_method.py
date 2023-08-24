import numpy as np

class InverseDistanceWeight2D:
    
    @staticmethod
    def get_weights(nodes_centroids, faces_of_nodes, faces_centroids, nodes):
        dtype = [('node_id', int), ('face_id', int), ('weight', float)]
        
        node_ids = []
        element_ids = []
        all_weights = []
        
        for node in nodes:
            elements = faces_of_nodes[node]
            centroid_node = nodes_centroids[node]
            centroids_elements = faces_centroids[elements]
            
            distances = np.linalg.norm(centroid_node - centroids_elements, axis=1)
            inv_distances = 1/distances
            weights = inv_distances/inv_distances.sum()
            
            node_ids.append(np.repeat(node, len(elements)))
            element_ids.append(elements)
            all_weights.append(weights)
        
        node_ids = np.concatenate(node_ids)
        element_ids = np.concatenate(element_ids)
        all_weights = np.concatenate(all_weights)

        nodes_weights = np.zeros(len(node_ids), dtype=dtype)
        nodes_weights['node_id'] = node_ids
        nodes_weights['face_id'] = element_ids
        nodes_weights['weight'] = weights
        
        return {
            'nodes_weights': nodes_weights
        }
    
    