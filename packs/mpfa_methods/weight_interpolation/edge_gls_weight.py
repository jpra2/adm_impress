import numpy as np
import scipy.sparse as sp
from packs.utils.utils_old import get_shortest_path

class EdgeGlsWeight:

    @staticmethod
    def get_Rmatrix():
        theta = np.pi/2

        R = np.array([
            [np.cos(theta), np.sin(theta)],
            [-np.sin(theta), np.cos(theta)]
        ])

        return R

    def get_T_face(self, permeability, K, L, unitary_normal, R):
        """
        Retorna a transformacao do gradiente em L em funcao do gradiente em K
        """
        permK = permeability[K]
        permL = permeability[L]

        M1 = np.array(
            [R.dot(unitary_normal).T[0],
            permL.dot(unitary_normal).T[0]]
        )

        M2 = np.array(
            [R.dot(unitary_normal).T[0],
            permK.dot(unitary_normal).T[0]]
        )

        M = np.linalg.inv(M1).dot(M2)
        return M
    
    def update_W_from_node(
            self,
            edges_of_nodes: np.ndarray,
            node,
            internal_edges,
            edge,
            adjacencies: np.ndarray,
            nodes_centroids,
            faces_centroids,
            edge_centroid,
            W,
            K_face,
            L_face,
            faces_to_weight,
            all_T_internal_edges,
            **kwargs
    ):
        
        adjacencies_internal_edges = adjacencies[internal_edges]

        edges_of_node = edges_of_nodes[node][np.isin(edges_of_nodes[node], internal_edges)]
        v1 = np.where(edges_of_node == edge)[0]
        # edges_of_node = np.roll(edges_of_node, edges_of_node.shape[0]-v1[0], axis=0)[1:]
        try:
            edges_of_node = np.roll(edges_of_node, edges_of_node.shape[0]-v1[0], axis=0)
        except IndexError:
            import pdb; pdb.set_trace()

        adjacencies_edges_of_node = adjacencies[edges_of_node]
        unique_adjacencies = np.unique(adjacencies_edges_of_node.flatten())
        map_local = np.arange(unique_adjacencies.shape[0])
        local_adjacencies = adjacencies_edges_of_node.copy()
        
        for i in unique_adjacencies:
            local_adjacencies[local_adjacencies==i] = map_local[unique_adjacencies==i]
        
        data = np.ones(local_adjacencies.shape[0])
        n = map_local.shape[0]
        local_graph = sp.csc_matrix((data, (local_adjacencies[:,0], local_adjacencies[:,1])), shape=(n ,n))
        
        test = np.isin(unique_adjacencies, [K_face, L_face])
        test = ~test
        to_iterate = unique_adjacencies[test]

        K_face_local = map_local[unique_adjacencies==K_face][0]
        for face in to_iterate:
            local_face = map_local[unique_adjacencies==face][0]
            local_path = get_shortest_path(local_graph, K_face_local, local_face)[1:]
            global_path = local_path.copy()
            for i in local_path:
                global_path[local_path==i] = unique_adjacencies[map_local==i]
            
            face2 = global_path[0]
            test_edge = (
                ((adjacencies_internal_edges[:, 0]==K_face) & (adjacencies_internal_edges[:, 1]==face2)) |
                ((adjacencies_internal_edges[:, 1]==K_face) & (adjacencies_internal_edges[:, 0]==face2))
            )

            T_total = all_T_internal_edges[test_edge][0]
            if K_face == adjacencies_internal_edges[test_edge, 0]:
                pass
            elif K_face == adjacencies_internal_edges[test_edge, 1]:
                T_total = np.linalg.inv(T_total)
            else:
                raise ValueError

            if global_path.shape[0] > 1:
                K_face2 = face2
                for face2 in global_path[1:]:
                    test_edge = (
                        ((adjacencies_internal_edges[:, 0]==K_face2) & (adjacencies_internal_edges[:, 1]==face2)) |
                        ((adjacencies_internal_edges[:, 1]==K_face2) & (adjacencies_internal_edges[:, 0]==face2))
                    )
                    T2 = all_T_internal_edges[test_edge][0]
                    if K_face2 == adjacencies_internal_edges[test_edge, 0]:
                        pass
                    elif K_face2 == adjacencies_internal_edges[test_edge, 1]:
                        T2 = np.linalg.inv(T2)
                    else:
                        raise ValueError
                    
                    T_total = T_total.dot(T2)
                    K_face2 = face2        

            ww = (nodes_centroids[node] - edge_centroid) + (faces_centroids[face] - nodes_centroids[node]).dot(T_total) 
            W.append(ww)
            faces_to_weight.append(face)
    
    def update_W_from_node_in_boundary(
            self,
            edges_of_nodes: np.ndarray,
            node,
            internal_edges,
            edge,
            adjacencies: np.ndarray,
            nodes_centroids,
            faces_centroids,
            edge_centroid,
            W,
            K_face,
            L_face,
            faces_to_weight,
            all_T_internal_edges,
            **kwargs
    ):
        
        adjacencies_internal_edges = adjacencies[internal_edges]

        edges_of_node = edges_of_nodes[node][np.isin(edges_of_nodes[node], internal_edges)]

        adjacencies_edges_of_node = adjacencies[edges_of_node]
        unique_adjacencies = np.unique(adjacencies_edges_of_node.flatten())
        map_local = np.arange(unique_adjacencies.shape[0])
        local_adjacencies = adjacencies_edges_of_node.copy()
        
        for i in unique_adjacencies:
            local_adjacencies[local_adjacencies==i] = map_local[unique_adjacencies==i]
        
        data = np.ones(local_adjacencies.shape[0])
        n = map_local.shape[0]
        local_graph = sp.csc_matrix((data, (local_adjacencies[:,0], local_adjacencies[:,1])), shape=(n ,n))
        
        test = np.isin(unique_adjacencies, [K_face, L_face])
        test = ~test
        to_iterate = unique_adjacencies[test]

        K_face_local = map_local[unique_adjacencies==K_face][0]
        for face in to_iterate:
            local_face = map_local[unique_adjacencies==face][0]
            local_path = get_shortest_path(local_graph, K_face_local, local_face)[1:]
            global_path = local_path.copy()
            for i in local_path:
                global_path[local_path==i] = unique_adjacencies[map_local==i]
            
            face2 = global_path[0]
            test_edge = (
                ((adjacencies_internal_edges[:, 0]==K_face) & (adjacencies_internal_edges[:, 1]==face2)) |
                ((adjacencies_internal_edges[:, 1]==K_face) & (adjacencies_internal_edges[:, 0]==face2))
            )

            T_total = all_T_internal_edges[test_edge][0]
            if K_face == adjacencies_internal_edges[test_edge, 0]:
                pass
            elif K_face == adjacencies_internal_edges[test_edge, 1]:
                T_total = np.linalg.inv(T_total)
            else:
                raise ValueError

            if global_path.shape[0] > 1:
                K_face2 = face2
                for face2 in global_path[1:]:
                    test_edge = (
                        ((adjacencies_internal_edges[:, 0]==K_face2) & (adjacencies_internal_edges[:, 1]==face2)) |
                        ((adjacencies_internal_edges[:, 1]==K_face2) & (adjacencies_internal_edges[:, 0]==face2))
                    )
                    T2 = all_T_internal_edges[test_edge][0]
                    if K_face2 == adjacencies_internal_edges[test_edge, 0]:
                        pass
                    elif K_face2 == adjacencies_internal_edges[test_edge, 1]:
                        T2 = np.linalg.inv(T2)
                    else:
                        raise ValueError
                    
                    T_total = T_total.dot(T2)
                    K_face2 = face2        

            ww = (nodes_centroids[node] - edge_centroid) + (faces_centroids[face] - nodes_centroids[node]).dot(T_total) 
            W.append(ww)
            faces_to_weight.append(face)

    def update_W_from_node_in_boundary_neumann(
            self,
            edges_of_nodes: np.ndarray,
            node,
            internal_edges,
            edge,
            adjacencies: np.ndarray,
            nodes_centroids,
            faces_centroids,
            edge_centroid,
            W,
            K_face,
            L_face,
            faces_to_weight,
            all_T_internal_edges,
            **kwargs
    ):
        
        adjacencies_internal_edges = adjacencies[internal_edges]

        edges_of_node = edges_of_nodes[node][np.isin(edges_of_nodes[node], internal_edges)]

        adjacencies_edges_of_node = adjacencies[edges_of_node]
        unique_adjacencies = np.unique(adjacencies_edges_of_node.flatten())
        map_local = np.arange(unique_adjacencies.shape[0])
        local_adjacencies = adjacencies_edges_of_node.copy()
        
        for i in unique_adjacencies:
            local_adjacencies[local_adjacencies==i] = map_local[unique_adjacencies==i]
        
        data = np.ones(local_adjacencies.shape[0])
        n = map_local.shape[0]
        local_graph = sp.csc_matrix((data, (local_adjacencies[:,0], local_adjacencies[:,1])), shape=(n ,n))
        
        test = np.isin(unique_adjacencies, [K_face, L_face])
        test = ~test
        to_iterate = unique_adjacencies[test]

        K_face_local = map_local[unique_adjacencies==K_face][0]
        for face in to_iterate:
            local_face = map_local[unique_adjacencies==face][0]
            local_path = get_shortest_path(local_graph, K_face_local, local_face)[1:]
            global_path = local_path.copy()
            for i in local_path:
                global_path[local_path==i] = unique_adjacencies[map_local==i]
            
            face2 = global_path[0]
            test_edge = (
                ((adjacencies_internal_edges[:, 0]==K_face) & (adjacencies_internal_edges[:, 1]==face2)) |
                ((adjacencies_internal_edges[:, 1]==K_face) & (adjacencies_internal_edges[:, 0]==face2))
            )

            T_total = all_T_internal_edges[test_edge][0]
            if K_face == adjacencies_internal_edges[test_edge, 0]:
                pass
            elif K_face == adjacencies_internal_edges[test_edge, 1]:
                T_total = np.linalg.inv(T_total)
            else:
                raise ValueError

            if global_path.shape[0] > 1:
                K_face2 = face2
                for face2 in global_path[1:]:
                    test_edge = (
                        ((adjacencies_internal_edges[:, 0]==K_face2) & (adjacencies_internal_edges[:, 1]==face2)) |
                        ((adjacencies_internal_edges[:, 1]==K_face2) & (adjacencies_internal_edges[:, 0]==face2))
                    )
                    T2 = all_T_internal_edges[test_edge][0]
                    if K_face2 == adjacencies_internal_edges[test_edge, 0]:
                        pass
                    elif K_face2 == adjacencies_internal_edges[test_edge, 1]:
                        T2 = np.linalg.inv(T2)
                    else:
                        raise ValueError
                    
                    T_total = T_total.dot(T2)
                    K_face2 = face2        

            ww = (nodes_centroids[node] - edge_centroid) + (faces_centroids[face] - nodes_centroids[node]).dot(T_total) 
            W.append(ww)
            faces_to_weight.append(face)

    def create_weights(
            self,
            edges,
            faces_centroids,
            nodes_centroids,
            edges_of_nodes,
            nodes_of_edges,
            unitary_normal_edges,
            adjacencies,
            bool_boundary_edges,
            permeability,
            edges_centroids,
            edges_dim,
            neumann_edges,
            neumann_edges_value,
            **kwargs
    ):
        
        bool_internal_edges = ~bool_boundary_edges
        internal_edges = edges[bool_internal_edges]
        # R = self.get_Rmatrix()

        edges_weight = []
        faces_weight = []
        all_weights = []
        all_neumann_weights = np.zeros(neumann_edges.shape[0])

        # all_T_internal_edges = np.array([self.get_T_face(permeability, adjacencies[i,0], adjacencies[i,1], unitary_normal_edges[i].reshape(2, 1), R) for i in internal_edges])
        all_T_internal_edges = np.array([self.get_T_face(permeability, adjacencies[i,0], adjacencies[i,1], unitary_normal_edges[i].reshape(2, 1), self.get_Rmatrix()) for i in internal_edges])
        
        for i, edge in enumerate(internal_edges):
            nodes_of_edge = nodes_of_edges[edge]
            faces_adjs_edge = adjacencies[edge]
            edge_centroid = edges_centroids[edge]

            node1 = nodes_of_edge[0]

            K_face = faces_adjs_edge[0]
            L_face = faces_adjs_edge[1]

            T_edge = all_T_internal_edges[i]

            W = [
                (nodes_centroids[node1] - edge_centroid) + (faces_centroids[K_face] - nodes_centroids[node1]),
                (nodes_centroids[node1] - edge_centroid) + (faces_centroids[L_face] - nodes_centroids[node1]).dot(T_edge)
            ]

            faces_to_weight = [K_face, L_face]

            for node in nodes_of_edge:
                self.update_W_from_node(
                    edges_of_nodes,
                    node,
                    internal_edges,
                    edge,
                    adjacencies,
                    nodes_centroids,
                    faces_centroids,
                    edge_centroid,
                    W,
                    K_face,
                    L_face,
                    faces_to_weight,
                    all_T_internal_edges
                )
            
            W = np.array(W)
            W = np.hstack([np.ones((W.shape[0], 1)), W])
            faces_to_weight = np.array(faces_to_weight)
            weights = np.array([1, 0, 0]).dot(np.linalg.inv(W.T.dot(W)).dot(W.T))

            edges_weight.append(np.repeat(edge, weights.shape[0]))
            faces_weight.append(faces_to_weight)
            all_weights.append(weights)
        
        bedges = edges[bool_boundary_edges]
        dirichlet_edges = np.setdiff1d(bedges, neumann_edges)

        L_face = -1
        for i, edge in enumerate(dirichlet_edges):
            nodes_of_edge = nodes_of_edges[edge]
            faces_adjs_edge = adjacencies[edge]
            edge_centroid = edges_centroids[edge]

            node1 = nodes_of_edge[0]
            K_face = faces_adjs_edge[0]

            W = [
                (nodes_centroids[node1] - edge_centroid) + (faces_centroids[K_face] - nodes_centroids[node1])
            ]

            faces_to_weight = [K_face]

            for node in nodes_of_edge:
                self.update_W_from_node_in_boundary(
                    edges_of_nodes,
                    node,
                    internal_edges,
                    edge,
                    adjacencies,
                    nodes_centroids,
                    faces_centroids,
                    edge_centroid,
                    W,
                    K_face,
                    L_face,
                    faces_to_weight,
                    all_T_internal_edges
                )

            W = np.array(W)
            W = np.hstack([np.ones((W.shape[0], 1)), W])
            faces_to_weight = np.array(faces_to_weight)
            weights = np.array([1, 0, 0]).dot(np.linalg.inv(W.T.dot(W)).dot(W.T))

            edges_weight.append(np.repeat(edge, weights.shape[0]))
            faces_weight.append(faces_to_weight)
            all_weights.append(weights)
        
        for i, edge in enumerate(neumann_edges):
            nodes_of_edge = nodes_of_edges[edge]
            faces_adjs_edge = adjacencies[edge]
            edge_centroid = edges_centroids[edge]

            node1 = nodes_of_edge[0]
            K_face = faces_adjs_edge[0]

            W = [
                (nodes_centroids[node1] - edge_centroid) + (faces_centroids[K_face] - nodes_centroids[node1])
            ]

            faces_to_weight = [K_face]

            for node in nodes_of_edge:
                self.update_W_from_node_in_boundary(
                    edges_of_nodes,
                    node,
                    internal_edges,
                    edge,
                    adjacencies,
                    nodes_centroids,
                    faces_centroids,
                    edge_centroid,
                    W,
                    K_face,
                    L_face,
                    faces_to_weight,
                    all_T_internal_edges
                )

            W = np.array(W)
            W = np.hstack([np.ones((W.shape[0], 1)), W])
            ni = W.shape[0]

            N = np.zeros((ni + 1, ni))
            N[range(ni), range(ni)] = 1

            permK_face = permeability[K_face]
            unitary_normal_edge = unitary_normal_edges[edge]
            edge_dim = edges_dim[edge]
            neumann_value = neumann_edges_value[i]

            vdot = unitary_normal_edge.dot(permK_face)
            presc_value = edge_dim*neumann_value
            # presc_value = neumann_value

            F = np.zeros(ni+1)
            F[-1] = presc_value

            W = np.vstack([W, np.zeros((1, 3))])
            W[-1, 1:] = vdot

            faces_to_weight = np.array(faces_to_weight)
            InvW = np.linalg.inv(W.T.dot(W))
            weights = np.array([1, 0, 0]).dot(InvW.dot(W.T).dot(N))

            edges_weight.append(np.repeat(edge, weights.shape[0]))
            faces_weight.append(faces_to_weight)
            all_weights.append(weights)

            neumann_weight = np.array([0, 0, 1]).dot(InvW.dot(W.T).dot(F))
            all_neumann_weights[i] = neumann_weight

        edges_weight = np.concatenate(edges_weight)
        faces_weight = np.concatenate(faces_weight)
        all_weights = np.concatenate(all_weights)

        return edges_weight, faces_weight, all_weights, all_neumann_weights


