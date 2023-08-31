import numpy as np
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation

class LpewWeight:
    """A LINEARITY PRESERVING CELL-CENTERED
        SCHEME FOR THE HETEROGENEOUS AND
        ANISOTROPIC DIFFUSION EQUATIONS ON
        GENERAL MESHES

    """

    @staticmethod
    def preprocess(
        nodes_centroids,
        unitary_normal_edges,
        nodes_of_edges,
        edges, 
        **kwargs
    ):

        return LsdsFluxCalculation.define_A_B_points_of_edges(
            nodes_centroids,
            unitary_normal_edges,
            nodes_of_edges,
            edges
        )
    
    def get_Rmatrix(self):
        theta = np.pi/2

        R = np.array([
            [np.cos(theta), np.sin(theta)],
            [-np.sin(theta), np.cos(theta)]
        ])

        return R

    def create_Tk_points(self, nodes_centroids, nodes_of_edges, tau=0.5, **kwargs):
        """create tk_points

        Args:
            nodes_centroids (_type_): _description_
            nodes_of_edges (_type_): _description_
            tau (float, optional): _description_. Defaults to 0.5.

        Returns:
            _type_: retorna o vetor tal que:
            tk_points[edge, 0] = Tk of node nodes_of_edges[edge, 0]
            tk_points[edge, 1] = Tk of node nodes_of_edges[edge, 1]
        """

        dim = nodes_of_edges.shape
        tk = np.zeros((dim[0], dim[1], 2))
        centroids_nodes_edges = nodes_centroids[nodes_of_edges]

        v1 = tau*centroids_nodes_edges[:, 0] + (1-tau)*centroids_nodes_edges[:, 1]
        v2 = tau*centroids_nodes_edges[:, 1] + (1-tau)*centroids_nodes_edges[:, 0]

        tk[:, 0] = v1
        tk[:, 1] = v2

        return {'tk_points': tk}

    def create_neta(self, tk_points, nodes_centroids, nodes_of_edges, h_dist, edges, adjacencies, bool_boundary_edges,  **kwargs):

        all_neta = []
        all_face_id = []
        all_edge_id = []
        all_vertice_id = []
        
        biedges = ~bool_boundary_edges

        for edge in edges[biedges]:
            nodes_edge = nodes_of_edges[edge]
            centroids_nodes_edge = nodes_centroids[nodes_edge]
            tk_edge = tk_points[edge]
            faces_adj_edge = adjacencies[edge]
            h_dist_edge = h_dist[edge]

            q0_tk = np.linalg.norm(centroids_nodes_edge - tk_edge, axis=1)

            netas = [
                q0_tk[0]/h_dist_edge[0], 
                q0_tk[0]/h_dist_edge[1],
                q0_tk[1]/h_dist_edge[0],
                q0_tk[1]/h_dist_edge[1]
            ]

            all_neta.append(netas)
            all_vertice_id.append([nodes_edge[0], nodes_edge[0], nodes_edge[1], nodes_edge[1]])
            all_face_id.append([faces_adj_edge[0], faces_adj_edge[1], faces_adj_edge[0], faces_adj_edge[1]])
            all_edge_id.append(np.repeat(edge, 4))
        
        for edge in edges[bool_boundary_edges]:
            nodes_edge = nodes_of_edges[edge]
            centroids_nodes_edge = nodes_centroids[nodes_edge]
            tk_edge = tk_points[edge]
            face_adj_edge = adjacencies[edge, 0]
            h_dist_edge = h_dist[edge, 0]
            q0_tk = np.linalg.norm(centroids_nodes_edge - tk_edge, axis=1)

            netas = [
                q0_tk[0]/h_dist_edge,
                q0_tk[1]/h_dist_edge
            ]

            all_neta.append(netas)
            all_vertice_id.append(nodes_edge)
            all_face_id.append([face_adj_edge, face_adj_edge])
            all_edge_id.append(np.repeat(edge, 2))
        
        all_neta = np.concatenate(np.array(all_neta, dtype='O'))
        all_vertice_id = np.concatenate(np.array(all_vertice_id, dtype='O').flatten())
        all_face_id = np.concatenate(np.array(all_face_id, dtype='O'))
        all_edge_id = np.concatenate(np.array(all_edge_id, dtype='O'))

        n = len(all_neta)

        dtype = [('node_id', np.int), ('face_id', np.int), ('edge_id', np.int), ('neta', np.float64)]
        str_array = np.zeros(n, dtype=dtype)

        str_array['node_id'][:] = all_vertice_id
        str_array['face_id'][:] = all_face_id
        str_array['edge_id'][:] = all_edge_id
        str_array['neta'][:] = all_neta

        return {'neta': str_array}

    def cosin_law(self, a, b, c):
        """
        Retorna o angulo oposto ao lado c em radianos 
        do triangulo de lados consecutivos abc

        Args:
            a (_type_): comprimento do lado a
            b (_type_): comprimento do lado b
            c (_type_): comprimento do lado c
        Returns:
            float: angulo em radianos
        """

        cosC = -(c**2 - (a**2 + b**2))/(2*a*b)
        C = np.arccos(cosC)
        return C

    def create_knt_barra_vef(self, tk_points, nodes_centroids, permeability, edges_of_nodes, nodes_of_edges, nodes, edges, adjacencies, faces_of_nodes, bool_boundary_nodes, faces_centroids, **kwargs):

        R_matrix = self.get_Rmatrix()

        ## ids kn and kt (k_barra)
        all_kface_id = []
        all_kedge_id1 = []
        all_kedge_id2 = []
        all_knode_id = []
        all_kn = []
        all_kt = []

        ## ids of angle v (v angle is build by triangle: node[idv], edge[idv0], edge[idv1])
        all_node_idv = []
        all_edge_idv0 = []
        all_edge_idv1 = []
        all_v_angle = []

        ## ids of kn and kt (k normal)
        alln_node_id = []
        alln_edge_id = []
        alln_face_id = []
        alln_kn = []
        alln_kt = []
        all_theta = []
        all_phi = []


        inodes = ~bool_boundary_nodes
        
        for node in nodes[inodes]:
            faces_node = faces_of_nodes[node]
            edges_node = edges_of_nodes[node]
            local_edge_sort_index = np.repeat(-1, edges_node.max()+1)
            local_edge_sort_index[edges_node] = np.arange(len(edges_node))
            node_centroid = nodes_centroids[node]

            for face in faces_node:
                perm_face = permeability[face]
                face_centroid = faces_centroids[face]
                edges_selected = edges[(adjacencies[:, 0] == face) | (adjacencies[:, 1] == face)]
                edges_selected = np.intersect1d(edges_node, edges_selected)
                
                if (local_edge_sort_index[edges_selected[0]] == 0) & (local_edge_sort_index[edges_selected[1]] == local_edge_sort_index.max()):
                    edges_selected[:] = edges_selected[::-1]
                elif local_edge_sort_index[edges_selected[0]] > local_edge_sort_index[edges_selected[1]]:
                    edges_selected[:] = edges_selected[::-1]
                
                tk_points_edges_selected = tk_points[edges_selected]
                nodes_edges_selected = nodes_of_edges[edges_selected]
                tk_points_node_edges_selected = tk_points_edges_selected[nodes_edges_selected==node]
                
                tk_vector = tk_points_node_edges_selected[1] - tk_points_node_edges_selected[0]
                norm_tk = np.linalg.norm(tk_vector)
                v1 = R_matrix.dot(tk_vector).reshape((2,1))

                ## set kn and kt barra
                kn = v1.T.dot(perm_face).dot(v1).flatten()[0]/(norm_tk**2)
                kt = v1.T.dot(perm_face).dot(tk_vector).flatten()[0]/(norm_tk**2)
                
                all_kface_id.append(face)
                all_kedge_id1.append(edges_selected[0])
                all_kedge_id2.append(edges_selected[1])
                all_knode_id.append(node)
                all_kn.append(kn)
                all_kt.append(kt)

                ## set vangles
                q0_tk1 = tk_points_node_edges_selected[1] - node_centroid
                q0_tk0 = tk_points_node_edges_selected[0] - node_centroid

                q0_tk1n = np.linalg.norm(q0_tk1)
                q0_tk0n = np.linalg.norm(q0_tk0)

                vangle1 = self.cosin_law(norm_tk, q0_tk1n, q0_tk0n)
                vangle0 = self.cosin_law(q0_tk0n, norm_tk, q0_tk1n)
                
                all_node_idv.append(node)
                all_edge_idv0.append(edges_selected[0])
                all_edge_idv1.append(edges_selected[1])
                all_v_angle.append(vangle0)

                all_node_idv.append(node)
                all_edge_idv0.append(edges_selected[1])
                all_edge_idv1.append(edges_selected[0])
                all_v_angle.append(vangle1)

                ## set kn and kt (k normal), theta and phi angles
                v2 = R_matrix.dot(q0_tk1).reshape((2,1))
                v3 = R_matrix.dot(q0_tk0).reshape((2,1))

                knn_tk1 = v2.T.dot(perm_face).dot(v2).flatten()[0]/(q0_tk1n**2)
                ktn_tk1 = v2.T.dot(perm_face).dot(q0_tk1).flatten()[0]/(q0_tk1n**2)
                
                knn_tk0 = v3.T.dot(perm_face).dot(v3).flatten()[0]/(q0_tk0n**2)
                ktn_tk0 = v3.T.dot(perm_face).dot(q0_tk0).flatten()[0]/(q0_tk0n**2)

                # theta and phi
                q0_okn = np.linalg.norm(face_centroid - node_centroid)
                tk1_okn = np.linalg.norm(face_centroid - tk_points_node_edges_selected[1])
                tk0_okn = np.linalg.norm(face_centroid - tk_points_node_edges_selected[0])
                
                theta_tk0 = self.cosin_law(q0_okn, q0_tk0n, tk0_okn)
                phi_tk0 = self.cosin_law(tk0_okn, q0_okn, q0_tk0n)
                theta_tk1 = self.cosin_law(q0_tk1n, q0_okn, tk1_okn)
                phi_tk1 = self.cosin_law(q0_okn, tk1_okn, q0_tk1n)

                alln_node_id.append(node)
                alln_edge_id.append(edges_selected[1])
                alln_face_id.append(face)
                alln_kn.append(knn_tk1)
                alln_kt.append(ktn_tk1)
                all_theta.append(theta_tk1)                
                all_phi.append(phi_tk1)                

                alln_node_id.append(node)
                alln_edge_id.append(edges_selected[0])
                alln_face_id.append(face)
                alln_kn.append(knn_tk0)
                alln_kt.append(ktn_tk0)
                all_theta.append(theta_tk0)                
                all_phi.append(phi_tk0)




        for node in nodes[bool_boundary_nodes]:
            faces_node = faces_of_nodes[node]
            edges_node = edges_of_nodes[node]
            local_edge_sort_index = np.repeat(-1, edges_node.max()+1)
            local_edge_sort_index[edges_node] = np.arange(len(edges_node))
            node_centroid = nodes_centroids[node]

            for face in faces_node:
                perm_face = permeability[face]
                edges_selected = edges[(adjacencies[:, 0] == face) | (adjacencies[:, 1] == face)]
                edges_selected = np.intersect1d(edges_node, edges_selected)
                
                if local_edge_sort_index[edges_selected[0]] > local_edge_sort_index[edges_selected[1]]:
                    edges_selected[:] = edges_selected[::-1]
                
                tk_points_edges_selected = tk_points[edges_selected]
                nodes_edges_selected = nodes_of_edges[edges_selected]
                tk_points_node_edges_selected = tk_points_edges_selected[nodes_edges_selected==node]
                
                tk_vector = tk_points_node_edges_selected[1] - tk_points_node_edges_selected[0]
                norm_tk = np.linalg.norm(tk_vector)
                v1 = R_matrix.dot(tk_vector).reshape((2,1))

                ## set kn and kt
                kn = v1.T.dot(perm_face).dot(v1).flatten()[0]/(norm_tk**2)
                kt = v1.T.dot(perm_face).dot(tk_vector).flatten()[0]/(norm_tk**2)
                
                all_kface_id.append(face)
                all_kedge_id1.append(edges_selected[0])
                all_kedge_id2.append(edges_selected[1])
                all_knode_id.append(node)
                all_kn.append(kn)
                all_kt.append(kt)

                ## set vangles
                q0_tk1 = tk_points_node_edges_selected[1] - node_centroid
                q0_tk0 = tk_points_node_edges_selected[0] - node_centroid

                q0_tk1n = np.linalg.norm(q0_tk1)
                q0_tk0n = np.linalg.norm(q0_tk0)

                vangle1 = self.cosin_law(norm_tk, q0_tk1n, q0_tk0n)
                vangle0 = self.cosin_law(q0_tk0n, norm_tk, q0_tk1n)
                
                all_node_idv.append(node)
                all_edge_idv0.append(edges_selected[0])
                all_edge_idv1.append(edges_selected[1])
                all_v_angle.append(vangle0)

                all_node_idv.append(node)
                all_edge_idv0.append(edges_selected[1])
                all_edge_idv1.append(edges_selected[0])
                all_v_angle.append(vangle1)

                ## set kn and kt (k normal), theta and phi angles
                v2 = R_matrix.dot(q0_tk1).reshape((2,1))
                v3 = R_matrix.dot(q0_tk0).reshape((2,1))

                knn_tk1 = v2.T.dot(perm_face).dot(v2).flatten()[0]/(q0_tk1n**2)
                ktn_tk1 = v2.T.dot(perm_face).dot(q0_tk1).flatten()[0]/(q0_tk1n**2)
                
                knn_tk0 = v3.T.dot(perm_face).dot(v3).flatten()[0]/(q0_tk0n**2)
                ktn_tk0 = v3.T.dot(perm_face).dot(q0_tk0).flatten()[0]/(q0_tk0n**2)

                # theta and phi
                q0_okn = np.linalg.norm(face_centroid - node_centroid)
                tk1_okn = np.linalg.norm(face_centroid - tk_points_node_edges_selected[1])
                tk0_okn = np.linalg.norm(face_centroid - tk_points_node_edges_selected[0])
                
                theta_tk0 = self.cosin_law(q0_okn, q0_tk0n, tk0_okn)
                phi_tk0 = self.cosin_law(tk0_okn, q0_okn, q0_tk0n)
                theta_tk1 = self.cosin_law(q0_tk1n, q0_okn, tk1_okn)
                phi_tk1 = self.cosin_law(q0_okn, tk1_okn, q0_tk1n)

                alln_node_id.append(node)
                alln_edge_id.append(edges_selected[1])
                alln_face_id.append(face)
                alln_kn.append(knn_tk1)
                alln_kt.append(ktn_tk1)
                all_theta.append(theta_tk1)                
                all_phi.append(phi_tk1)                

                alln_node_id.append(node)
                alln_edge_id.append(edges_selected[0])
                alln_face_id.append(face)
                alln_kn.append(knn_tk0)
                alln_kt.append(ktn_tk0)
                all_theta.append(theta_tk0)                
                all_phi.append(phi_tk0)


        all_kface_id = np.array(all_kface_id)
        all_kedge_id1 = np.array(all_kedge_id1)
        all_kedge_id2 = np.array(all_kedge_id2)
        all_knode_id = np.array(all_knode_id)
        all_kn = np.array(all_kn)
        all_kt = np.array(all_kt)
        dtype1 = [
            ('face_id', np.int),
            ('edge_id0', np.int),
            ('edge_id1', np.int),
            ('node_id', np.int),
            ('kn', np.float64),
            ('kt', np.float64)
        ]
        array1 = np.zeros(len(all_kface_id), dtype=dtype1)
        array1['face_id'][:] = all_kface_id
        array1['edge_id0'][:] = all_kedge_id1
        array1['edge_id1'][:] = all_kedge_id2
        array1['node_id'][:] = all_knode_id
        array1['kn'][:] = all_kn
        array1['kt'][:] = all_kt

        all_node_idv = np.array(all_node_idv)
        all_edge_idv0 = np.array(all_edge_idv0)
        all_edge_idv1 = np.array(all_edge_idv1)
        all_v_angle = np.array(all_v_angle)
        dtype2 = [
            ('edge_id0', np.int),
            ('edge_id1', np.int),
            ('node_id', np.int),
            ('v_angle', np.float64)
        ]
        array2 = np.zeros(len(all_node_idv), dtype=dtype2)
        array2['node_id'][:] = all_node_idv
        array2['edge_id1'][:] = all_edge_idv1
        array2['edge_id0'][:] = all_edge_idv0
        array2['v_angle'][:] = all_v_angle

        alln_node_id = np.array(alln_node_id)
        alln_edge_id = np.array(alln_edge_id)
        alln_face_id = np.array(alln_face_id)
        alln_kn = np.array(alln_kn)
        alln_kt = np.array(alln_kt)
        all_theta = np.array(all_theta)
        all_phi = np.array(all_phi)
        dtype3 = [
            ('edge_id', np.int),
            ('face_id', np.int),
            ('node_id', np.int),
            ('kn', np.float64),
            ('kt', np.float64),
            ('theta', np.float64),
            ('phi', np.float64)
        ]
        array3 = np.zeros(len(alln_node_id), dtype=dtype3)
        array3['node_id'][:] = alln_node_id
        array3['edge_id'][:] = alln_edge_id
        array3['face_id'][:] = alln_face_id
        array3['kn'][:] = alln_kn
        array3['kt'][:] = alln_kt
        array3['theta'][:] = all_theta
        array3['phi'][:] = all_phi

        resp = {
            'kn_kt_barra': array1,
            'v_angle': array2,
            'kn_kt_theta_phi': array3
        }

        return resp

    def create_zeta(self, kn_kt_barra, v_angle, kn_kt_theta_phi, adjacencies, edges, nodes_of_edges, bool_boundary_edges, neta, bool_boundary_nodes, nodes, edges_of_nodes, **kwargs):
        
        terms = np.zeros(12)
        l_terms = np.array([3, 4, 6, 9, 10, 12]) - 1
        r_terms = np.array([1, 2, 5, 7, 8, 11]) - 1
        
        all_zeta = []
        all_node_id = []
        all_edge_id = []
        
        biedges = ~bool_boundary_edges
        binodes = ~bool_boundary_nodes
        
        bound_edges = edges[bool_boundary_edges]
        
        for node in nodes:
            edges_node = edges_of_nodes[node]
            local_edge_sort_index = np.arange(len(edges_node))
            # local_edge_sort_index[edges_node] = np.arange(len(edges_node))
            
            local_bedges = np.intersect1d(edges_node, bound_edges)
            local_iedges = np.setdiff1d(edges_node, local_bedges)
            
            for edge in local_bedges:
                pass
            
            for edge in local_iedges:
                index_edge = local_edge_sort_index[edges_node==edge][0]
                edge0 = edges_node[local_edge_sort_index == index_edge-1]
                edge1 = edges_node[local_edge_sort_index == index_edge+1]
                
                
            
                import pdb; pdb.set_trace()
            
            
            
            
            
            
            
            
            






