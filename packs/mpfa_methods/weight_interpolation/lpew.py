import numpy as np
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
import pandas as pd
from packs.manager.meshmanager import MeshProperty
from packs.utils import calculate_face_properties

class LpewWeight:
    """A LINEARITY PRESERVING CELL-CENTERED
        SCHEME FOR THE HETEROGENEOUS AND
        ANISOTROPIC DIFFUSION EQUATIONS ON
        GENERAL MESHES

    """

    datas = ['tk_points', 'neta', 'kn_kt_barra', 'kn_kt_theta_phi_vangle',
             'zeta', 'lambda_barra']
    
    data_weights = ['nodes_weights']

    @property
    def Ok_values(self):
        return ['k-1', 'k']

    @staticmethod
    def preprocess(
        nodes_centroids,
        unitary_normal_edges,
        nodes_of_edges,
        edges,
        nodes,
        nodes_of_nodes,
        edges_of_nodes,
        faces_centroids,
        faces_of_nodes,
        **kwargs
    ):
        
        lsds = LsdsFluxCalculation()

        return lsds.preprocess(
            nodes_centroids,
            unitary_normal_edges,
            nodes_of_edges,
            edges,
            nodes,
            nodes_of_nodes,
            edges_of_nodes,
            faces_centroids,
            faces_of_nodes
        )
    
    def cot(self, theta):
        """get cotangent

        Args:
            theta (_type_): angle
        """

        return (np.cos(theta)/np.sin(theta))

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

    def insert_knt_barra_vef_data(
            self,
            R_matrix,
            tk_points,
            edges_selected,
            nodes_of_edges,
            node,
            face,
            perm_face,
            node_centroid,
            face_centroid,
            all_kface_id,
            all_kedge_id1,
            all_kedge_id2,
            all_knode_id,
            all_kn,
            all_kt,
            alln_node_id,
            alln_face_id,
            alln_edge_id,
            alln_kn,
            alln_kt,
            all_theta,
            all_phi,
            all_v_angle,
            **kwargs            
    ):

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

        # vangle1 = np.arccos(
        #     (np.dot(-q0_tk1, -tk0_tk1))/(q0_tk1n*norm_tk)
        # )
        # vangle0 = np.arccos(
        #     (np.dot(-q0_tk0, tk0_tk1))/(q0_tk0n*norm_tk)
        # )

        ## set kn and kt (k normal), theta and phi angles
        v2 = R_matrix.dot(q0_tk1).reshape((2,1))
        v3 = R_matrix.dot(q0_tk0).reshape((2,1))

        knn_tk1 = v2.T.dot(perm_face).dot(v2).flatten()[0]/(q0_tk1n**2)
        ktn_tk1 = v2.T.dot(perm_face).dot(q0_tk1).flatten()[0]/(q0_tk1n**2)
        
        knn_tk0 = v3.T.dot(perm_face).dot(v3).flatten()[0]/(q0_tk0n**2)
        ktn_tk0 = v3.T.dot(perm_face).dot(q0_tk0).flatten()[0]/(q0_tk0n**2)

        # theta and phi
        q0_ok = face_centroid - node_centroid
        tk1_ok = face_centroid - tk_points_node_edges_selected[1]
        tk0_ok = face_centroid - tk_points_node_edges_selected[0]
        q0_okn = np.linalg.norm(q0_ok)
        tk1_okn = np.linalg.norm(tk1_ok)
        tk0_okn = np.linalg.norm(tk0_ok)
        
        theta_tk0 = self.cosin_law(q0_okn, q0_tk0n, tk0_okn)
        phi_tk0 = self.cosin_law(tk0_okn, q0_okn, q0_tk0n)
        theta_tk1 = self.cosin_law(q0_tk1n, q0_okn, tk1_okn)
        phi_tk1 = self.cosin_law(q0_okn, tk1_okn, q0_tk1n)

        # theta_tk02 = np.arccos(
        #     (np.dot(q0_ok, q0_tk0))/(q0_okn*q0_tk0n)
        # )
        # phi_tk02 = np.arccos(
        #     (np.dot(-q0_ok, -tk0_ok))/(q0_okn*tk0_okn)
        # )
        # theta_tk12 = np.arccos(
        #     (np.dot(q0_ok, q0_tk1))/(q0_okn*q0_tk1n)
        # )
        # phi_tk12 = np.arccos(
        #     (np.dot(-q0_ok, -tk1_ok))/(q0_okn*tk1_okn)
        # )

        alln_node_id.append(node)
        alln_edge_id.append(edges_selected[1])
        alln_face_id.append(face)
        alln_kn.append(knn_tk1)
        alln_kt.append(ktn_tk1)
        all_theta.append(theta_tk1)                
        all_phi.append(phi_tk1)
        all_v_angle.append(vangle1)             

        alln_node_id.append(node)
        alln_edge_id.append(edges_selected[0])
        alln_face_id.append(face)
        alln_kn.append(knn_tk0)
        alln_kt.append(ktn_tk0)
        all_theta.append(theta_tk0)                
        all_phi.append(phi_tk0)
        all_v_angle.append(vangle0)

    def create_knt_barra_vef(self, tk_points, nodes_centroids, permeability, edges_of_nodes, nodes_of_edges, nodes, edges, adjacencies, faces_of_nodes, bool_boundary_nodes, faces_centroids, **kwargs):

        R_matrix = self.get_Rmatrix()

        ## ids kn and kt (k_barra)
        all_kface_id = []
        all_kedge_id1 = []
        all_kedge_id2 = []
        all_knode_id = []
        all_kn = []
        all_kt = []

        ## ids of kn and kt (k normal), theta, phi, v_angle
        alln_node_id = []
        alln_edge_id = []
        alln_face_id = []
        alln_kn = []
        alln_kt = []
        all_theta = []
        all_phi = []
        all_v_angle = []

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
                edges_selected = (edges[(adjacencies[:, 0] == face) | (adjacencies[:, 1] == face)]).copy()
                edges_selected = np.intersect1d(edges_node, edges_selected)
                
                if (local_edge_sort_index[edges_selected[0]] == 0) & (local_edge_sort_index[edges_selected[1]] == local_edge_sort_index.max()):
                    edges_selected[:] = edges_selected[::-1]
                elif local_edge_sort_index[edges_selected[0]] > local_edge_sort_index[edges_selected[1]]:
                    edges_selected[:] = edges_selected[::-1]

                self.insert_knt_barra_vef_data(
                    R_matrix,
                    tk_points,
                    edges_selected,
                    nodes_of_edges,
                    node,
                    face,
                    perm_face,
                    node_centroid,
                    face_centroid,
                    all_kface_id,
                    all_kedge_id1,
                    all_kedge_id2,
                    all_knode_id,
                    all_kn,
                    all_kt,
                    alln_node_id,
                    alln_face_id,
                    alln_edge_id,
                    alln_kn,
                    alln_kt,
                    all_theta,
                    all_phi,
                    all_v_angle
                )

        for node in nodes[bool_boundary_nodes]:
            faces_node = faces_of_nodes[node]
            edges_node = edges_of_nodes[node]
            local_edge_sort_index = np.repeat(-1, edges_node.max()+1)
            local_edge_sort_index[edges_node] = np.arange(len(edges_node))
            node_centroid = nodes_centroids[node]

            for face in faces_node:
                perm_face = permeability[face]
                face_centroid = faces_centroids[face]
                edges_selected = (edges[(adjacencies[:, 0] == face) | (adjacencies[:, 1] == face)]).copy()
                edges_selected = np.intersect1d(edges_node, edges_selected)
                
                if local_edge_sort_index[edges_selected[0]] > local_edge_sort_index[edges_selected[1]]:
                    edges_selected[:] = edges_selected[::-1]
                
                self.insert_knt_barra_vef_data(
                    R_matrix,
                    tk_points,
                    edges_selected,
                    nodes_of_edges,
                    node,
                    face,
                    perm_face,
                    node_centroid,
                    face_centroid,
                    all_kface_id,
                    all_kedge_id1,
                    all_kedge_id2,
                    all_knode_id,
                    all_kn,
                    all_kt,
                    alln_node_id,
                    alln_face_id,
                    alln_edge_id,
                    alln_kn,
                    alln_kt,
                    all_theta,
                    all_phi,
                    all_v_angle
                )                

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

        alln_node_id = np.array(alln_node_id)
        alln_edge_id = np.array(alln_edge_id)
        alln_face_id = np.array(alln_face_id)
        alln_kn = np.array(alln_kn)
        alln_kt = np.array(alln_kt)
        all_theta = np.array(all_theta)
        all_phi = np.array(all_phi)
        all_v_angle = np.array(all_v_angle)
        dtype3 = [
            ('edge_id', np.int),
            ('face_id', np.int),
            ('node_id', np.int),
            ('kn', np.float64),
            ('kt', np.float64),
            ('theta', np.float64),
            ('phi', np.float64),
            ('v_angle', np.float64)
        ]
        array3 = np.zeros(len(alln_node_id), dtype=dtype3)
        array3['node_id'][:] = alln_node_id
        array3['edge_id'][:] = alln_edge_id
        array3['face_id'][:] = alln_face_id
        array3['kn'][:] = alln_kn
        array3['kt'][:] = alln_kt
        array3['theta'][:] = all_theta
        array3['phi'][:] = all_phi
        array3['v_angle'][:] = all_v_angle

        resp = {
            'kn_kt_barra': array1,
            'kn_kt_theta_phi_vangle': array3
        }

        return resp

    def get_zeta_from_terms(self, terms):
        t1 = terms[0]*terms[1] + terms[2]*terms[3] + terms[4] - terms[5]
        t2 = terms[6]*terms[7] + terms[8]*terms[9] - terms[10] + terms[11]
        return t1/t2

    def get_Ok_value_by_edge(self, reference_edge, reference_node, edge_face, nodes_of_edges, nodes_centroids, **kwargs):
        """
            verificar se o edge da face esta a esquerda ou direita do edge de referencia
        """
        z_direction = np.array([0, 0, 1])
        nodes_edge_reference = (nodes_of_edges[reference_edge]).copy()
        nodes_edge_face = (nodes_of_edges[edge_face]).copy()

        if nodes_edge_reference[1] == reference_node:
            nodes_edge_reference[:] = nodes_edge_reference[::-1]
        
        if nodes_edge_face[1] == reference_node:
            nodes_edge_face[:] = nodes_edge_face[::-1]

        c1 = nodes_centroids[nodes_edge_reference]
        v1 = np.zeros(3)
        v1[0:2] = c1[1] - c1[0]
        c2 = nodes_centroids[nodes_edge_face]
        v2 = v1.copy()
        v2[0:2] = c2[1] - c2[0]

        v3 = np.cross(v1, v2)
        test = np.dot(v3, z_direction)

        if test > 0:
            return self.Ok_values[1]
        elif test < 0:
            return self.Ok_values[0]
        else:
            raise ValueError

    def get_zeta_face_terms(self, kn_kt_barra, kn_kt_theta_phi_vangle, node, face, edge_reference, edge_face, **kwargs):
        terms = np.zeros(6)
        try:
            terms[0] = kn_kt_barra['kn'][
                (kn_kt_barra['node_id']==node) & 
                (kn_kt_barra['face_id']==face)
            ]
        except:
            import pdb; pdb.set_trace()

        terms[1] = self.cot(kn_kt_theta_phi_vangle['v_angle'][
            (kn_kt_theta_phi_vangle['node_id']==node) & 
            (kn_kt_theta_phi_vangle['face_id']==face) &
            (kn_kt_theta_phi_vangle['edge_id']==edge_face)
        ])

        terms[2] = kn_kt_barra['kt'][
            (kn_kt_barra['node_id']==node) & 
            (kn_kt_barra['face_id']==face)
        ]
        
        terms[3] = kn_kt_theta_phi_vangle['kn'][
            (kn_kt_theta_phi_vangle['node_id']==node) & 
            (kn_kt_theta_phi_vangle['face_id']==face) &
            (kn_kt_theta_phi_vangle['edge_id']==edge_reference)
        ]

        terms[4] = self.cot(kn_kt_theta_phi_vangle['theta'][
            (kn_kt_theta_phi_vangle['node_id']==node) & 
            (kn_kt_theta_phi_vangle['face_id']==face) &
            (kn_kt_theta_phi_vangle['edge_id']==edge_reference)
        ])

        terms[5] = kn_kt_theta_phi_vangle['kt'][
            (kn_kt_theta_phi_vangle['node_id']==node) & 
            (kn_kt_theta_phi_vangle['face_id']==face) &
            (kn_kt_theta_phi_vangle['edge_id']==edge_reference)
        ]
        
        return terms

    def create_zeta(self, kn_kt_barra, kn_kt_theta_phi_vangle, adjacencies, edges, bool_boundary_edges, nodes, edges_of_nodes, nodes_centroids, nodes_of_edges, **kwargs):
        
        terms = np.zeros(12)
        l_terms = np.array([3, 4, 6, 9, 10, 12]) - 1 # left terms
        r_terms = np.array([1, 2, 5, 7, 8, 11]) - 1 # right terms
        
        all_zeta = []
        all_node_id = []
        all_edge_id = []
        
        bound_edges = edges[bool_boundary_edges]
        for node in nodes:
            
            edges_node = edges_of_nodes[node]
            local_edge_sort_index = np.arange(len(edges_node))
            
            local_bedges = np.intersect1d(edges_node, bound_edges)
            local_iedges = np.setdiff1d(edges_node, local_bedges)

            for edge in local_bedges:
                terms[:] = 0
                index_edge = local_edge_sort_index[edges_node==edge][0]
                faces_adj = adjacencies[edge]

                for face in faces_adj:
                    if face == -1:
                        continue
                    edges_face = edges[(adjacencies[:, 0] == face) | (adjacencies[:, 1] == face)]
                    edge_face = np.intersect1d(edges_node, edges_face)
                    edge_face = edge_face[edge_face != edge][0]                   

                    Ok = self.get_Ok_value_by_edge(
                        edge,
                        node,
                        edge_face,
                        nodes_of_edges,
                        nodes_centroids
                    )

                    if Ok == self.Ok_values[0]: ## face k
                        terms[r_terms] = self.get_zeta_face_terms(
                            kn_kt_barra,
                            kn_kt_theta_phi_vangle,
                            node,
                            face,
                            edge,
                            edge_face
                        )
                    
                    elif Ok == self.Ok_values[1]: ## face k+1
                        terms[l_terms] = self.get_zeta_face_terms(
                            kn_kt_barra,
                            kn_kt_theta_phi_vangle,
                            node,
                            face,
                            edge,
                            edge_face
                        )
                    else:
                        raise ValueError
                
                all_zeta.append(self.get_zeta_from_terms(terms))
                all_edge_id.append(edge)
                all_node_id.append(node)
            
            for edge in local_iedges:
                terms[:] = 0
                
                faces_adj = adjacencies[edge]

                edges_face0 = edges[(adjacencies[:, 0] == faces_adj[0]) | (adjacencies[:, 1] == faces_adj[0])]
                edges_face1 = edges[(adjacencies[:, 0] == faces_adj[1]) | (adjacencies[:, 1] == faces_adj[1])]

                edge_face0 = np.intersect1d(edges_node, edges_face0)
                edge_face0 = edge_face0[edge_face0 != edge][0]
                edge_face1 = np.intersect1d(edges_node, edges_face1)
                edge_face1 = edge_face1[edge_face1 != edge][0]
                edges_faces = [edge_face0, edge_face1]
                
                index_edge = local_edge_sort_index[edges_node==edge][0]
                edge0 = edges_node[local_edge_sort_index[index_edge-1]]
                try:
                    edge1 = edges_node[local_edge_sort_index[index_edge+1]]
                except IndexError:
                    edge1 = edges_node[0]
                
                if len(faces_adj[edges_faces==edge1]) == 0:
                    import pdb; pdb.set_trace()
                elif len(faces_adj[edges_faces==edge0]) == 0:
                    import pdb; pdb.set_trace()
                
            

                terms[l_terms] = self.get_zeta_face_terms(
                    kn_kt_barra,
                    kn_kt_theta_phi_vangle,
                    node,
                    faces_adj[edges_faces==edge1],
                    edge,
                    edge1
                )
                
                terms[r_terms] = self.get_zeta_face_terms(
                    kn_kt_barra,
                    kn_kt_theta_phi_vangle,
                    node,
                    faces_adj[edges_faces==edge0],
                    edge,
                    edge0
                )
                
                all_zeta.append(self.get_zeta_from_terms(terms))
                all_edge_id.append(edge)
                all_node_id.append(node)
        
        all_zeta = np.array(all_zeta)
        all_node_id = np.array(all_node_id)
        all_edge_id = np.array(all_edge_id)

        dtype = [
            ('node_id', np.int),
            ('edge_id', np.int),
            ('zeta', np.float64)
        ]

        array = np.zeros(len(all_zeta), dtype=dtype)
        array['zeta'] = all_zeta
        array['node_id'] = all_node_id
        array['edge_id'] = all_edge_id

        return {'zeta': array}
    
    def get_lambda_barra_terms(self, kn_kt_theta_phi_vangle, neta, zeta, node, face, edge):
        
        kn1 = kn_kt_theta_phi_vangle['kn'][
            (kn_kt_theta_phi_vangle['node_id']==node) &
            (kn_kt_theta_phi_vangle['edge_id']==edge) &
            (kn_kt_theta_phi_vangle['face_id']==face)
        ]

        neta1 = neta['neta'][
            (neta['node_id']==node) &
            (neta['face_id']==face) &
            (neta['edge_id']==edge)
        ]
        
        zeta_edge1 = zeta['zeta'][
            (zeta['node_id']==node) &
            (zeta['edge_id']==edge)
        ]
        
        return kn1*neta1*zeta_edge1
    
    def create_lambda_barra(self, kn_kt_theta_phi_vangle, neta, zeta, nodes, faces_of_nodes, adjacencies, edges_of_nodes, edges, **kwargs):

        all_lambda = []
        all_node_id = []
        all_face_id = []

        for node in nodes:
            faces_node = faces_of_nodes[node]
            edges_node = edges_of_nodes[node]
            for face in faces_node:
                edges_face = edges[(adjacencies[:, 0] == face) | (adjacencies[:, 1] == face)]
                edges_face = np.intersect1d(edges_face, edges_node)
                
                term1 = self.get_lambda_barra_terms(
                    kn_kt_theta_phi_vangle,
                    neta,
                    zeta,
                    node,
                    face,
                    edges_face[0]
                )
                
                term2 = self.get_lambda_barra_terms(
                    kn_kt_theta_phi_vangle,
                    neta,
                    zeta,
                    node,
                    face,
                    edges_face[1]
                )
                
                all_lambda.append(term1+term2)
                all_node_id.append(node)
                all_face_id.append(face)
        
        all_lambda = np.array(all_lambda).flatten()
        all_node_id = np.array(all_node_id)
        all_face_id = np.array(all_face_id)
        
        dtype = [
            ('node_id', np.int),
            ('face_id', np.int),
            ('lambda_barra', np.float64)
        ]
        
        array = np.zeros(len(all_lambda), dtype=dtype)
        array['node_id'] = all_node_id
        array['face_id'] = all_face_id
        array['lambda_barra'] = all_lambda
        
        return {'lambda_barra': array}

    def create_lpew2_weights(self, lambda_barra: np.ndarray, **kwargs):

        df = pd.DataFrame({
            'node_id': lambda_barra['node_id'],
            'face_id': lambda_barra['face_id'],
            'lambda_barra': lambda_barra['lambda_barra']
        })

        df2 = df.groupby(['node_id']).sum()['lambda_barra']
        
        df['weight'] = df.apply(lambda row: row['lambda_barra']/df2[df2.index==row['node_id']].values[0],axis=1)
    
        dtype = [
            ('node_id', np.int),
            ('face_id', np.int),
            ('weight', np.float64)
        ]
        array = np.zeros(df.shape[0], dtype=dtype)
        array['node_id'] = df['node_id'].values
        array['face_id'] = df['face_id'].values
        array['weight'] = df['weight'].values

        return {'nodes_weights': array}


def intermediate(mesh_properties: MeshProperty):
    mesh_properties.export_data()

def preprocess(mesh_properties: MeshProperty):

    lpew2 = LpewWeight()
    resp = lpew2.preprocess(**mesh_properties.get_all_data())
    mesh_properties.insert_or_update_data(resp)
    intermediate(mesh_properties)

def create_Tk(mesh_properties: MeshProperty):
    """criar k normal e tangente
    """
    k = 0
    lpew = LpewWeight()
    if mesh_properties.verify_name_in_data_names(lpew.datas[k]):
        return   
    resp = lpew.create_Tk_points(
        mesh_properties['nodes_centroids'],
        mesh_properties['nodes_of_edges']
    )
    mesh_properties.insert_data(resp)
    intermediate(mesh_properties)

def create_neta(mesh_properties: MeshProperty):
    
    k = 1
    lpew = LpewWeight()
    if mesh_properties.verify_name_in_data_names(lpew.datas[k]):
        return 
    resp = lpew.create_neta(**mesh_properties.get_all_data())
    mesh_properties.insert_data(resp)
    intermediate(mesh_properties)

def create_knt_vef(mesh_properties: MeshProperty):
    
    k = 2
    lpew = LpewWeight()
    if mesh_properties.verify_name_in_data_names(lpew.datas[k]):
        return 
    resp = lpew.create_knt_barra_vef(**mesh_properties.get_all_data())
    mesh_properties.insert_data(resp)
    intermediate(mesh_properties)

def create_zeta(mesh_properties: MeshProperty):
    
    k = 4
    lpew = LpewWeight()
    if mesh_properties.verify_name_in_data_names(lpew.datas[k]):
        return 
    resp = lpew.create_zeta(**mesh_properties.get_all_data())
    mesh_properties.insert_data(resp)
    intermediate(mesh_properties)

def create_lambda_barra(mesh_properties: MeshProperty):
    
    k = 5
    lpew = LpewWeight()
    if mesh_properties.verify_name_in_data_names(lpew.datas[k]):
        return 
    resp = lpew.create_lambda_barra(**mesh_properties.get_all_data())
    mesh_properties.insert_data(resp)
    intermediate(mesh_properties)

def create_lpew2_weights(mesh_properties: MeshProperty):
    
    k = 0
    lpew = LpewWeight()
    if mesh_properties.verify_name_in_data_names(lpew.data_weights[k]):
        return 
    resp = lpew.create_lpew2_weights(**mesh_properties.get_all_data())
    mesh_properties.insert_data(resp)
    intermediate(mesh_properties)

def get_lpew2_weights(mesh_properties: MeshProperty, **kwargs):

    create_Tk(mesh_properties)
    create_neta(mesh_properties)
    create_knt_vef(mesh_properties)
    create_zeta(mesh_properties)
    create_lambda_barra(mesh_properties)
    create_lpew2_weights(mesh_properties)




    


