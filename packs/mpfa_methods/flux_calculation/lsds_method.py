import numpy as np
import scipy.sparse as sp

class LsdsFluxCalculation:
    """
        Calulate least squares mpfad flux

        paper: A least squares based diamond scheme for anisotropic
                diffusion problems on polygonal meshes
                
                doi: 10.1002/fld.5031
    
    """
    def preprocess(
            self,
            nodes_centroids,
            unitary_normal_edges,
            nodes_of_edges,
            edges,
            **kwargs
    ):
        resp = dict()

        resp.update(
            self.define_A_B_points_of_edges(
                nodes_centroids,
                unitary_normal_edges,
                nodes_of_edges,
                edges
            )
        )

        return resp
        
    def define_A_B_points_of_edges(
            self,
            nodes_centroids,
            unitary_normal_edges,
            nodes_of_edges,
            edges,
            **kwargs
    ):
        """define os pontos A e B de cada edge da malha
            o ponto B deve esta a esquerda do vetor normal do edge
            e o ponto A deve estar a direita
            de modo que:
                nodes_of_edges[edge, 0] = B
                nodes_of_edges[edge, 1] = A

        Args:
            nodes_centroids (_type_): _description_
            unitary_normal_edges (_type_): _description_
            nodes_of_edges (_type_): _description_
            edges (_type_): _description_

        Returns:
            nodes_of_edges reorganizado
        """
        
        resp = nodes_of_edges.copy()

        theta = np.pi/2
        cossin = np.array([np.cos(theta), np.sin(theta)])
        R = np.array([
            [cossin[0], -cossin[1]],
            [cossin[1],  cossin[0]]
        ])

        for edge in edges:
            B = nodes_of_edges[edge, 0]
            A = nodes_of_edges[edge, 1]
            AB = nodes_centroids[B] - nodes_centroids[A]
            unitary_normal_vector_rotated = R.dot(unitary_normal_edges[edge])
            proj = AB.dot(unitary_normal_vector_rotated)
            if proj > 0:
                pass
            else:
                resp[edge,:] = [A, B]
        
        return {'nodes_of_edges': resp}

    def get_Skl(
            self,
            edges,
            bool_boundary_edges,
            nodes_of_edges,
            adjacencies,
            unitary_normal_edges,
            nodes_centroids,
            permeability,
            **kwargs
    ):
        
        bool_internal_edges = ~bool_boundary_edges
        S_alpha_sigma = self.get_local_S_alpha_sigma()
        Skl = np.zeros((edges.shape[0], 3, 3))

        for edge in edges[bool_internal_edges]:
            edge_nodes = nodes_of_edges[edge]
            faces_adj = adjacencies[edge]
            Skl[edge] = self.get_local_Skl(
                nodes_centroids[edge_nodes],
                unitary_normal_edges[edge],
                permeability[faces_adj],
                S_alpha_sigma
            )

        return Skl

    def get_local_Skl(
            self,
            edge_nodes_centroids,
            unitary_normal_edge,
            faces_permeabilities,
            S_alpha_sigma,
            **kwargs
    ):
        
        # nk_k = np.dot(unitary_normal_edge, faces_permeabilities[0])
        # nk_l = np.dot(unitary_normal_edge, faces_permeabilities[1])
        # A_centroid = edge_nodes_centroids[1]
        # B_centroid = edge_nodes_centroids[0]

        # S_alpha_sigma[:, 0, 0:2] = A_centroid
        # S_alpha_sigma[:, 1, 0:2] = B_centroid
        # S_alpha_sigma[0, 2, 0:2] = nk_k
        # S_alpha_sigma[1, 2, 0:2] = nk_l
        
        self.update_local_S_alpha_sigma(
            S_alpha_sigma,
            unitary_normal_edge,
            faces_permeabilities,
            edge_nodes_centroids
        )

        Skl = np.linalg.inv(S_alpha_sigma[1]).dot(S_alpha_sigma[0])
        
        return Skl
    
    def update_local_S_alpha_sigma(
            self,
            S_alpha_sigma,
            unitary_normal_edge,
            faces_permeabilities,
            edge_nodes_centroids
    ):
        
        nk_k = np.dot(unitary_normal_edge, faces_permeabilities[0])
        nk_l = np.dot(unitary_normal_edge, faces_permeabilities[1])
        A_centroid = edge_nodes_centroids[1]
        B_centroid = edge_nodes_centroids[0]

        S_alpha_sigma[:, 0, 0:2] = A_centroid
        S_alpha_sigma[:, 1, 0:2] = B_centroid
        S_alpha_sigma[0, 2, 0:2] = nk_k
        S_alpha_sigma[1, 2, 0:2] = nk_l
        
    def get_local_S_alpha_sigma(self):

        S_alpha_sigma = np.zeros((2, 3, 3))
        S_alpha_sigma[:, [0,1], [2]] = 1
        return S_alpha_sigma

    def get_x_and_y_k_sigma(
            self,
            edges,
            edges_dim,
            unitary_normal_edges,
            permeability,
            adjacencies,
            **kwargs
    ):
        
        x_and_y_k_sigma = np.zeros(adjacencies.shape)

        for edge in edges:
            edge_dim = edges_dim[edge]
            unitary_normal_edge = unitary_normal_edges[edge]
            face_k = adjacencies[edge][0]
            permeability_face_k = permeability[face_k]
            x_and_y_k_sigma[edge] = -edge_dim*np.dot(unitary_normal_edge, permeability_face_k)
        
        return x_and_y_k_sigma

    def get_D_and_mi(
            self,
            adjacencies,
            nodes_centroids,
            faces_centroids,
            nodes_of_edges,
            edges,
            **kwargs
    ):
        xk = faces_centroids[:, 0][adjacencies[:, 0]]
        yk = faces_centroids[:, 1][adjacencies[:, 0]]
        xl = faces_centroids[:, 0][adjacencies[:, 1]]
        yl = faces_centroids[:, 1][adjacencies[:, 1]]
        xA = nodes_centroids[:, 0][nodes_of_edges[edges][:,1]]
        yA = nodes_centroids[:, 1][nodes_of_edges[edges][:,1]]
        xB = nodes_centroids[:, 0][nodes_of_edges[edges][:,0]]
        yB = nodes_centroids[:, 1][nodes_of_edges[edges][:,0]]

        x = np.vstack([xk, xl, xA, xB]).T
        y = np.vstack([yk, yl, yA, yB]).T

        Rx = x.sum(axis=1)
        Ry = y.sum(axis=1)
        Qxx = np.power(x, 2).sum(axis=1)
        Qyy = np.power(y, 2).sum(axis=1)
        Qxy = np.diag(np.dot(x, y.T))

        D = 4*(Qxx*Qyy - np.power(Qxy, 2)) + Rx*(Qxy*Ry - Qyy*Rx) + Ry*(Qxy*Rx - Qxx*Ry)
        mi_xx = 4*Qxx - np.power(Rx, 2)
        mi_xy = Rx*Ry - 4*Qxy
        mi_yy = 4*Qyy - np.power(Ry, 2)
        mi_x = Qxy*Rx - Qxx*Ry
        mi_y = Qxy*Ry - Qyy*Rx 

        result = {
            'D': D,
            'mi_xx': mi_xx,
            'mi_xy': mi_xy,
            'mi_yy': mi_yy,
            'mi_x': mi_x,
            'mi_y': mi_y
        }

        dtype = [(x, np.float64) for x in list(result.keys())]
        
        resp = np.zeros(len(result['D']), dtype=dtype)
        for i in list(result.keys()):
            resp[i] = result[i]
        
        dtype_x_alpha = [('xk', np.float64), ('xl', np.float64), ('xA', np.float64), ('xB', np.float64)]
        dtype_y_alpha = [('yk', np.float64), ('yl', np.float64), ('yA', np.float64), ('yB', np.float64)]

        x_alpha = np.zeros(x.shape[0], dtype=dtype_x_alpha)
        y_alpha = np.zeros(y.shape[0], dtype=dtype_y_alpha)

        for i in range(4):
            x_alpha[dtype_x_alpha[i][0]][:] = x[:, i]
            y_alpha[dtype_y_alpha[i][0]][:] = y[:, i]
        
        return resp, x_alpha, y_alpha
    
    def get_xi_alpha_internal_edges(
            self,
            adjacencies,
            nodes_centroids,
            faces_centroids,
            nodes_of_edges,
            edges,
            edges_dim,
            unitary_normal_edges,
            permeability,
            bool_boundary_edges,
            **kwargs
    ):
        
        bool_internal_edges = ~bool_boundary_edges
        
        D_and_mi, x_alpha, y_alpha = self.get_D_and_mi(
            adjacencies,
            nodes_centroids,
            faces_centroids,
            nodes_of_edges,
            edges,
        )

        xy_k_sigma = self.get_x_and_y_k_sigma(
            edges,
            edges_dim,
            unitary_normal_edges,
            permeability,
            adjacencies
        )

        epsilon_alpha = np.zeros((len(D_and_mi['D']), 4))

        dtype_x_alpha = x_alpha.dtype.names
        dtype_y_alpha = y_alpha.dtype.names
        xy_alpha_local = np.zeros(3)
        xy_alpha_local[2] = 1
        m = np.zeros((2, 3))

        for edge in edges[bool_internal_edges]:
            
            data = D_and_mi[edge]
            m[:] = np.array([
                [data['mi_yy'], data['mi_xy'], data['mi_y']],
                [data['mi_xy'], data['mi_xx'], data['mi_x']]
            ])
            D = data['D']
            xksigma = xy_k_sigma[edge]
            
            for alpha in range(4):
                alpha_x = dtype_x_alpha[alpha]
                alpha_y = dtype_y_alpha[alpha]
                xy_alpha_local[0:2] = [x_alpha[alpha_x][edge], y_alpha[alpha_y][edge]]
                value = (1/D)*xksigma.dot(m.dot(xy_alpha_local))
                epsilon_alpha[edge,alpha] = value
        
        dtype_epsilon_alpha = [('k', np.float64), ('L', np.float64), ('A', np.float64), ('B', np.float64)]
        resp = np.zeros(epsilon_alpha.shape[0], dtype=dtype_epsilon_alpha)
        for i in range(4):
            alpha = dtype_epsilon_alpha[i][0]
            resp[alpha] = epsilon_alpha[:,i]
        
        return resp

    def get_M_matrix(
            self,
            faces_centroids,
            bool_boundary_edges,
            nodes_centroids,
            nodes_of_edges,
            adjacencies,
            faces,
            edges,
            unitary_normal_edges,
            permeability,
            **kwargs
    ):
        
        theta = np.pi/2
        cossin = np.array([np.cos(theta), np.sin(theta)])
        R = np.array([
            [cossin[0], -cossin[1]],
            [cossin[1],  cossin[0]]
        ])

        n_edges = len(edges)
        bool_internal_edges = ~bool_boundary_edges
        M_matrix = np.zeros((n_edges, 4, 3))
        Skl = self.get_Skl(
            edges,
            bool_boundary_edges,
            nodes_of_edges,
            adjacencies,
            unitary_normal_edges,
            nodes_centroids,
            permeability
        )

        for edge in edges[bool_internal_edges]:
            faces_adj = adjacencies[edge]
            nodes_edge = nodes_of_edges[edge]
            
            A = nodes_edge[1]
            B = nodes_edge[0]
            K = faces_adj[0]
            L = faces_adj[1]

            # AB = nodes_centroids[B] - nodes_centroids[A]
            # unitary_edge_rotated = R.dot(unitary_normal_edges[edge])
            # proj = AB.dot(unitary_edge_rotated)
            # if proj > 0:
            #     pass
            # else:
            #     A = nodes_edge[1]
            #     B = nodes_edge[0]

            local_M_matrix = np.array([
                np.concatenate([faces_centroids[K], [1.]]),
                np.concatenate([faces_centroids[L], [1.]]),
                np.concatenate([nodes_centroids[A], [1.]]),
                np.concatenate([nodes_centroids[B], [1.]])
            ])

            local_M_matrix[1,:] = local_M_matrix[1,:].dot(Skl[edge])
            M_matrix[edge] = local_M_matrix

        return M_matrix

    def get_internal_edges_flux_params(
            self,
            faces_centroids,
            bool_boundary_edges,
            nodes_centroids,
            nodes_of_edges,
            adjacencies,
            faces,
            edges,
            unitary_normal_edges,
            permeability,
            edges_dim,
            **kwargs
    ):
        
        xy_ksigma = self.get_x_and_y_k_sigma(
            edges,
            edges_dim,
            unitary_normal_edges,
            permeability,
            adjacencies
        )

        M_matrix = self.get_M_matrix(
            faces_centroids,
            bool_boundary_edges,
            nodes_centroids,
            nodes_of_edges,
            adjacencies,
            faces,
            edges,
            unitary_normal_edges,
            permeability
        )

        I_23 = np.eye(3)[0:2]

        k_params = np.zeros((len(edges), 4))
        bool_internal_edges = ~bool_boundary_edges

        for edge in edges[bool_internal_edges]:
            xy_kigma_edge = xy_ksigma[edge]
            M = M_matrix[edge]
            MM = np.linalg.inv(M.T.dot(M)).dot(M.T)
            params = xy_kigma_edge.dot(I_23).dot(MM)
            k_params[edge] = params
        
        return k_params

    def get_boundary_edges_flux_params(
            self,
            edges,
            edges_dim,
            unitary_normal_edges,
            permeability,
            adjacencies,
            bool_boundary_edges,
            nodes_of_edges,
            nodes_centroids,
            faces_centroids,
            **kwargs
    ):
        
        xy_k_sigma = self.get_x_and_y_k_sigma(
            edges,
            edges_dim,
            unitary_normal_edges,
            permeability,
            adjacencies
        )

        boundary_weights = np.zeros((bool_boundary_edges.sum(), 4))

        for edge in edges[bool_boundary_edges]:
            xy_k_sigma_edge = xy_k_sigma[edge]
            edge_nodes = nodes_of_edges[edge]
            faces_adj = adjacencies[edge]
            A_centroid = nodes_centroids[edge_nodes[1]]
            B_centroid = nodes_centroids[edge_nodes[0]]
            K_centroid = faces_centroids[faces_adj[0]]
            
            ak = K_centroid - A_centroid
            bk = K_centroid - B_centroid
            ab = B_centroid - A_centroid

            S_k_sigma_2 = np.linalg.det(np.array([ak, bk]).T)

            mi_k = (1/S_k_sigma_2)*np.dot(xy_k_sigma_edge, np.array([-ab[1], ab[0]]))
            mi_A = (1/S_k_sigma_2)*np.dot(xy_k_sigma_edge, np.array([-bk[1], bk[0]]))
            mi_B = (1/S_k_sigma_2)*np.dot(xy_k_sigma_edge, np.array([ak[1], -ak[0]]))

            internal_weights[edge,:] = [mi_k, 0, mi_A, mi_B]
        
        return boundary_weights

    def get_all_edges_flux_params(
            self,
            faces_centroids,
            bool_boundary_edges,
            nodes_centroids,
            nodes_of_edges,
            adjacencies,
            faces,
            edges,
            unitary_normal_edges,
            permeability,
            edges_dim,
            **kwargs
    ):
        """Returns the xi_params for all aldges 
            includes the xi_alpha (internal edges) and mi_alpha (boundary edges)
        
        Args:
            faces_centroids (_type_): _description_
            bool_boundary_edges (_type_): _description_
            nodes_centroids (_type_): _description_
            nodes_of_edges (_type_): _description_
            adjacencies (_type_): _description_
            faces (_type_): _description_
            edges (_type_): _description_
            unitary_normal_edges (_type_): _description_
            permeability (_type_): _description_
            edges_dim (_type_): _description_

        Returns:
            _type_: _description_
        """
        

        internal_edges_params = self.get_internal_edges_flux_params(
            faces_centroids,
            bool_boundary_edges,
            nodes_centroids,
            nodes_of_edges,
            adjacencies,
            faces,
            edges,
            unitary_normal_edges,
            permeability,
            edges_dim
        )

        boundary_edges_params = self.get_boundary_edges_flux_params(
            edges,
            edges_dim,
            unitary_normal_edges,
            permeability,
            adjacencies,
            bool_boundary_edges,
            nodes_of_edges,
            nodes_centroids,
            faces_centroids
        )

        internal_edges_params[bool_boundary_edges] = boundary_edges_params[bool_boundary_edges]

        return internal_edges_params

    def mount_problem(
            self,
            nodes_weights,
            xi_params,
            faces,
            edges,
            bool_boundary_edges,
            adjacencies,
            boundary_conditions,
            **kwargs
    ):
        """ Returns the dict with transmissibility matrix and source term

        Args:
            nodes_weights (_type_): _description_
            xi_params (_type_): _description_
            faces (_type_): _description_
            edges (_type_): _description_
            bool_boundary_edges (_type_): _description_
            adjacencies (_type_): _description_
            boundary_conditions (_type_): _description_

        Returns:
            _type_: _description_
        """



        resp = dict()
        
        
        return resp








        

        
