import numpy as np

class LsdsFluxCalculation:
    """
        Calulate least squares mpfad flux

        paper: A least squares based diamond scheme for anisotropic
                diffusion problems on polygonal meshes
                
                doi: 10.1002/fld.5031
    
    """

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
        S_alpha_sigma = np.zeros((2, 3, 3))
        S_alpha_sigma[:, [0,1], [2]] = 1
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
        
        nk_k = np.dot(unitary_normal_edge, faces_permeabilities[0])
        nk_l = np.dot(unitary_normal_edge, faces_permeabilities[1])

        S_alpha_sigma[:, 0, 0:2] = edge_nodes_centroids[0]
        S_alpha_sigma[:, 1, 0:2] = edge_nodes_centroids[1]
        S_alpha_sigma[0, 2, 0:2] = nk_k
        S_alpha_sigma[1, 2, 0:2] = nk_l
        Skl = np.linalg.inv(S_alpha_sigma[1]).dot(S_alpha_sigma[0])
        
        return Skl

    def get_x_L_barra(
            self,
            adjacencies,
            edges,
            bool_boundary_edges,
            faces_centroids,
            Skl,
            **kwargs

    ):
        
        bool_internal_edges = ~bool_boundary_edges
        faces_L = adjacencies[bool_internal_edges][:, 1]
        centroids_faces_L = faces_centroids[faces_L]
        Skl_internal_edges = Skl[bool_internal_edges]
        xl_barra = np.zeros()

        for xl, Skl_local in zip(centroids_faces_L, Skl_internal_edges):
            pass
    
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
        xA = nodes_centroids[:, 0][nodes_of_edges[edges][:,0]]
        yA = nodes_centroids[:, 1][nodes_of_edges[edges][:,0]]
        xB = nodes_centroids[:, 0][nodes_of_edges[edges][:,1]]
        yB = nodes_centroids[:, 1][nodes_of_edges[edges][:,1]]

        x = np.vstack([xk, xA, xB]).T
        y = np.vstack([yk, yA, yB]).T

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
        
        dtype_x_alpha = [('xk', np.float64), ('xA', np.float64), ('xB', np.float64)]
        dtype_y_alpha = [('yk', np.float64), ('yA', np.float64), ('yB', np.float64)]

        x_alpha = np.zeros(x.shape[0], dtype=dtype_x_alpha)
        y_alpha = np.zeros(y.shape[0], dtype=dtype_y_alpha)

        for i in range(3):
            x_alpha[dtype_x_alpha[i][0]][:] = x[:, i]
            y_alpha[dtype_y_alpha[i][0]][:] = y[:, i]
        
        return resp, x_alpha, y_alpha
    
    def get_xi_alpha(
            self,
            adjacencies,
            nodes_centroids,
            faces_centroids,
            nodes_of_edges,
            edges,
            edges_dim,
            unitary_normal_edges,
            permeability,
            **kwargs
    ):
        
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

        epsilon_alpha = np.zeros((len(D_and_mi['D']), 3))

        dtype_x_alpha = x_alpha.dtype.names
        dtype_y_alpha = y_alpha.dtype.names
        xy_alpha_local = np.zeros(3)
        xy_alpha_local[2] = 1
        m = np.zeros((2, 3))

        for edge in edges:
            
            data = D_and_mi[edge]
            m[:] = np.array([
                [data['mi_yy'], data['mi_xy'], data['mi_y']],
                [data['mi_xy'], data['mi_xx'], data['mi_x']]
            ])
            D = data['D']
            xksigma = xy_k_sigma[edge]
            
            for alpha in range(3):
                alpha_x = dtype_x_alpha[alpha]
                alpha_y = dtype_y_alpha[alpha]
                xy_alpha_local[0:2] = [x_alpha[alpha_x][edge], y_alpha[alpha_y][edge]]
                value = (1/D)*xksigma.dot(m.dot(xy_alpha_local))
                epsilon_alpha[edge,alpha] = value
        
        dtype_epsilon_alpha = [('k', np.float64), ('A', np.float64), ('B', np.float64)]
        resp = np.zeros(epsilon_alpha.shape[0], dtype=dtype_epsilon_alpha)
        for i in range(3):
            alpha = dtype_epsilon_alpha[i][0]
            resp[alpha] = epsilon_alpha[:,i]
        
        return resp

