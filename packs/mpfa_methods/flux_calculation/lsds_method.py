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

    def get_Q_and_R(
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

        import pdb; pdb.set_trace()

        return {
            'Rx': Rx,
            'Ry': Ry,
            'Qxx': Qxx,
            'Qyy': Qyy,
            'Qxy': Qxy
        }
