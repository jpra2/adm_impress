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