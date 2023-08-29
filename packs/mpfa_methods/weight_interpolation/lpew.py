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

    def create_Tk_points(self, nodes_centroids, nodes_of_edges, tau=0.5):
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





    def create_knt_vef(self, tk_points, nodes_centroids, faces_centroids, permeability, edges_of_nodes, nodes_of_edges, nodes, edges, adjacencies, **kwargs):

        R_matrix = self.get_Rmatrix()

        all_face_id = []
        all_vertice_id = []
        all_edge_id = []
        all_kn = []
        all_kt = []

        for edge in edges:
            nodes_edge = nodes_of_edges[edge]



