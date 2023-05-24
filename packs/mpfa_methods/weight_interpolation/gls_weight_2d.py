import numpy as np
import scipy.sparse as sp


class CalculateGlsWeight2D:
    """ Calulate vertices weights from gls method

    paper: A least squares based diamond scheme for anisotropic
            diffusion problems on polygonal meshes
            
            doi: 10.1002/fld.5031

            tese de tulio: vertice no contorno
    """

    def get_weights_internal_nodes(
            self,
            nodes,
            bool_boundary_nodes,
            nodes_of_nodes,
            edges_of_nodes,
            nodes_centroids,
            faces_of_nodes,
            faces_centroids,
            adjacencies,
            unitary_normal_edges,
            permeability,
            nodes_to_calculate,
            **kwargs
    ):

        nodes_ids = []
        faces_ids = []
        all_weight = []
        bool_internal_nodes = ~bool_boundary_nodes
        nodes_to_iterate = np.intersect1d(nodes[bool_internal_nodes], nodes_to_calculate)

        for node in nodes_to_iterate:
            nodes_adj = nodes_of_nodes[node]
            edges_adj = edges_of_nodes[node]
            centroid_node = nodes_centroids[node]
            centroids_nodes_adj = nodes_centroids[nodes_adj]
            faces_adj = faces_of_nodes[node]
            centroids_faces_adj = faces_centroids[faces_adj]

            n = len(faces_adj)

            M = self.M(
                self.Mv(
                    n,
                    self.mfaces(n, centroids_faces_adj, centroid_node),
                    self.mnodes(len(nodes_adj), centroids_nodes_adj, centroid_node),
                    self.mnormalperm(n, faces_adj, edges_adj, adjacencies, unitary_normal_edges, permeability)
                ),
                self.Nv(n)
            )

            weights = self.eT(n).dot(M)

            nodes_ids.append(np.repeat(node, n))
            faces_ids.append(faces_adj)
            all_weight.append(weights)

        nodes_ids = np.concatenate(nodes_ids)
        faces_ids = np.concatenate(faces_ids)
        all_weight = np.concatenate(all_weight)

        return nodes_ids, faces_ids, all_weight

    @staticmethod
    def eT(n):
        eT = np.zeros(2 * n + 1)
        eT[-1] = 1
        return eT

    @staticmethod
    def mfaces(n_mfaces: int, centroids_faces_adj: np.ndarray, centroid_node: np.ndarray):
        mfaces = np.zeros((n_mfaces, 2 * n_mfaces))
        f_dists = centroids_faces_adj - centroid_node
        for i in range(n_mfaces):
            mfaces[i, 2 * i:2 * i + 2] = f_dists[i]

        return mfaces

    @staticmethod
    def mnodes(n_mnodes, centroids_nodes_adj, centroid_node):

        mnodes = np.zeros((n_mnodes, 2 * n_mnodes))
        nodes_dists = centroids_nodes_adj - centroid_node
        for i in range(n_mnodes):
            mnodes[i, 2 * i:2 * i + 2] = nodes_dists[i]

        mnodes[0, 2 * n_mnodes - 2:] = -nodes_dists[0]

        for i in range(1, n_mnodes):
            mnodes[i, 2 * i - 2:2 * i] = -nodes_dists[i]

        return mnodes

    @staticmethod
    def mnormalperm(n_mfaces, faces_adj, edges_adj, adjacencies, unitary_normal_edges, permeability):

        n_mnormalperm = n_mfaces
        mnormalperm = np.zeros((n_mnormalperm, 2 * n_mnormalperm))
        for i in range(n_mnormalperm):
            face = faces_adj[i]
            edge = edges_adj[i]
            faces_adj_edge = adjacencies[edge]
            unitary_normal_edge = unitary_normal_edges[edge]
            face_perm = permeability[face]
            if face == faces_adj_edge[1]:
                normal = -unitary_normal_edge
            else:
                normal = unitary_normal_edge

            mnormalperm[i, 2 * i:2 * i + 2] = np.dot(normal, face_perm)

        for i in range(1, n_mnormalperm):
            face = faces_adj[i]
            edge = edges_adj[i]
            faces_adj_edge = adjacencies[edge]
            unitary_normal_edge = unitary_normal_edges[edge]
            face_perm = permeability[faces_adj_edge[faces_adj_edge != face][0]]
            if face == faces_adj_edge[1]:
                normal = -unitary_normal_edge
            else:
                normal = unitary_normal_edge

            mnormalperm[i, 2 * i - 2:2 * i] = -np.dot(normal, face_perm)

        i = 0
        face = faces_adj[i]
        edge = edges_adj[i]
        faces_adj_edge = adjacencies[edge]
        unitary_normal_edge = unitary_normal_edges[edge]
        face_perm = permeability[faces_adj_edge[faces_adj_edge != face][0]]
        if face == faces_adj_edge[1]:
            normal = -unitary_normal_edge
        else:
            normal = unitary_normal_edge

        mnormalperm[0, 2 * n_mnormalperm - 2:] = -np.dot(normal, face_perm)

        return mnormalperm

    @staticmethod
    def Mv(n, mfaces, mnodes, mnormalperm):

        n_mfaces = n
        Mv = np.zeros((3 * n_mfaces, 2 * n_mfaces + 1))
        Mv[0:n_mfaces, -1] = 1
        Mv[0:n_mfaces, 0:2 * n_mfaces] = mfaces
        Mv[n_mfaces:2 * n_mfaces, 0:2 * n_mfaces] = mnodes
        Mv[2 * n_mfaces:3 * n_mfaces, 0:2 * n_mfaces] = mnormalperm

        return Mv

    @staticmethod
    def Nv(n):

        n_mfaces = n
        Nv = np.zeros((3 * n_mfaces, n_mfaces))
        Nv[np.arange(n_mfaces), np.arange(n_mfaces)] = 1

        return Nv

    @staticmethod
    def M(Mv, Nv):
        return np.linalg.inv(Mv.T.dot(Mv)).dot(Mv.T).dot(Nv)

    @staticmethod
    def weights(eT, M):
        return eT.dot(M)

    def get_weights_bnodes(
            self,
            nodes,
            bool_boundary_nodes,
            nodes_of_nodes,
            edges_of_nodes,
            faces_of_nodes,
            nodes_centroids,
            faces_centroids,
            adjacencies,
            unitary_normal_edges,
            permeability,
            nodes_to_calculate,
            **kwargs
    ):

        nodes_ids = []
        faces_ids = []
        all_weight = []

        bool_internal_nodes = ~bool_boundary_nodes
        internal_nodes = nodes[bool_internal_nodes]
        nodes_to_iterate = np.intersect1d(nodes[bool_boundary_nodes], nodes_to_calculate)

        for node in nodes_to_iterate:
            nodes_adj = nodes_of_nodes[node]
            edges_adj = edges_of_nodes[node]
            faces_adj = faces_of_nodes[node]
            local_internal_nodes = np.intersect1d(internal_nodes, nodes_adj)
            centroid_node = nodes_centroids[node]
            centroids_faces_adj = faces_centroids[faces_adj]
            map_faces_adj = np.repeat(-1, faces_adj.max() + 3)
            map_faces_adj[faces_adj] = np.arange(faces_adj.shape[0])

            M = self.bM(
                faces_adj.shape[0],
                local_internal_nodes.shape[0],
                edges_adj.shape[0],
                self.mfaces(
                    faces_adj.shape[0],
                    centroids_faces_adj,
                    centroid_node
                ),
                self.bmnodes(
                    local_internal_nodes.shape[0],
                    faces_adj.shape[0],
                    local_internal_nodes,
                    edges_adj,
                    centroid_node,
                    nodes_centroids,
                    map_faces_adj,
                    adjacencies,
                    nodes_adj
                ),
                self.bnormalperm(
                    edges_adj.shape[0],
                    faces_adj.shape[0],
                    edges_adj,
                    adjacencies,
                    map_faces_adj,
                    unitary_normal_edges,
                    permeability
                )
            )

            N = self.bN(
                faces_adj.shape[0],
                local_internal_nodes.shape[0],
                edges_adj.shape[0]
            )

            weights = self.eT(faces_adj.shape[0]).dot(np.linalg.inv(M.T.dot(M)).dot(M.T).dot(N))

            nodes_ids.append(np.repeat(node, weights.shape[0]))
            faces_ids.append(faces_adj)
            all_weight.append(weights)

        nodes_ids = np.concatenate(nodes_ids)
        faces_ids = np.concatenate(faces_ids)
        all_weight = np.concatenate(all_weight)

        return nodes_ids, faces_ids, all_weight

    @staticmethod
    def bmnodes(n_nodes, n_faces, internal_nodes_adj, edges_adj, centroid_node, nodes_centroids, map_faces_adj,
                adjacencies, nodes_adj):

        if n_nodes == 0:
            return np.array([])
        mnodes = np.zeros((n_nodes, 2 * n_faces))
        for i, node_adj in enumerate(internal_nodes_adj):
            centroid_node_adj = nodes_centroids[node_adj]
            node_dist_v = centroid_node_adj - centroid_node
            edge_corresp = edges_adj[nodes_adj == node_adj][0]
            faces_of_edge_corresp = adjacencies[edge_corresp]

            face_r = map_faces_adj[faces_of_edge_corresp[1]]
            face_l = map_faces_adj[faces_of_edge_corresp[0]]

            mnodes[i, 2 * face_r:2 * face_r + 2] = -node_dist_v
            mnodes[i, 2 * face_l:2 * face_l + 2] = node_dist_v

        return mnodes

    @staticmethod
    def bnormalperm(n_edges, n_faces, edges_adj, adjacencies, map_faces_adj, unitary_normal_edges, permeability):
        mnormalperm = np.zeros((n_edges, 2 * n_faces))

        for i, edge in enumerate(edges_adj):
            faces_edge = adjacencies[edge]
            unitary_normal_edge = unitary_normal_edges[edge]
            id_local_faces_edge = map_faces_adj[faces_edge]

            value_l = np.dot(unitary_normal_edge, permeability[faces_edge[0]])

            mnormalperm[i, 2 * id_local_faces_edge[0]:2 * id_local_faces_edge[0] + 2] = value_l

            if faces_edge[1] == -1:
                continue

            value_r = np.dot(unitary_normal_edge, permeability[faces_edge[1]])
            mnormalperm[i, 2 * id_local_faces_edge[1]:2 * id_local_faces_edge[1] + 2] = -value_r

        return mnormalperm

    @staticmethod
    def bM(n_faces, n_nodes, n_edges, mfaces, mnodes, mnormalperm):
        
        M = np.zeros((n_faces + n_nodes + n_edges, 2 * n_faces + 1))
        M[0:n_faces, 0:2 * n_faces] = mfaces
        M[0:n_faces, -1] = 1
        if n_nodes > 0:
            M[n_faces:n_faces + n_nodes, 0:2 * n_faces] = mnodes
        M[n_faces + n_nodes:n_faces + n_nodes + n_edges, 0:2 * n_faces] = mnormalperm
            
        return M

    @staticmethod
    def bN(n_faces, n_nodes, n_edges):
        N = np.zeros((n_faces + n_nodes + n_edges, n_faces))
        vec = np.arange(n_faces)
        N[vec, vec] = 1

        return N

    def get_nodes_weights(self, **kwargs):

        """
            paper: A least squares based diamond scheme for anisotropic
                diffusion problems on polygonal meshes
                
                doi: 10.1002/fld.5031

                tese de tulio: vertice no contorno

            @param kwargs: dict with the keys:
                adjacencies: faces adjacencies of edges,
                faces: global ids of faces,
                nodes_of_nodes: nodes adjacencies of nodes
                edges_of_nodes: edges adjacencies of nodes
                nodes_of_edges: nodes adjacencies of edges
                faces_of_nodes: faces adjacencies of nodes
                nodes_centroids: centroids of mesh nodes,
                edges: global ids of mesh edges,
                bool_boundary_edges: bool vector of len(edges) with true in boundary edges,
                bool_boundary_nodes: bool vector of len(nodes) with true in boundary nodes
                nodes: global id of mesh nodes
                faces_centroids: centroids of faces
                permeability: mesh permeability
                unitary_normal_edges: unitary normal vector of edges. this vector point 
                                      to outward of left face from adjacencies vector
                nodes_to_calculate: nodes for weight calculation
        """

        dtype = [('node_id', int), ('face_id', int), ('weight', float)]

        nodes_ids, faces_ids, weights = self.get_weights_internal_nodes(**kwargs)
        nodes_ids2, faces_ids2, weights2 = self.get_weights_bnodes(**kwargs)

        nodes_ids = np.concatenate([nodes_ids, nodes_ids2])
        faces_ids = np.concatenate([faces_ids, faces_ids2])
        weights = np.concatenate([weights, weights2])

        nodes_weights = np.zeros(len(nodes_ids), dtype=dtype)
        nodes_weights['node_id'] = nodes_ids
        nodes_weights['face_id'] = faces_ids
        nodes_weights['weight'] = weights

        return nodes_weights

def mount_weight_matrix(nodes_weights):
    n_faces = len(np.unique(nodes_weights['face_id']))
    n_nodes = len(np.unique(nodes_weights['node_id']))
    mweight = np.zeros((n_nodes, n_faces))

    lines = np.array([], dtype=int)
    cols = lines.copy()
    data = np.array([], dtype=np.float64)

    for node in np.unique(nodes_weights['node_id']):
        faces = nodes_weights['face_id'][nodes_weights['node_id'] == node]
        weights = nodes_weights['weight'][nodes_weights['node_id'] == node]
        lines = np.append(lines, np.repeat(node, faces.shape[0]))
        cols = np.append(cols, faces)
        data = np.append(data, weights)
    
    mweight[lines, cols] = data
    
    return mweight

def mount_sparse_weight_matrix(nodes_weights):
    n_faces = len(np.unique(nodes_weights['face_id']))
    n_nodes = len(np.unique(nodes_weights['node_id']))

    lines = np.array([], dtype=int)
    cols = lines.copy()
    data = np.array([], dtype=np.float64)

    for node in np.unique(nodes_weights['node_id']):
        faces = nodes_weights['face_id'][nodes_weights['node_id'] == node]
        weights = nodes_weights['weight'][nodes_weights['node_id'] == node]
        lines = np.append(lines, np.repeat(node, faces.shape[0]))
        cols = np.append(cols, faces)
        data = np.append(data, weights)
    
    mweight = sp.csr_matrix((data, (lines, cols)), shape=(n_nodes, n_faces))
    
    return mweight


def get_gls_nodes_weights(**kwargs):
    """
    paper: A least squares based diamond scheme for anisotropic
            diffusion problems on polygonal meshes
            
            doi: 10.1002/fld.5031

            tese de tulio: vertice no contorno

        @param kwargs: dict with the keys:
            adjacencies: faces adjacencies of edges,
            faces: global ids of faces,
            nodes_of_nodes: nodes adjacencies of nodes by edges
            edges_of_nodes: edges adjacencies of nodes
            nodes_of_edges: nodes adjacencies of edges
            faces_of_nodes: faces adjacencies of nodes
            nodes_centroids: centroids of mesh nodes,
            edges: global ids of mesh edges,
            bool_boundary_edges: bool vector of len(edges) with true in boundary edges,
            bool_boundary_nodes: bool vector of len(nodes) with true in boundary nodes
            nodes: global id of mesh nodes
            faces_centroids: centroids of faces
            permeability: mesh permeability
            unitary_normal_edges: unitary normal vector of edges.
                                  this vector point to outward of left face from adjacencies vector
            nodes_to_calculate: nodes for weight calculation
    """

    calculate_weights = CalculateGlsWeight2D()
    nodes_weights = calculate_weights.get_nodes_weights(**kwargs)

    return nodes_weights
