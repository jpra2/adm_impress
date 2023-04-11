from packs.manager.arraydatamanager import SuperArrayManager
import numpy as np

class CalculateGlsWeight2D:

    def __init__(self, **kwargs):
        """ Initialize the class
            
            paper: A least squares based diamond scheme for anisotropic
            diffusion problems on polygonal meshes
            
            doi: 10.1002/fld.5031

            tese de tulio: nos no contorno

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
            unitary_normal_edges: unitary normal vector of edges
        """
        # TODO ainda falta implementar os pesos para a condicao de fluxo prescrito

        self.adjacencies = kwargs.get('adjacencies')
        self.faces = kwargs.get('faces')
        self.nodes_of_nodes = kwargs.get('nodes_of_nodes')
        self.edges_of_nodes = kwargs.get('edges_of_nodes')
        self.nodes_of_edges = kwargs.get('nodes_of_edges')
        self.faces_of_nodes = kwargs.get('faces_of_nodes')
        self.nodes_centroids = kwargs.get('nodes_centroids')
        self.edges = kwargs.get('edges')
        self.bool_boundary_edges = kwargs.get('bool_boundary_edges')
        self.bool_boundary_nodes = kwargs.get('bool_boundary_nodes')
        self.bool_internal_nodes = ~self.bool_boundary_nodes
        self.nodes = kwargs.get('nodes')
        self.faces_centroids = kwargs.get('faces_centroids')
        self.permeability = kwargs.get('permeability')
        self.unitary_normal_edges = kwargs.get('unitary_normal_edges')

    def get_weights_internal_nodes(self):
        
        nodes_ids = []
        faces_ids = []
        all_weight = []

        for node in self.nodes[self.bool_internal_nodes]:
            nodes_adj = self.nodes_of_nodes[node]
            edges_adj = self.edges_of_nodes[node]
            centroid_node = self.nodes_centroids[node]
            centroids_nodes_adj = self.nodes_centroids[nodes_adj]
            faces_adj = self.faces_of_nodes[node]
            centroids_faces_adj = self.faces_centroids[faces_adj]

            n = len(faces_adj)
            # mfaces = self.mfaces(n, centroids_faces_adj, centroid_node)
            # mnodes = self.mnodes(len(nodes_adj), centroids_nodes_adj, centroid_node)
            # mnormalperm = self.mnormalperm(n, faces_adj, edges_adj, self.adjacencies, self.unitary_normal_edges, self.permeability)
            # Mv = self.Mv(n, mfaces, mnodes, mnormalperm)
            # Nv = self.Nv(n)            
            # M = self.M(Mv, Nv)

            M = self.M(
                self.Mv(
                    n,
                    self.mfaces(n, centroids_faces_adj, centroid_node),
                    self.mnodes(len(nodes_adj), centroids_nodes_adj, centroid_node),
                    self.mnormalperm(n, faces_adj, edges_adj, self.adjacencies, self.unitary_normal_edges, self.permeability)
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

    def eT(self, n):
        eT = np.zeros(2*n+1)
        eT[-1] = 1
        return eT

    def mfaces(self, n_mfaces: int, centroids_faces_adj: np.ndarray, centroid_node: np.ndarray):
        mfaces = np.zeros((n_mfaces, 2*n_mfaces))
        f_dists = centroids_faces_adj - centroid_node
        for i in range(n_mfaces):
            mfaces[i,2*i:2*i+2] = f_dists[i]
        
        return mfaces

    def mnodes(self, n_mnodes, centroids_nodes_adj, centroid_node):
        
        mnodes = np.zeros((n_mnodes, 2*n_mnodes))
        nodes_dists = centroids_nodes_adj - centroid_node
        for i in range(n_mnodes):
            mnodes[i,2*i:2*i+2] = nodes_dists[i]
            
        mnodes[0, 2*n_mnodes-2:] = -nodes_dists[0]
        
        for i in range(1, n_mnodes):
            mnodes[i, 2*i-2:2*i] = -nodes_dists[i]
        
        return mnodes

    def mnormalperm(self, n_mfaces, faces_adj, edges_adj, adjacencies, unitary_normal_edges, permeability):

        n_mnormalperm = n_mfaces
        mnormalperm = np.zeros((n_mnormalperm, 2*n_mnormalperm))
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
            
            mnormalperm[i, 2*i:2*i+2] = np.dot(normal, face_perm)
        
        for i in range(1, n_mnormalperm):
            face = faces_adj[i]
            edge = edges_adj[i]
            faces_adj_edge = adjacencies[edge]
            unitary_normal_edge = unitary_normal_edges[edge]
            face_perm = permeability[faces_adj_edge[faces_adj_edge!=face][0]]
            if face == faces_adj_edge[1]:
                normal = -unitary_normal_edge
            else:
                normal = unitary_normal_edge
            
            mnormalperm[i, 2*i-2:2*i] = -np.dot(normal, face_perm)
        
        i = 0
        face = faces_adj[i]
        edge = edges_adj[i]
        faces_adj_edge = adjacencies[edge]
        unitary_normal_edge = unitary_normal_edges[edge]
        face_perm = permeability[faces_adj_edge[faces_adj_edge!=face][0]]
        if face == faces_adj_edge[1]:
            normal = -unitary_normal_edge
        else:
            normal = unitary_normal_edge
        
        mnormalperm[0, 2*n_mnormalperm-2:] = -np.dot(normal, face_perm)

        return mnormalperm

    def Mv(self, n, mfaces, mnodes, mnormalperm):
        
        n_mfaces = n
        Mv = np.zeros((3*n_mfaces, 2*n_mfaces+1))
        Mv[0:n_mfaces,-1] = 1
        Mv[0:n_mfaces, 0:2*n_mfaces] = mfaces
        Mv[n_mfaces:2*n_mfaces, 0:2*n_mfaces] = mnodes
        Mv[2*n_mfaces:3*n_mfaces, 0:2*n_mfaces] = mnormalperm
        
        return Mv

    def Nv(self, n):
        
        n_mfaces = n
        Nv = np.zeros((3*n_mfaces, n_mfaces))
        Nv[np.arange(n_mfaces), np.arange(n_mfaces)] = 1

        return Nv

    def M(self, Mv, Nv):
        return np.linalg.inv(Mv.T.dot(Mv)).dot(Mv.T).dot(Nv)
    
    def weights(self, eT, M):
        return eT.dot(M)
    
    def get_weights_bnodes(self):
        
        nodes_ids = []
        faces_ids = []
        all_weight = []

        internal_nodes = self.nodes[self.bool_internal_nodes]

        for node in self.nodes[self.bool_boundary_nodes]:
            nodes_adj = self.nodes_of_nodes[node]
            edges_adj = self.edges_of_nodes[node]
            faces_adj = self.faces_of_nodes[node]
            local_internal_nodes = np.intersect1d(internal_nodes, nodes_adj)
            centroid_node = self.nodes_centroids[node]
            centroids_faces_adj = self.faces_centroids[faces_adj]
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
                    self.nodes_centroids,
                    map_faces_adj,
                    self.adjacencies,
                    nodes_adj
                ),
                self.bnormalperm(
                    edges_adj.shape[0],
                    faces_adj.shape[0],
                    edges_adj,
                    self.adjacencies,
                    map_faces_adj,
                    self.unitary_normal_edges,
                    self.permeability
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
            
    def bmnodes(self, n_nodes, n_faces, internal_nodes_adj, edges_adj, centroid_node, nodes_centroids, map_faces_adj, adjacencies, nodes_adj):
        
        if n_nodes == 0:
            return np.array([])
        mnodes = np.zeros((n_nodes, 2*n_faces))
        for i, node_adj in enumerate(internal_nodes_adj):
            centroid_node_adj = nodes_centroids[node_adj]
            node_dist_v = centroid_node_adj - centroid_node
            edge_corresp = edges_adj[nodes_adj == node_adj][0]
            faces_of_edge_corresp = adjacencies[edge_corresp]
            
            face_r = map_faces_adj[faces_of_edge_corresp[1]]
            face_l = map_faces_adj[faces_of_edge_corresp[0]]
            # if map_faces_adj[faces_of_edge_corresp[1]] > map_faces_adj[faces_of_edge_corresp[0]]:
            #     face_r = map_faces_adj[faces_of_edge_corresp[0]]
            #     face_l = map_faces_adj[faces_of_edge_corresp[1]]
            
            mnodes[i, 2*face_r:2*face_r+2] = -node_dist_v
            mnodes[i, 2*face_l:2*face_l+2] = node_dist_v
        
        return mnodes

    def bnormalperm(self, n_edges, n_faces, edges_adj, adjacencies, map_faces_adj, unitary_normal_edges, permeability):
        mnormalperm = np.zeros((n_edges, 2*n_faces))

        for i, edge in enumerate(edges_adj):
            faces_edge = adjacencies[edge]
            unitary_normal_edge = unitary_normal_edges[edge]
            id_local_faces_edge = map_faces_adj[faces_edge]

            value_l = np.dot(unitary_normal_edge, permeability[faces_edge[0]])

            mnormalperm[i, 2*id_local_faces_edge[0]:2*id_local_faces_edge[0]+2] = value_l

            if faces_edge[1] == -1:
                continue

            value_r = np.dot(unitary_normal_edge, permeability[faces_edge[1]])
            mnormalperm[i, 2*id_local_faces_edge[1]:2*id_local_faces_edge[1]+2] = -value_r

        
        return mnormalperm

    def bM(self, n_faces, n_nodes, n_edges, mfaces, mnodes, mnormalperm):
        M = np.zeros((n_faces+n_nodes+n_edges, 2*n_faces + 1))
        M[0:n_faces,0:2*n_faces] = mfaces
        M[0:n_faces,-1] = 1
        M[n_faces:n_faces+n_nodes, 0:2*n_faces] = mnodes
        M[n_faces+n_nodes:n_faces+n_nodes+n_edges, 0:2*n_faces] = mnormalperm

        return M

    def bN(self, n_faces, n_nodes, n_edges):
        N = np.zeros((n_faces+n_nodes+n_edges, n_faces))
        vec = np.arange(n_faces)
        N[vec, vec] = 1

        return N

    def get_nodes_weights(self):

        dtype = [('node_id', int), ('face_id', int), ('weight', float)]

        nodes_ids, faces_ids, weights = self.get_weights_internal_nodes()
        nodes_ids2, faces_ids2, weights2 = self.get_weights_bnodes()

        nodes_ids = np.concatenate([nodes_ids, nodes_ids2])
        faces_ids = np.concatenate([faces_ids, faces_ids2])
        weights = np.concatenate([weights, weights2])

        nodes_weights = np.zeros(len(nodes_ids), dtype=dtype)
        nodes_weights['node_id'] = nodes_ids
        nodes_weights['face_id'] = faces_ids
        nodes_weights['weight'] = weights

        return nodes_weights


def get_gls_nodes_weights(**kwargs):

    """
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
            unitary_normal_edges: unitary normal vector of edges
    """

    calculate_weights = CalculateGlsWeight2D(**kwargs)
    nodes_weights = calculate_weights.get_nodes_weights()

    return nodes_weights






