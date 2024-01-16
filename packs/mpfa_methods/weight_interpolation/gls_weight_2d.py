import numpy as np
import scipy.sparse as sp
from packs.manager.boundary_conditions import BoundaryConditions
from packs import defnames



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
            nodes_of_edges,
            **kwargs
    ):

        nodes_ids = []
        faces_ids = []
        all_weight = []
        bool_internal_nodes = ~bool_boundary_nodes
        nodes_to_iterate = np.intersect1d(nodes[bool_internal_nodes], nodes_to_calculate)

        for node in nodes_to_iterate:
            edges_adj = edges_of_nodes[node]
            nodes_edges_adj = nodes_of_edges[edges_adj]
            nodes_adj = nodes_edges_adj[nodes_edges_adj!=node]
            centroid_node = nodes_centroids[node]
            faces_adj = faces_of_nodes[node]
            centroids_faces_adj = faces_centroids[faces_adj]

            n = len(faces_adj)
            n_nodes = len(nodes_adj)
            n_edges = len(edges_adj)
            n_faces = n

            local_index_faces = np.arange(n_faces)
            local_index_edges = np.arange(n_edges)

            mfaces = np.zeros((n_faces, 2 * n_faces))
            f_dists = centroids_faces_adj - centroid_node
            for i in range(n_faces):
                mfaces[i, 2 * i:2 * i + 2] = f_dists[i]
            
            mnodes = np.zeros((n_nodes, 2*n_faces))
            mnormal_perm = np.zeros((n_edges, 2*n_faces))


            for i in range(n_edges):
                edge = edges_adj[i]
                node_adj = nodes_adj[i]
                unitary_normal_edge = unitary_normal_edges[edge]
                id_edge = local_index_edges[edges_adj==edge]
                id_node = id_edge
                faces_adj_edge = adjacencies[edge]
                tau = nodes_centroids[node_adj] - centroid_node
                id_face0 = local_index_faces[faces_adj == faces_adj_edge[0]][0]
                id_face1 = local_index_faces[faces_adj == faces_adj_edge[1]][0]

                mnodes[id_node, 2*id_face0: 2*id_face0+2] = +tau
                mnodes[id_node, 2*id_face1: 2*id_face1+2] = -tau

                mnormal_perm[id_edge, 2*id_face0: 2*id_face0+2] = np.dot(unitary_normal_edge, permeability[faces_adj_edge[0]])
                mnormal_perm[id_edge, 2*id_face1: 2*id_face1+2] = np.dot(-unitary_normal_edge, permeability[faces_adj_edge[1]])

            M = self.M(
                self.Mv(
                    n,
                    mfaces,
                    mnodes,
                    mnormal_perm
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

        # for node in np.unique(nodes_ids):
        #     test = nodes_ids == node
        #     print(all_weight[test].sum())
        #     all_weight[test] = all_weight[test]/all_weight[test].sum()
        #     print(all_weight[test].sum())
        #     import pdb; pdb.set_trace()

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
            edges,
            bool_boundary_nodes,
            nodes_of_edges,
            edges_of_nodes,
            faces_of_nodes,
            nodes_centroids,
            faces_centroids,
            adjacencies,
            unitary_normal_edges,
            permeability,
            nodes_to_calculate,
            bool_boundary_edges,
            neumann_edges = np.array([]),
            neumann_edges_value = np.array([]),
            **kwargs
    ):

        nodes_ids = []
        faces_ids = []
        all_weight = []

        nodes_to_iterate = np.intersect1d(nodes[bool_boundary_nodes], nodes_to_calculate)
        bool_biedges = ~bool_boundary_edges
        all_neumann_weights = []
        
        for node in nodes_to_iterate:
            edges_adj = edges_of_nodes[node]
            faces_adj = faces_of_nodes[node]
            internal_edges = np.intersect1d(edges_adj, edges[bool_biedges])
            nodes_edges_internal = nodes_of_edges[internal_edges]
            nodes_adj = nodes_edges_internal[nodes_edges_internal!=node]

            centroid_node = nodes_centroids[node]
            centroids_faces_adj = faces_centroids[faces_adj]

            n = len(faces_adj)
            n_nodes = len(nodes_adj)
            n_edges = len(edges_adj)
            n_faces = n

            local_index_nodes = np.arange(n_nodes)
            local_index_faces = np.arange(n_faces)
            local_index_edges = np.arange(n_edges)

            mfaces = np.zeros((n_faces, 2 * n_faces))
            f_dists = centroids_faces_adj - centroid_node
            for i in range(n_faces):
                mfaces[i, 2 * i:2 * i + 2] = f_dists[i]

            if n_nodes == 0:
                mnodes = np.array([])
            else:
                mnodes = np.zeros((n_nodes, 2*n_faces))
                for i in local_index_nodes:
                    tau = nodes_centroids[nodes_adj[i]] - centroid_node
                    face0 = adjacencies[internal_edges[i], 0]
                    id_face0 = local_index_faces[faces_adj == face0][0]
                    mnodes[i, 2*id_face0: 2*id_face0+2] = +tau
                    face1 = adjacencies[internal_edges[i], 1]                    
                    id_face1 = local_index_faces[faces_adj == face1][0]
                    mnodes[i, 2*id_face1: 2*id_face1+2] = -tau

            mnormal_perm = np.zeros((n_edges, 2*n_faces))

            for i in local_index_edges:
                edge = edges_adj[i]
                unitary_normal_edge = unitary_normal_edges[edge]
                id_edge = i
                faces_adj_edge = adjacencies[edge]
                face0 = faces_adj_edge[0]
                id_face0 = local_index_faces[faces_adj == face0][0]
                mnormal_perm[id_edge, 2*id_face0: 2*id_face0+2] = np.dot(unitary_normal_edge, permeability[faces_adj_edge[0]])
                face1 = faces_adj_edge[1]
                if face1 == -1:
                    continue
                id_face1 = local_index_faces[faces_adj == face1][0]
                mnormal_perm[id_edge, 2*id_face1: 2*id_face1+2] = np.dot(-unitary_normal_edge, permeability[faces_adj_edge[1]])

            M = self.bM(
                n_faces,
                n_nodes,
                n_edges,
                mfaces,
                mnodes,
                mnormal_perm
            )

            N = self.bN(
                n_faces,
                n_nodes,
                n_edges
            )    

            weights = self.eT(n_faces).dot(np.linalg.inv(M.T.dot(M)).dot(M.T).dot(N))
            
            neummann_edges_node = np.intersect1d(neumann_edges, edges_adj)
            if neummann_edges_node.shape[0] > 0:
                F = np.zeros((n_faces + n_nodes + n_edges, 1))
                local_map_edges_node = np.arange(n_faces+n_nodes, n_faces+n_nodes+n_edges)
                
                for edge in neummann_edges_node:
                    loc = edges_adj == edge
                    F[local_map_edges_node[loc]] = neumann_edges_value[neumann_edges == edge]
                
                neumann_weight = self.eT(faces_adj.shape[0]).dot(np.linalg.inv(M.T.dot(M)).dot(M.T).dot(F))
                all_neumann_weights.append([node, neumann_weight[0]])
                
            nodes_ids.append(np.repeat(node, weights.shape[0]))
            faces_ids.append(faces_adj)
            all_weight.append(weights)
        
        nodes_ids = np.concatenate(nodes_ids)
        faces_ids = np.concatenate(faces_ids)
        all_weight = np.concatenate(all_weight)
        all_neumann_weights = np.array(all_neumann_weights)

        return nodes_ids, faces_ids, all_weight, all_neumann_weights

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
                permeability: mesh faces permeability
                unitary_normal_edges: unitary normal vector of edges. this vector point 
                                      to outward of left face from adjacencies vector
                nodes_to_calculate: nodes for weight calculation
                neumann_edges: edges to set neumann boundary conditions (set [] if not exists)
                neumann_edges_value: value of neumann boundary condition (set [] if not exists) 
        """

        dtype = [('node_id', int), ('face_id', int), ('weight', float)]
        neumann_dtype = [('node_id', int), ('nweight', float)]

        nodes_ids, faces_ids, weights = self.get_weights_internal_nodes(**kwargs)
        nodes_ids2, faces_ids2, weights2, all_neumann_weights = self.get_weights_bnodes(**kwargs)

        nodes_ids = np.concatenate([nodes_ids, nodes_ids2])
        faces_ids = np.concatenate([faces_ids, faces_ids2])
        weights = np.concatenate([weights, weights2])

        nodes_weights = np.zeros(len(nodes_ids), dtype=dtype)
        nodes_weights['node_id'] = nodes_ids
        nodes_weights['face_id'] = faces_ids
        nodes_weights['weight'] = weights
        
        neumann_weights = np.zeros(len(all_neumann_weights), dtype=neumann_dtype)
        if len(all_neumann_weights) > 0:
            
            neumann_weights['node_id'] = all_neumann_weights[:, 0]
            neumann_weights['nweight'] = all_neumann_weights[:, 1]
        
        resp = {
            'nodes_weights': nodes_weights,
            'neumann_weights': neumann_weights
        }
        
        return resp
        
        
        
        


    
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

    lines = nodes_weights['node_id']
    cols = nodes_weights['face_id']
    data = nodes_weights['weight']

    # for node in np.unique(nodes_weights['node_id']):
    #     faces = nodes_weights['face_id'][nodes_weights['node_id'] == node]
    #     weights = nodes_weights['weight'][nodes_weights['node_id'] == node]
    #     lines = np.append(lines, np.repeat(node, faces.shape[0]))
    #     cols = np.append(cols, faces)
    #     data = np.append(data, weights)
    
    mweight = sp.csr_matrix((data, (lines, cols)), shape=(n_nodes, n_faces))
    
    return mweight

def get_weight_matrix_structure(nodes_weights):
    weights_matrix = mount_sparse_weight_matrix(nodes_weights)
    structure = sp.find(weights_matrix)
    n = len(structure[0])
    dtype = [('row', np.int), ('col', np.int), ('data', np.float64)]
    array = np.zeros(n, dtype=dtype)
    array['row'] = structure[0]
    array['col'] = structure[1]
    array['data'] = structure[2]

    return {defnames.nodes_weights_matrix_structure: array}
    
def mount_sparse_matrix_from_structure(nodes_weights_matrix_structure, n_nodes, n_faces):
    row = nodes_weights_matrix_structure['row']
    col = nodes_weights_matrix_structure['col']
    data = nodes_weights_matrix_structure['data']

    M = sp.csr_matrix((data, (row, col)), shape=(n_nodes, n_faces))
    return M

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
            permeability: mesh faces permeability
            unitary_normal_edges: unitary normal vector of edges.
                                  this vector point to outward of left face from adjacencies vector
            nodes_to_calculate: nodes for weight calculation
            neumann_edges: edges to set neumann boundary conditions (set [] if not exists)
            neumann_edges_value: value of neumann boundary condition (set [] if not exists) 
    """

    calculate_weights = CalculateGlsWeight2D()
    nodes_weights = calculate_weights.get_nodes_weights(**kwargs)

    return nodes_weights
