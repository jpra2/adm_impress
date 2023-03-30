from packs.manager.arraydatamanager import SuperArrayManager
import numpy as np

class CalculateGlsWeight2D:

    def __init__(self, **kwargs):
        """ Initialize the class

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
        
        for node in self.nodes[self.bool_internal_nodes]:
            nodes_adj = self.nodes_of_nodes[node]
            edges_adj = self.edges_of_nodes[node]
            centroid_node = self.nodes_centroids[node]
            centroids_nodes_adj = self.nodes_centroids[nodes_adj]
            faces_adj = self.faces_of_nodes[node]
            adj_edges_adj = self.adjacencies[edges_adj]
            centroids_faces_adj = self.faces_centroids[faces_adj]
            
            #################################################
            n_mfaces = len(faces_adj)
            mfaces = np.zeros((n_mfaces, 2*n_mfaces))
            f_dists = centroids_faces_adj - centroid_node
            for i in range(n_mfaces):
                mfaces[i,2*i:2*i+2] = f_dists[i]
            #################################################
            
            #################################################
            n_mnodes = len(nodes_adj)
            mnodes = np.zeros((n_mnodes, 2*n_mnodes))
            nodes_dists = centroids_nodes_adj - centroid_node
            for i in range(n_mnodes):
                mnodes[i,2*i:2*i+2] = nodes_dists[i]
                
            mnodes[0, 2*n_mnodes-2:] = -nodes_dists[0]
            
            for i in range(1, n_mnodes):
                mnodes[i, 2*i-2:2*i] = -nodes_dists[i]
            #################################################
            
            ##################################################
            n_mnormalperm = n_mfaces
            mnormalperm = np.zeros((n_mnormalperm, 2*n_mnormalperm))
            for i in range(n_mnormalperm):
                face = faces_adj[i]
                edge = edges_adj[i]
                faces_adj_edge = self.adjacencies[edge]
                unitary_normal_edge = self.unitary_normal_edges[edge]
                face_perm = self.permeability[face]
                if face == faces_adj_edge[1]:
                    normal = -unitary_normal_edge
                else:
                    normal = unitary_normal_edge
                
                mnormalperm[i, 2*i:2*i+2] = np.dot(normal, face_perm)
            
            for i in range(1, n_mnormalperm):
                face = faces_adj[i]
                edge = edges_adj[i]
                faces_adj_edge = self.adjacencies[edge]
                unitary_normal_edge = self.unitary_normal_edges[edge]
                face_perm = self.permeability[faces_adj_edge[faces_adj_edge!=face][0]]
                if face == faces_adj_edge[1]:
                    normal = -unitary_normal_edge
                else:
                    normal = unitary_normal_edge
                
                mnormalperm[i, 2*i-2:2*i] = -np.dot(normal, face_perm)
            
            i = 0
            face = faces_adj[i]
            edge = edges_adj[i]
            faces_adj_edge = self.adjacencies[edge]
            unitary_normal_edge = self.unitary_normal_edges[edge]
            face_perm = self.permeability[faces_adj_edge[faces_adj_edge!=face][0]]
            if face == faces_adj_edge[1]:
                normal = -unitary_normal_edge
            else:
                normal = unitary_normal_edge
            
            mnormalperm[0, 2*n_mnormalperm-2:] = -np.dot(normal, face_perm)
            ########################################################

            ###################################################
            Mv = np.zeros((3*n_mfaces, 2*n_mfaces+1))
            Mv[0:n_mfaces,-1] = 1
            Mv[0:n_mfaces, 0:2*n_mfaces] = mfaces
            Mv[n_mfaces:2*n_mfaces, 0:2*n_mfaces] = mnodes
            Mv[2*n_mfaces:3*n_mfaces, 0:2*n_mfaces] = mnormalperm
            ###################################################

            #####################################
            Nv = np.zeros((3*n_mfaces, n_mfaces))
            Nv[np.arange(n_mfaces), np.arange(n_mfaces)] = 1
            #####################################

            ##################################
            V = np.linalg.inv(np.dot(Mv, Mv.T))
            ##################################



            
            import pdb; pdb.set_trace()



                


            
                
            
            import pdb; pdb.set_trace()








