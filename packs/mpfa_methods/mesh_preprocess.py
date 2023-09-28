from packs.manager.meshmanager import MeshProperty
from packs.utils import calculate_face_properties
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists
import numpy as np
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation


class MpfaPreprocess:
    
    def define_A_B_points_of_edges(self, mesh_properties: MeshProperty):
        """define os pontos A e B de cada edge da malha
            o ponto B deve estar a esquerda do vetor normal do edge
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
        
        resp = mesh_properties.nodes_of_edges.copy()

        theta = np.pi/2
        cossin = np.array([np.cos(theta), np.sin(theta)])
        R = np.array([
            [cossin[0], -cossin[1]],
            [cossin[1],  cossin[0]]
        ])

        for edge in mesh_properties.edges:
            B = mesh_properties.nodes_of_edges[edge, 0]
            A = mesh_properties.nodes_of_edges[edge, 1]
            AB = mesh_properties.nodes_centroids[B] - mesh_properties.nodes_centroids[A]
            unitary_normal_vector_rotated = R.dot(mesh_properties.unitary_normal_edges[edge])
            proj = AB.dot(unitary_normal_vector_rotated)
            if proj >= 0:
                pass
            else:
                resp[edge,:] = [A, B]
        
        mesh_properties.insert_or_update_data({'nodes_of_edges': resp})
        mesh_properties.export_data()
        
    def calculate_areas(self, mesh_properties: MeshProperty):
        

        if not mesh_properties.verify_name_in_data_names('areas'):

            centroids_nodes = mesh_properties.nodes_centroids
            z_centroids = np.zeros((len(centroids_nodes), 1))
            centroids_nodes = np.hstack([centroids_nodes, z_centroids])
            nodes_of_faces = mesh_properties.nodes_of_faces
            cnodes_faces = centroids_nodes[nodes_of_faces]
            n_faces = len(mesh_properties.faces)

            areas = np.zeros(n_faces)
            for i in range(n_faces):
                areas[i] = calculate_face_properties.polygon_area(cnodes_faces[i])
            mesh_properties.insert_data({'areas': areas})
            mesh_properties.export_data()

    def calculate_h_dist(self, mesh_properties: MeshProperty):
        if not mesh_properties.verify_name_in_data_names('h_dist'):
            
            h_dist = calculate_face_properties.create_face_to_edge_distances(
                mesh_properties.faces_centroids,
                mesh_properties.adjacencies,
                mesh_properties.nodes_of_edges,
                mesh_properties.edges,
                mesh_properties.nodes_centroids,
                mesh_properties.bool_boundary_edges
            )
            mesh_properties.insert_data({'h_dist': h_dist})

            bedges = mesh_properties.bool_boundary_edges
            iedges = ~bedges

            m1 = np.mean(h_dist[iedges, 0])
            m2 = np.mean(h_dist[iedges, 1])
            m3 = np.mean(h_dist[bedges, 0])
            m_hdist = np.mean([m1, m2, m3])
            mesh_properties.insert_data({'m_hdist': np.array([m_hdist])})
            mesh_properties.export_data()

    def create_properties_if_not_exists(self, mesh_name, mesh_properties_name):
        create_properties_if_not_exists(mesh_name, mesh_properties_name)
    
    def test_unitary_normal_edges(self, mesh_properties: MeshProperty):
        
        faces_centroids = mesh_properties.faces_centroids
        nodes_of_edges = mesh_properties.nodes_of_edges
        nodes_centroids = mesh_properties.nodes_centroids
        adjacencies = mesh_properties.adjacencies
        edges = mesh_properties.edges
        unitary_normal_edges = mesh_properties.unitary_normal_edges
        
        edges_centroids = (nodes_centroids[nodes_of_edges[:, 1]] + nodes_centroids[nodes_of_edges[:, 0]])/2 
        
        for edge in edges:
            unitary_normal = unitary_normal_edges[edge]
            edge_centroid = edges_centroids[edge]
            face_centroid = faces_centroids[adjacencies[edge, 0]]
            direction = edge_centroid - face_centroid
            test = np.dot(direction, unitary_normal)
            if test >=0:
                pass
            else:
                raise ValueError('Os vetores normais dos edges devem estar apontando da face da esquerda do vetor de adjacencias para o centroide do edge')
            




