from packs.manager.meshmanager import MeshProperty
from packs.utils import calculate_face_properties
from packs.mpfa_methods.weight_interpolation.test.test_gls_weights import create_properties_if_not_exists


class MpfaPreprocess:

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
        pass




