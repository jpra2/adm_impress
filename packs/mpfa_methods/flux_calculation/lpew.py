from packs.manager.meshmanager import MeshProperty
from packs.mpfa_methods.weight_interpolation.lpew import preprocess
from packs.mpfa_methods.test.test_monophasic_lsds_dong_paper import calculate_h_dist, calculate_areas
from packs.mpfa_methods.weight_interpolation.lpew import LpewWeight
import numpy as np

class LpewFluxCalculation:

    data_names = ['kn_kt_dsflux']

    def preprocess(self, mesh_properties: MeshProperty):
        """
         O no B2 esta a esquerda da normal do edge e o no B1 a direita
         A face B esta a esquerda da do edge e a face A a direita
         B2 - B; B1 - A;
        """
        preprocess(mesh_properties)
        calculate_h_dist(mesh_properties)
        calculate_areas(mesh_properties)
    
    def create_kn_and_kt(
            self,
            edges,
            nodes_of_edges,
            nodes_centroids,
            adjacencies,
            permeability,
            bool_boundary_edges,
            **kwargs
    ):
        R_matrix = LpewWeight.get_Rmatrix()

        edges_ids = []
        faces_ids = []
        all_kn = []
        all_kt = []

        for edge in edges[bool_boundary_edges]:
            nodes_of_edge = nodes_of_edges[edge]
            faces_adj = adjacencies[edge]
            centroids_nodes_edge = nodes_centroids[nodes_of_edge]
            b1b2 = (centroids_nodes_edge[0] - centroids_nodes_edge[1]).reshape((2, 1))
            rb1b2 = R_matrix.dot(b1b2)
            edge_dim2 = np.linalg.norm(b1b2)**2
            for face in faces_adj:
                perm_face = permeability[face]
                kn = rb1b2.T.dot(perm_face).dot(rb1b2)/edge_dim2
                kt = rb1b2.T.dot(perm_face).dot(b1b2)/edge_dim2
                edges_ids.append(edge)
                faces_ids.append(face)
                all_kn.append(kn[0, 0])
                all_kt.append(kt[0, 0])







        
def create_kn_and_kt(mesh_properties: MeshProperty):
    lpewflux = LpewFluxCalculation()

    if not mesh_properties.verify_name_in_data_names(lpewflux.data_names[0]):
        lpewflux.create_kn_and_kt(**mesh_properties.get_all_data())
        mesh_properties.export_data()



    
    




def get_xi_params_lpew_flux(mesh_properties: MeshProperty):
    """
    Na etapa de preprocessamento:
    O no B2 esta a esquerda da normal do edge e o no B1 a direita
    A face B esta a esquerda da do edge e a face A a direita
    """

    lpewflux = LpewFluxCalculation()

    lpewflux.preprocess(mesh_properties)
    create_kn_and_kt(mesh_properties)





