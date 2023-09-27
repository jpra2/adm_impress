from packs.manager.meshmanager import MeshProperty
from packs.mpfa_methods.weight_interpolation.lpew import preprocess
from packs.mpfa_methods.test.test_monophasic_lsds_dong_paper import calculate_h_dist, calculate_areas
from packs.mpfa_methods.weight_interpolation.lpew import LpewWeight
import numpy as np
from packs.manager.boundary_conditions import BoundaryConditions

class DiamondFluxCalculation:

    data_names = ['kn_kt_dsflux', 'kappa_D_dsflux', 'xi_params_dsflux']

    @staticmethod
    def get_edges_dim(nodes_centroids, nodes_of_edges, edges, **kwargs):
        edges_dim = np.linalg.norm(
            nodes_centroids[nodes_of_edges[edges, 0]] - nodes_centroids[nodes_of_edges[edges, 1]],
            axis=1
        )
        return edges_dim

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

        iedges = ~bool_boundary_edges

        edges_dim = self.get_edges_dim(nodes_centroids, nodes_of_edges, edges)

        for edge in edges[iedges]:
            nodes_of_edge = nodes_of_edges[edge]
            faces_adj = adjacencies[edge]
            centroids_nodes_edge = nodes_centroids[nodes_of_edge]
            b1b2 = (centroids_nodes_edge[0] - centroids_nodes_edge[1]).reshape((2, 1))
            rb1b2 = R_matrix.dot(b1b2)
            edge_dim2 = edges_dim[edge]**2
            for face in faces_adj:
                perm_face = permeability[face]
                kn = rb1b2.T.dot(perm_face).dot(rb1b2)/edge_dim2
                kt = rb1b2.T.dot(perm_face).dot(b1b2)/edge_dim2
                edges_ids.append(edge)
                faces_ids.append(face)
                all_kn.append(kn[0, 0])
                all_kt.append(kt[0, 0])
        
        for edge in edges[bool_boundary_edges]:
            nodes_of_edge = nodes_of_edges[edge]
            face = adjacencies[edge, 0]
            centroids_nodes_edge = nodes_centroids[nodes_of_edge]
            b1b2 = (centroids_nodes_edge[0] - centroids_nodes_edge[1]).reshape((2, 1))
            rb1b2 = R_matrix.dot(b1b2)
            edge_dim2 = np.linalg.norm(b1b2)**2
            perm_face = permeability[face]
            kn = rb1b2.T.dot(perm_face).dot(rb1b2)/edge_dim2
            kt = rb1b2.T.dot(perm_face).dot(b1b2)/edge_dim2
            edges_ids.append(edge)
            faces_ids.append(face)
            all_kn.append(kn[0, 0])
            all_kt.append(kt[0, 0])
        
        edges_ids = np.array(edges_ids)
        faces_ids = np.array(faces_ids)
        all_kn = np.array(all_kn)
        all_kt = np.array(all_kt)

        dtype = [('edge_id', np.int), ('face_id', np.int), ('kn', np.float64), ('kt', np.float64)]

        array = np.zeros(len(edges_ids), dtype=dtype)

        array['edge_id'] = edges_ids
        array['face_id'] = faces_ids
        array['kn'] = all_kn
        array['kt'] = all_kt

        return {self.data_names[0]: array}
    
    def create_kappa_and_D(self, kn_kt_dsflux, h_dist, faces_centroids, nodes_of_edges, nodes_centroids, adjacencies, bool_boundary_edges, edges, **kwargs):
        iedges = ~bool_boundary_edges
        
        dtype = [('kappa', np.float64), ('D', np.float64)]
        array = np.zeros(len(edges), dtype=dtype)

        edges_dim = self.get_edges_dim(nodes_centroids, nodes_of_edges, edges)

        for edge in edges[iedges]:
            faces_adj = adjacencies[edge]
            nodes_of_edge = nodes_of_edges[edge]
            centroids_nodes_edge = nodes_centroids[nodes_of_edge]
            ObOa = faces_centroids[faces_adj[1]] - faces_centroids[faces_adj[0]]
            b1b2 = centroids_nodes_edge[0] - centroids_nodes_edge[1]
            edge_dim = edges_dim[edge]

            Kan = kn_kt_dsflux['kn'][(kn_kt_dsflux['edge_id']==edge) & (kn_kt_dsflux['face_id']==faces_adj[1])]
            Kbn = kn_kt_dsflux['kn'][(kn_kt_dsflux['edge_id']==edge) & (kn_kt_dsflux['face_id']==faces_adj[0])]
            Kat = kn_kt_dsflux['kt'][(kn_kt_dsflux['edge_id']==edge) & (kn_kt_dsflux['face_id']==faces_adj[1])]
            Kbt = kn_kt_dsflux['kt'][(kn_kt_dsflux['edge_id']==edge) & (kn_kt_dsflux['face_id']==faces_adj[0])]
            ha = h_dist[edge, 1]
            hb = h_dist[edge, 0]
            
            array[edge]['kappa'] = (Kan*Kbn)/(Kan*hb + Kbn*ha)
            array[edge]['D'] = (b1b2.dot(ObOa))/edge_dim**2 - (1/edge_dim)*(Kbt*hb/Kbn + Kat*ha/Kan)
        
        return {self.data_names[1]: array}
        
    def create_xi_param_dsflux(self, kappa_D_dsflux, nodes_centroids, nodes_of_edges, edges, bool_boundary_edges, faces_centroids, adjacencies, kn_kt_dsflux, h_dist, **kwargs):
        """
        xi_params[:, 0] = B
        xi_param[:, 1] = A
        xi_params[:, 2] = B1
        xi_params[:, 3] = B2
        """

        iedges = ~bool_boundary_edges
        edges_dim = self.get_edges_dim(nodes_centroids, nodes_of_edges, edges)

        xi_params_dsflux = np.zeros((len(edges), 4))

        for edge in edges[iedges]:
            kappa = kappa_D_dsflux['kappa'][edge]
            D = kappa_D_dsflux['D'][edge]
            edge_dim = edges_dim[edge]
            xi_B = +kappa*edge_dim
            xi_A = -xi_B
            xi_B2 = xi_B*D
            xi_B1 = -xi_B2

            xi_params_dsflux[edge][:] = [xi_B, xi_A, xi_B1, xi_B2]
        
        for edge in edges[bool_boundary_edges]:
            face = adjacencies[edge, 0]
            nodes_of_edge = nodes_of_edges[edge]
            centroids_nodes_edge = nodes_centroids[nodes_of_edge]
            b1b2 = centroids_nodes_edge[0] - centroids_nodes_edge[1]
            b2b1 = -b1b2
            b2Ob = faces_centroids[face] - centroids_nodes_edge[1]
            b1Ob = faces_centroids[face] - centroids_nodes_edge[0]
            edge_dim = edges_dim[edge]
            Kbn = kn_kt_dsflux['kn'][(kn_kt_dsflux['edge_id']==edge) & (kn_kt_dsflux['face_id']==face)][0]
            Kbt = kn_kt_dsflux['kt'][(kn_kt_dsflux['edge_id']==edge) & (kn_kt_dsflux['face_id']==face)][0]
            hb = h_dist[edge, 0]
            t0 = Kbn/(hb*edge_dim)
            coef_B1 = -t0*b2Ob.dot(b2b1) + Kbt
            coef_B2 = -t0*b1Ob.dot(b1b2) - Kbt
            coef_B = t0*edge_dim**2

            xi_params_dsflux[edge][:] = [coef_B, 0, coef_B1, coef_B2]
        
        return {self.data_names[2]: xi_params_dsflux}











        pass
            
    def mount_problem(self, edges_dim, boundary_conditions: BoundaryConditions, xi_params_dsflux, neumann_weights, nodes_weights, adjacencies, faces, nodes_of_edges, edges, bool_boundary_edges, edges_of_nodes, **kwargs):
        
        
        n_faces = faces.shape[0]
        source = np.zeros(n_faces)
        
        lines = []
        cols = []
        data = []
        
        dirichet_nodes = boundary_conditions['dirichlet_nodes']['id']
        values_dirichlet = boundary_conditions['dirichlet_nodes']['value']
        
        edges_of_dirichlet_nodes = np.unique(np.concatenate(edges_of_nodes[dirichet_nodes]))
        neumann_edges = boundary_conditions['neumann_edges']['id']
        edges_of_dirichlet_nodes = np.setdiff1d(edges_of_dirichlet_nodes, neumann_edges)
        
        
        
        import pdb; pdb.set_trace()
        
        ## adicionando o fluxo prescrito nos edges
        for edge, value in zip(boundary_conditions['neumann_edges']['id'], boundary_conditions['neumann_edges']['value']):
            face_adj = adjacencies[edge, 0]
            edge_dim = edges_dim[edge]
            source[face_adj] += value*edge_dim
        
        for edge, value in zip(boundary_conditions['dirichlet_edges']['id'], boundary_conditions['dirichlet_edges']['value']):
            nodes_edge = nodes_of_edges[edge]
            pass
        
        
        
        
        
        

        


      
def create_kn_and_kt(mesh_properties: MeshProperty):
    lpewflux = DiamondFluxCalculation()

    if not mesh_properties.verify_name_in_data_names(lpewflux.data_names[0]):
        resp = lpewflux.create_kn_and_kt(**mesh_properties.get_all_data())
        mesh_properties.insert_data(resp)
        mesh_properties.export_data()

def create_kappa_and_D(mesh_properties: MeshProperty):
    lpewflux = DiamondFluxCalculation()

    if not mesh_properties.verify_name_in_data_names(lpewflux.data_names[1]):
        resp = lpewflux.create_kappa_and_D(**mesh_properties.get_all_data())
        mesh_properties.insert_data(resp)
        mesh_properties.export_data()

def create_xi_param_dsflux(mesh_properties: MeshProperty):
    lpewflux = DiamondFluxCalculation()

    if not mesh_properties.verify_name_in_data_names(lpewflux.data_names[2]):
        resp = lpewflux.create_xi_param_dsflux(**mesh_properties.get_all_data())
        mesh_properties.insert_data(resp)
        mesh_properties.export_data()

def get_xi_params_ds_flux(mesh_properties: MeshProperty):
    """
    Na etapa de preprocessamento:
    O no B2 esta a esquerda da normal do edge e o no B1 a direita
    A face B esta a esquerda da do edge e a face A a direita
    """

    lpewflux = DiamondFluxCalculation()

    lpewflux.preprocess(mesh_properties)
    create_kn_and_kt(mesh_properties)
    create_kappa_and_D(mesh_properties)
    create_xi_param_dsflux(mesh_properties)





