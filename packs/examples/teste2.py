from packs import defpaths
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from ..directories import data_loaded
import numpy as np
import os
import scipy.sparse as sp

'Rename function as prep_2D'

def load_mesh_2d(mesh_name):
    mesh_name2 = mesh_name.split('/')
    mesh_name3 = mesh_name2[-1].split('.')[0]
    fine_mesh_path = mesh_name #os.path.join(rel_path, mesh_name)
    fine_mesh_properties_name = mesh_name3 
    
    fine_mesh_path_v4 = fine_mesh_path

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)
    M = Mesh(fine_properties) #convert mesh 2D properties to 3D language type
    return M

class Mesh():
    global z
    z = data_loaded['z_2D']
            
    def __init__(self, mesh):
        self.volumes = self.Volumes(mesh)
        self.faces = self.Faces(mesh)
        self.edges = self.Edges(mesh)
        self.nodes = self.Nodes(mesh)
        self.bool_boundary_edges = mesh.bool_boundary_edges
        self.adjacencies = mesh.adjacencies
        #self.data = self.Data()
        self.areas = mesh.areas
        #import pdb;pdb.set_trace()
    
    class Volumes():
        def __init__(self, mesh):
            self.all = mesh.faces
            self.adjacencies = mesh.adjacencies
            self.volumes_centroids = mesh.faces_centroids
            self.areas = mesh.areas

        def center(self, arg):
            #import pdb; pdb.set_trace()
            centroids = np.ones((len(self.volumes_centroids[:,0]),3))*(-z/2)
            centroids[:,[0,1]] = self.volumes_centroids
            return centroids

        def bridge_adjacencies(self, mesh, arg2, arg3):
            lines = np.array([self.adjacencies[:, 0], self.adjacencies[:, 1], self.adjacencies[:, 0], self.adjacencies[:, 1]]).flatten()
            cols = np.array([self.adjacencies[:, 1], self.adjacencies[:, 0], self.adjacencies[:, 0], self.adjacencies[:, 1]]).flatten()
            data = np.array([np.ones(len(self.adjacencies[:, 0])), np.ones(len(self.adjacencies[:, 0])),
                            np.zeros(len(self.adjacencies[:, 0])), np.zeros(len(self.adjacencies[:, 0]))]).flatten()
            all_neig = sp.csc_matrix((data, (lines, cols)), shape = (n_volumes, n_volumes)).toarray()
            all_neig = all_neig.astype(int)
            all_neig2 = all_neig + np.identity(n_volumes)
            allneig_and_vol = all_neig2.astype(int)
            print("IN")
            import pdb; pdb.set_trace()
            return self.adjacencies #são as faces que compõem o volume (ou seja, as edges, preciso CORRIGIR, não é esse vetor - ver com JOÃO)
        
        def volume(self, vols):
            return self.areas[vols]* z
            
    class Faces():
        def __init__(self, mesh):
            self.all = mesh.edges
            self.bool_boundary_edges = mesh.bool_boundary_edges
            self.internal = mesh.edges[~mesh.bool_boundary_edges]
            self.normal = mesh.unitary_normal_edges
            self.adjacencies = mesh.adjacencies
            self.areas = mesh.areas
            self.boundary = mesh.edges[mesh.bool_boundary_edges]
            self.faces = mesh.edges      
            self.faces_centroids = mesh.edges_centroids  
            self.internal_faces_unitary_normal = np.empty((len(self.internal),3))
            self.internal_faces_unitary_normal[:,[0,1]] = mesh.unitary_normal_edges[~mesh.bool_boundary_edges]    
            self.internal_faces_unitary_normal[...,-1] = 0
            self.faces_of_volumes = mesh.edges_of_faces

            self.normal = np.empty((len(self.faces),3))
            self.normal[:,[0,1]] = mesh.unitary_normal_edges    
            self.normal[abs(self.normal)<1e-13] = 0
            self.normal[:,-1] = 0
            
            
            #self.normal = self.normal * np.min(abs(self.faces_centroids[:,0]))*2 #correct this later
            
            #import pdb; pdb.set_trace()
        
        def center(self, faces_ids):
            centroids = np.ones((len(self.faces_centroids[:,0]),3))*(-z/2)
            centroids[:,[0,1]] = self.faces_centroids
            return centroids[faces_ids]
        
        def bridge_adjacencies(self, faces, arg2, arg3):
            if arg3 == 3:
                adjacencies= self.adjacencies[faces]
            else:
                'I did not tried this with non quadrilateral meshes!!!'
                adjacencies = self.faces_of_volumes #only internal faces
                '''adj0 = self.adjacencies[:,0]
                adj0 = adj0[adj0>=0]
                adj1 = self.adjacencies[:,1]
                adj1 = adj1[adj1>=0]

                edg0 = self.all[self.adjacencies[:,0]>=0]
                edg1 = self.all[self.adjacencies[:,1]>=0]
                lines = np.array([np.concatenate((adj0, adj1))]).flatten()
                cols = np.array([np.concatenate((edg0, edg1))]).flatten()
                data = np.array([np.concatenate((np.ones(len(edg0),dtype=bool),np.ones(len(edg1),dtype=bool)))]).flatten()
                faces_edges_all_bool = sp.csc_matrix((data, (lines, cols)), shape = (len(self.faces),len(self.all))).toarray()
                faces_edges_all_idxs = faces_edges_all_bool * self.all[np.newaxis,:]
                faces_edges = np.zeros([len(self.faces),4])
                faces_edges_bool = np.ones_like(faces_edges,dtype=bool)
                faces_edges[faces_edges_bool] = faces_edges_all_idxs[faces_edges_all_bool] 
                adjacencies= faces_edges.astype(int)'''
            return adjacencies

        def _connectivities(self, faces_ids):
            pass
            
            
        def area(self, faces):
            return self.areas
            
    class Edges():
        def __init__(self, mesh):
            self.all = mesh.edges
            self.nodes_of_edges = mesh.nodes_of_edges
        
        def center(self, faces_centroids):
            return faces_centroids #see this later, its supposed to be edges_centroids
        
        def bridge_adjacencies(self, mesh, arg2, arg3):
            return self.nodes_of_edges
    
    class Nodes():
        def __init__(self, mesh):
            self.all = mesh.nodes
        
        def center(self, nodes_centroids):
            return nodes_centroids 
    
    