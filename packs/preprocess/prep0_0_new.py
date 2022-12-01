from .. import directories as direc
from . import directories_impress as direc_impress
from ..directories import data_loaded
from ..utils.utils_old import get_box
from math import pi
import numpy as np
from .prep0_0 import Preprocess0
# from .preprocess1 import set_saturation_regions

class PreprocessUnfied(Preprocess0):
    '''
    objetivo: setar k_harmonico, pretransmissibilidade, transmissibilidade, area, dist_centroid, u_normal nas faces
    permeabilidade, volume, centroide, porosidade nos volumes
    '''

    def check_struct(self, M):
        faces = M.faces.all
        n_faces = len(faces)
        areas = M.faces.area(faces)
        edges_nodes = M.edges.bridge_adjacencies(M.edges.all,1,0)

        # Generalizando para faces triangulares e quadrilateras
        faces_edges = M.faces.bridge_adjacencies(faces, 2, 1)

        try:
            faces_quadr = (faces_edges.shape[1]==4)
            return (faces_edges.shape[1]==4)
        except:
            return False

    def set_faces_normals(self, M):
        """ get faces normals """
        #import time
        #t0=time.time()
        faces = M.faces.all
        n_faces = len(faces)
        areas = M.faces.area(faces)
        edges_nodes = M.edges.bridge_adjacencies(M.edges.all,1,0)

        # Generalizando para faces triangulares e quadrilateras
        faces_edges = M.faces.bridge_adjacencies(faces, 2, 1)
        #t1 = time.time()
        #print(t1-t0)
        #import pdb; pdb.set_trace()
        faces_quadr = np.zeros_like(faces,dtype=bool)
        list_faces_edges = faces_edges.tolist()
        max_n_edges_per_face = max(map(len,list_faces_edges))
        faces_edges = np.ones((len(faces),max_n_edges_per_face),dtype=int)*-1
        i=0
        for element in list_faces_edges:
            stop = len(element)
            faces_edges[i,0:stop] = element
            faces_quadr[i] = (len(element)==4)
            i+=1

        nodes_coord = M.nodes.center(M.nodes.all)
        edges_vectors = abs(nodes_coord[edges_nodes[:,1]]-nodes_coord[edges_nodes[:,0]])
        faces_vectors = edges_vectors[faces_edges]

        prod1 = abs(np.cross(faces_vectors[:,0],faces_vectors[:,1],axis=-1))
        prod2 = abs(np.cross(faces_vectors[:,0],faces_vectors[:,2],axis=-1))
        normals = np.copy(prod1)
        normals[prod1.sum(axis=-1)==0] = prod2[prod1.sum(axis=-1)==0]
        normals/=np.linalg.norm(normals,axis=-1)[:,np.newaxis]

        #faces_vecs = np.empty((n_faces,max_n_edges_per_face,3),dtype=float)
        #faces_vecs[faces_quadr] = edges_hs_by_axis[faces_edges[faces_quadr]]

        #faces_hs = np.empty_like(faces_edges,dtype=float)
        #faces_hs[faces_quadr] = edges_hs[faces_edges[faces_quadr]]
        #any_quad = np.any(faces_quadr)
        #expr=(-1)*any_quad+max_n_edges_per_face #auxiliary expression to make this automatized when there is no quadrangle
        #faces_hs[~faces_quadr,:expr] = edges_hs[faces_edges[~faces_quadr,:expr]]
        #faces_hs[~faces_quadr,expr]=-1 #não sei para que vai ser util ainda...
        return normals

    def set_area_unstruct(self,M):
        faces = M.faces.all
        n_faces = len(faces)
        areas = M.faces.area(faces)
        cent_vols = M.volumes.center(M.volumes.all)
        faces_nodes = M.faces.bridge_adjacencies(faces,2,0)

        faces_quadr = np.zeros_like(faces,dtype=bool)
        list_faces_nodes = faces_nodes.tolist()
        max_nodes_per_face = max(map(len,list_faces_nodes))
        faces_nodes = np.ones((len(faces),max_nodes_per_face),dtype=int)*-1
        i=0
        for element in list_faces_nodes:
            stop = len(element)
            faces_nodes[i,0:stop] = element
            faces_quadr[i] = (len(element)==4)
            i+=1

        faces_node_point = M.nodes.center(faces_nodes[:,0])

        vols_neig_internal_faces = M.faces.bridge_adjacencies(M.faces.internal,2,3)
        vols_neig_boundary_faces = M.faces.bridge_adjacencies(M.faces.boundary,2,3)
        normals_faces = self.set_faces_normals(M)

        d_fpoint_cent_internal = abs(cent_vols[vols_neig_internal_faces] - faces_node_point[M.faces.internal][:,np.newaxis,:])
        d_fpoint_cent_boundary = abs(cent_vols[vols_neig_boundary_faces][:,0,:] - faces_node_point[M.faces.boundary])

        dist_cent = np.empty((n_faces,2))
        dist_cent[M.faces.internal,:] = np.sum(d_fpoint_cent_internal*abs(normals_faces[M.faces.internal][:,np.newaxis,:]),axis=-1)
        dist_cent[M.faces.boundary,0] = np.sum(d_fpoint_cent_boundary*abs(normals_faces[M.faces.boundary]),axis=-1)
        dist_cent[M.faces.boundary,1] = dist_cent[M.faces.boundary,0] #?? not sure about this
        dist_between_centers = np.sum(dist_cent,axis=1) #isso é a soma dos hs "h_face_cent_vol"

        M.data[M.data.variables_impress['area']] = areas
        M.data[M.data.variables_impress['dist_cent']] = dist_between_centers
        M.data[M.data.variables_impress['volume']] = M.volumes.volume(M.volumes.all)
        M.data[M.data.variables_impress['NODES']] = M.nodes.center(M.nodes.all)
        M.data['h_face_cent_vol'] = dist_cent #util para o MPFA-D
        M.data['faces_nodes'] = faces_nodes
        M.data['faces_normals'] = normals_faces
        #M.data[M.data.variables_impress['hs']] = faces_hs
        #t1 = time.time()
        #print(t1-t0)
        #import pdb; pdb.set_trace()


    def run(self, M):
        self.update_centroids_and_unormal(M)
        self.set_permeability_and_phi(M)
        self.set_area_unstruct(M)
        #import pdb; pdb.set_trace()
        '''if self.check_struct(M):
            self.set_area_hex_structured(M)
            self.set_k_harm_hex_structured(M)
            self.set_pretransmissibility(M)
            self.set_transmissibility_monophasic(M)
            self.initial_gama(M)
        else: self.set_area_unstruct(M)'''

        #self.set_pretransmissibility(M)
        #self.set_transmissibility_monophasic(M)
        #self.initial_gama(M)
