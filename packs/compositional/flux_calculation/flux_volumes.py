import numpy as np
import scipy.sparse as sp
from packs.utils import constants as ctes

class Flux:
    """ Class created for computing flux accondingly to the First Order Upwind \
    Method - actually, it's only Flux because of the choice of mobilities and \
    other properties to be taken as function of the potencial gradient. But here \
    flux is computed at the interfaces, and then, volumes, only that."""

    def update_flux(self, M, fprop, Ft_internal_faces, rho_j_internal_faces,
        mobilities_internal_faces):
        ''' Main function that calls others '''

        self.Nk = fprop.Nk

        Fj_internal_faces = self.update_Fj_internal_faces(Ft_internal_faces,
            rho_j_internal_faces, mobilities_internal_faces, fprop.Pcap[:,ctes.v0],
            ctes.z[ctes.v0], ctes.pretransmissibility_internal_faces)
        Fk_internal_faces = self.update_Fk_internal_faces(
            fprop.xkj_internal_faces, fprop.Csi_j_internal_faces, Fj_internal_faces)
        Fk_vols_total = self.update_flux_volumes(Fk_internal_faces)
        return Fk_vols_total

    def update_Fj_internal_faces(self, Ft_internal_faces, rho_j_internal_faces,
        mobilities_internal_faces, Pcap_face, z_face,
        pretransmissibility_internal_faces):
        ''' Function to calculate phase flux '''

        frj = mobilities_internal_faces[0,...] / \
            np.sum(mobilities_internal_faces[0,...], axis = 0)

        Fj_internal_faces = frj[np.newaxis,...] * (Ft_internal_faces +
            pretransmissibility_internal_faces * (np.sum(mobilities_internal_faces *
            (Pcap_face[:,:,1] - Pcap_face[:,:,0] - ctes.g * rho_j_internal_faces *
            (z_face[:,1] - z_face[:,0])), axis=1) - np.sum(mobilities_internal_faces,
            axis=1) * (Pcap_face[:,:,1] - Pcap_face[:,:,0] - ctes.g *
            rho_j_internal_faces * (z_face[:,1] - z_face[:,0]))))
        return Fj_internal_faces

    def update_Fk_internal_faces(self, xkj_internal_faces, Csi_j_internal_faces, Fj_internal_faces):
        ''' Function to compute component molar flux '''
        Fk_internal_faces = np.sum(xkj_internal_faces * Csi_j_internal_faces *
            Fj_internal_faces, axis = 1)
        return Fk_internal_faces

    def update_flux_volumes(self, Fk_internal_faces):
        ''' Function to compute component molar flux balance through the control \
        volume interfaces '''
        cx = np.arange(ctes.n_components)

        Fk_internal_faces *= ctes.N_sign
        lines = np.array([np.repeat(cx,len(ctes.in_vols_pairs[:,0])), np.repeat(cx,len(ctes.in_vols_pairs[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(ctes.in_vols_pairs[:,0],ctes.n_components), np.tile(ctes.in_vols_pairs[:,1], ctes.n_components)]).flatten()
        data = np.array([-Fk_internal_faces, Fk_internal_faces]).flatten()
        Fk_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_components, ctes.n_volumes)).toarray()

        #lines = np.array([np.repeat(cx,len(ctes.v0[:,0])), np.repeat(cx,len(ctes.v0[:,1]))]).astype(int).flatten()
        #cols = np.array([np.tile(ctes.v0[:,0],ctes.n_components), np.tile(ctes.v0[:,1], ctes.n_components)]).flatten()
        #data = np.array([-Fk_internal_faces, Fk_internal_faces]).flatten()
        #Fk_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_components, ctes.n_volumes)).toarray()

        '''
        inds = np.array([0,-1])
        a = (self.Nk[:,inds] - abs(self.Nk[:,inds]))/2
        a = a[:,self.v0]
        if a>=0:
            Fk_contour_face = self.Nk[0,-1] **2/2*len(self.Nk[-1])
        else:
            Fk_contour_face = self.Nk[0,0]**2/2*len(self.Nk[0])
        Fk_contour_face = self.Nk[0,-1]**2/2*len(self.Nk[0])
        Fk_vols_total[:,0] = -(Fk_internal_faces[:,0] - Fk_contour_face)
        Fk_vols_total[:,-1] = -(Fk_contour_face - Fk_internal_faces[:,1])'''

        return Fk_vols_total
