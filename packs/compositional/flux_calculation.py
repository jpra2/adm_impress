import numpy as np
from ..directories import data_loaded
from ..utils import constants as ctes
import scipy.sparse as sp

class Flux:

    def update_flux(self, fprop, kprop, total_flux_internal_faces):
        self.update_phase_flux_internal_faces(fprop, total_flux_internal_faces)
        self.update_flux_volumes(fprop, kprop)

    def update_phase_flux_internal_faces(self, fprop, total_flux_internal_faces):
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        #Pot_hid = fprop.P + fprop.Pcap
        #Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        #Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]

        frj = fprop.mobilities_internal_faces[0,:,:] / np.sum(fprop.mobilities_internal_faces[0,:,:], axis = 0)

        self.phase_flux_internal_faces = frj * (total_flux_internal_faces + (np.sum(fprop.mobilities_internal_faces *
                                        ctes.pretransmissibility_internal_faces * (fprop.Pcap[:,ctes.v0[:,1]] -
                                        fprop.Pcap[:,ctes.v0[:,0]] - ctes.g * fprop.phase_densities_internal_faces *
                                        (z_up - z)), axis=1) - np.sum(fprop.mobilities_internal_faces *
                                         ctes.pretransmissibility_internal_faces, axis=1) * (fprop.Pcap[:,ctes.v0[:,1]] -
                                         fprop.Pcap[:,ctes.v0[:,0]] - ctes.g *fprop.phase_densities_internal_faces * (z_up - z))))

        #phase_flux_internal_faces = - (fprop.mobilities_internal_faces * ctes.pretransmissibility_internal_faces
        #                             * (Pot_hidj_up - Pot_hidj - ctes.g * fprop.phase_densities_internal_faces
        #                             * (z_up - z)))
        # M.flux_faces[M.faces.internal] = total_flux_internal_faces * M.faces.normal[M.faces.internal].T

    def update_flux_volumes(self, fprop, kprop):
        component_flux_internal_faces = np.sum(fprop.component_molar_fractions_internal_faces * fprop.phase_molar_densities_internal_faces *
                                self.phase_flux_internal_faces, axis = 1)
        cx = np.arange(kprop.n_components)
        lines = np.array([np.repeat(cx,len(ctes.v0[:,0])), np.repeat(cx,len(ctes.v0[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(ctes.v0[:,0],kprop.n_components), np.tile(ctes.v0[:,1], kprop.n_components)]).flatten()
        data = np.array([-component_flux_internal_faces, component_flux_internal_faces]).flatten()
        fprop.component_flux_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (kprop.n_components, ctes.n_volumes)).toarray()
