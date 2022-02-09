import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from scipy import linalg
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

class TPFASolver:
    def __init__(self, dVjdNk, dVjdP):
        self.dVt_derivatives(dVjdNk, dVjdP)

    def get_pressure(self, M, wells, fprop, Pold, delta_t):
        T = self.update_transmissibility(M, wells, fprop, delta_t)
        D = self.update_independent_terms(M, fprop, Pold, wells, delta_t)
        Pnew = self.update_pressure(T, D)
        Ft_internal_faces = self.update_total_flux_internal_faces(M, fprop, Pnew)
        self.update_flux_wells(fprop, Pnew, wells, delta_t)
        return Pnew, Ft_internal_faces, self.q

    def dVt_derivatives(self, dVjdNk, dVjdP):
        self.dVtP = dVjdP.sum(axis=1)[0]
        self.dVtk = dVjdNk.sum(axis=1)

    def update_transmissibility(self, M, wells, fprop, delta_t):
        self.t0_internal_faces_prod = fprop.xkj_internal_faces * \
                                      fprop.Csi_j_internal_faces * \
                                      fprop.mobilities_internal_faces

        ''' Transmissibility '''
        t0 = (self.t0_internal_faces_prod).sum(axis = 1)
        t0 = t0 * ctes.pretransmissibility_internal_faces
        T = sp.csr_matrix((ctes.n_volumes, ctes.n_volumes))
        # Look for a way of doing this not using a loop!!!
        for i in range(ctes.n_components):
            lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
            cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
            data = np.array([-t0[i,:], -t0[i,:], +t0[i,:], +t0[i,:]]).flatten()

            Ta = (sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes)))#.toarray()
            T += Ta.multiply(self.dVtk[i, :, np.newaxis])

        T *= delta_t
        ''' Transmissibility diagonal term '''
        #diag = np.diag((ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP))
        #T += diag
        T.setdiag(T.diagonal().flatten() + (ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP))
        self.T_noCC = np.copy(T.toarray())#    .toarray() #Transmissibility without contour conditions

        ''' Includding contour conditions '''
        T[wells['ws_p'],:] = 0
        T[wells['ws_p'], wells['ws_p']] = 1
        return T

    def pressure_independent_term(self, fprop, Pold):
        vector = ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP
        pressure_term = vector * Pold
        return pressure_term

    def capillary_and_gravity_independent_term(self, fprop):

        t0_j = self.t0_internal_faces_prod * \
            ctes.pretransmissibility_internal_faces
        t0_k = ctes.g * np.sum(fprop.rho_j_internal_faces * t0_j, axis=1)

        # Look for a better way to do this !!!
        cap = np.zeros([ctes.n_volumes])
        grav = np.zeros([ctes.n_volumes,ctes.n_volumes])
        if any((ctes.z - ctes.z[0]) != 0):
            for i in range(ctes.n_components):
                lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                data = np.array([t0_k[i,:], t0_k[i,:], -t0_k[i,:], -t0_k[i,:]]).flatten()
                t0_rho = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes)).toarray()
                grav += t0_rho * self.dVtk[i,:]

                for j in range(ctes.n_phases):
                    lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                    cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                    data = np.array([t0_j[i,j,:], t0_j[i,j,:], -t0_j[i,j,:], -t0_j[i,j,:]]).flatten()
                    t0 = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes))*self.dVtk[i,:]
                    cap += t0 @ fprop.Pcap[j,:]

        gravity_term = grav @ ctes.z
        #import pdb; pdb.set_trace()
        # capillary_term = np.sum(self.dVtk * np.sum (fprop.xkj *
        #         fprop.Csi_j * fprop.mobilities * fprop.Pcap, axis = 1), axis = 0)
        return cap, gravity_term

    def volume_discrepancy_independent_term(self, fprop):
        volume_discrepancy_term = fprop.Vp - fprop.Vt
        if np.max(abs(volume_discrepancy_term)) > 5e-4:
            #import pdb; pdb.set_trace()
            print('hit: ', np.max(abs(volume_discrepancy_term)))
        return volume_discrepancy_term

    def well_term(self, fprop, wells):
        self.q = np.zeros([ctes.n_components, ctes.n_volumes])
        well_term = np.zeros(ctes.n_volumes)
        if len(wells['ws_q']) > 0:
            self.q[:,wells['ws_q']] =  wells['values_q']
            well_term[wells['ws_q']] = np.sum(self.dVtk[:,wells['ws_q']] *
                self.q[:,wells['ws_q']], axis = 0)
        return well_term

    def update_independent_terms(self, M, fprop, Pold, wells, delta_t):
        self.pressure_term = self.pressure_independent_term(fprop, Pold)
        self.capillary_term, self.gravity_term = self.capillary_and_gravity_independent_term(fprop)
        self.volume_term = self.volume_discrepancy_independent_term(fprop)
        well_term = self.well_term(fprop, wells)
        independent_terms = self.pressure_term - self.volume_term  + delta_t * \
            well_term - delta_t * (self.capillary_term + self.gravity_term)
        independent_terms[wells['ws_p']] = wells['values_p'] + ctes.g * \
            fprop.rho_j[0,0,wells['ws_p']] * (ctes.z[wells['ws_p']] - ctes.z[ctes.bhp_ind])
        return independent_terms

    def update_pressure(self, T, D):
        #P = linalg.solve(T,D)
        P = spsolve(T,D)
        return P

    def update_total_flux_internal_faces(self, M, fprop, Pnew):
        Pot_hid = Pnew + fprop.Pcap
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        Ft_internal_faces = - np.sum(fprop.mobilities_internal_faces
            * ctes.pretransmissibility_internal_faces * ((Pot_hidj_up - Pot_hidj) -
            ctes.g * fprop.rho_j_internal_faces * (z_up - z)), axis = 1)
        return Ft_internal_faces

    def update_flux_wells(self, fprop, Pnew, wells, delta_t):
        wp = wells['ws_p']

        if len(wp)>=1:
            #if Pnew[0]<Pnew[1]: import pdb; pdb.set_trace()
            well_term =  (self.T_noCC[wp,:] @ Pnew - self.pressure_term[wp] +
                self.volume_term[wp]) / delta_t  + self.capillary_term[wp] + \
                self.gravity_term[wp]
            mob_ratio = fprop.mobilities[:,:,wp] / np.sum(fprop.mobilities[:,:,wp], axis = 1)

            q_term = fprop.xkj[:,:,wp] * mob_ratio * fprop.Csi_j[:,:,wp]

            ws_p_inj = np.argwhere(wells['ws_p']==wells['ws_inj']).flatten()
            q_term[...,ws_p_inj] = wells['inj_p_term']
            #import pdb; pdb.set_trace()

            self.q[:,wp] = np.sum(q_term * well_term, axis = 1)
            fprop.qk_prod = self.q[:,wells['ws_prod']]
            fprop.q_phase = mob_ratio * well_term
