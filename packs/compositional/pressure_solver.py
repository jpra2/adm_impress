import numpy as np
from ..directories import data_loaded
from ..utils import constants as ctes
from ..solvers.solvers_scipy.solver_sp import SolverSp
from scipy import linalg
import scipy.sparse as sp
from .partial_derivatives_new import PartialDerivatives

class TPFASolver:

    def get_pressure(self, M, wells, fprop, kprop, delta_t,r):
        if r==0.8: self.dVt_derivatives(fprop, kprop)
        T = self.update_transmissibility(M, wells, fprop, kprop, delta_t)
        D = self.update_independent_terms(M, fprop, kprop, wells, delta_t)
        self.update_pressure(T, D, fprop)
        self.update_total_flux_internal_faces(M, fprop, kprop)
        self.update_flux_wells(fprop, kprop, wells, delta_t)
        return self.P, self.total_flux_internal_faces, self.q

    def dVt_derivatives(self, fprop, kprop):
        self.dVtk = np.zeros([kprop.n_components, ctes.n_volumes])
        if kprop.load_k:

            if not kprop.compressible_k:
                dVtP = np.zeros(ctes.n_volumes)
                self.dVtk[0:kprop.Nc,:] = 1 / fprop.phase_molar_densities[0,0,:]
            else: self.dVtk[0:kprop.Nc,:], dVtP = PartialDerivatives(fprop).get_all_derivatives(kprop, fprop)

        else: dVtP = np.zeros(ctes.n_volumes)

        if kprop.load_w:

            self.dVtk[kprop.n_components-1,:] = 1 / fprop.phase_molar_densities[:,kprop.n_phases-1,:]
            dVwP = -fprop.component_mole_numbers[kprop.Nc,:] * fprop.ksi_W0 * ctes.Cw / \
                        fprop.phase_molar_densities[0,kprop.n_phases-1,:]**2

        else: dVwP = np.zeros(ctes.n_volumes)

        self.dVtP = dVtP + dVwP


    def update_transmissibility(self, M, wells, fprop, kprop, delta_t):
        self.t0_internal_faces_prod = fprop.component_molar_fractions_internal_faces * fprop.phase_molar_densities_internal_faces \
                                * fprop.mobilities_internal_faces
        ''' Transmissibility '''
        t0 = (self.t0_internal_faces_prod).sum(axis = 1)
        t0 = t0 * ctes.pretransmissibility_internal_faces
        T = np.zeros([ctes.n_volumes, ctes.n_volumes])

        # Look for a way of doing this not using a loop
        for i in range(kprop.n_components):
            lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
            cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
            data = np.array([-t0[i,:], -t0[i,:], +t0[i,:], +t0[i,:]]).flatten()

            Ta = (sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes))).toarray()
            T += Ta * self.dVtk[i,:]

        #T = (T * self.dVtk.T[:,np.newaxis, :]).sum(axis = 2)

        T = T * delta_t
        ''' Transmissibility diagonal term '''
        diag = np.diag((ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP))
        T += diag
        self.T_noCC = np.copy(T)

        ''' Includding contour conditions '''
        T[wells['ws_p'],:] = 0
        T[wells['ws_p'], wells['ws_p']] = 1
        return T

    def pressure_independent_term(self, fprop):
        vector = ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP
        pressure_term = vector * fprop.P
        return pressure_term

    def capillary_and_gravity_independent_term(self, fprop, kprop):

        t0_j = self.t0_internal_faces_prod * ctes.pretransmissibility_internal_faces
        t0_k = ctes.g * np.sum(fprop.phase_densities_internal_faces * t0_j, axis=1)

        # Look for a better way to do this
        cap = np.zeros([ctes.n_volumes])
        grav = np.zeros([ctes.n_volumes,ctes.n_volumes])
        if any((ctes.z - ctes.z[0]) != 0):
            for i in range(kprop.n_components):
                lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                data = np.array([t0_k[i,:], t0_k[i,:], -t0_k[i,:], -t0_k[i,:]]).flatten()
                t0_rho = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes)).toarray()
                grav += t0_rho * self.dVtk[i,:]

                for j in range(kprop.n_phases):
                    lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                    cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                    data = np.array([t0_j[i,j,:], t0_j[i,j,:], -t0_j[i,j,:], -t0_j[i,j,:]]).flatten()
                    t0 = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes))*self.dVtk[i,:]
                    cap += t0 @ fprop.Pcap[j,:]

        gravity_term = grav @ ctes.z

        # capillary_term = np.sum(self.dVtk * np.sum (fprop.component_molar_fractions *
        #         fprop.phase_molar_densities * fprop.mobilities * fprop.Pcap, axis = 1), axis = 0)
        return cap, gravity_term

    def volume_discrepancy_independent_term(self, fprop):
        volume_discrepancy_term = fprop.Vp - fprop.Vt
        if np.max(volume_discrepancy_term/fprop.Vp) > 5e-4:
            raise ValueError('diminuir delta_t')
        return volume_discrepancy_term

    def well_term(self, kprop, fprop, wells):
        self.q = np.zeros([kprop.n_components, ctes.n_volumes]) #for now
        well_term = np.zeros(ctes.n_volumes)
        #self.phase_existance = np.zeros([1, kprop.n_phases, ctes.n_volumes])

        #if kprop.load_k:
        #    self.phase_existance[0,0,:] = np.sign(fprop.L)
        #    self.phase_existance[0,1,:] = np.sign(fprop.V)**2
        #if kprop.load_w:
        #    self.phase_existance[0,kprop.n_phases-1,:] = 1

        if len(wells['ws_q']) > 0:
            #self.injected_fluid_molar_density = data_loaded['Wells']['P1']['ksi_inj']
            #self.q[:,wells['ws_q']] = (wells['values_q'] * self.injected_fluid_molar_density).T
            volume_q = np.argwhere(wells['value_type'][0] == 'volumetric').ravel()
            molar_q = np.argwhere(wells['value_type'] == 'molar').ravel()
            self.q[:,wells['ws_q'][volume_q]] = wells['values_q'][volume_q].T / self.dVtk[:,wells['ws_q'][volume_q]]
            self.q[:,wells['ws_q'][molar_q]] = np.sum(wells['values_q'][molar_q].T,axis=0)
            well_term[wells['ws_q'][volume_q]] = np.sum(wells['values_q'][volume_q].T,axis=0)
            well_term[wells['ws_q'][molar_q]] = np.sum(self.dVtk[:,wells['ws_q'][molar_q]] * self.q[:,wells['ws_q'][molar_q]], axis = 0)

        return well_term

    def update_independent_terms(self, M, fprop, kprop, wells, delta_t):
        self.pressure_term = self.pressure_independent_term(fprop)
        self.capillary_term, self.gravity_term = self.capillary_and_gravity_independent_term(fprop, kprop)
        self.volume_term = self.volume_discrepancy_independent_term(fprop)
        well_term = self.well_term(kprop, fprop, wells)
        independent_terms = self.pressure_term  - self.volume_term  +  delta_t * well_term - delta_t * (self.capillary_term + self.gravity_term)
        if len(wells['ws_p'])>1:
            bhp_ind = np.argwhere(M.volumes.center[wells['ws_p']][:,2] == min(M.volumes.center[wells['ws_p']][:,2])).ravel()
        else: bhp_ind = wells['ws_p']
        independent_terms[wells['ws_p']] = wells['values_p'] + ctes.g * fprop.phase_densities[0,0,wells['ws_p']]*(ctes.z[wells['ws_p']] - ctes.z[bhp_ind])
        return independent_terms

    def update_pressure(self, T, D, fprop):
        self.P = linalg.solve(T,D)

    def update_total_flux_internal_faces(self, M, fprop, kprop):
        Pot_hid = self.P + fprop.Pcap
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]

        self.total_flux_internal_faces = - np.sum(fprop.mobilities_internal_faces * ctes.pretransmissibility_internal_faces
                                         * ((Pot_hidj_up - Pot_hidj) - ctes.g * fprop.phase_densities_internal_faces
                                         * (z_up - z)), axis = 1)

    def update_flux_wells(self, fprop, kprop, wells, delta_t):
        """Tentar melhorar isso daqui"""

        wp = wells['ws_p']
        well_term = np.zeros([kprop.n_components,1,len(wp)])
        if len(wp)>=1:
            well_term[0,0,0:len(wp)] = (self.T_noCC[wp,:] @ self.P - self.pressure_term[wp] + self.volume_term[wp] ) / delta_t \
                            + self.capillary_term[wp] + self.gravity_term[wp]
            n_components = kprop.n_components
            ref = 0 # componente de refererencia para pegar a mobilidade

            for i in range(len(wp)):
                mob_k = np.sum(fprop.mobilities[:,:,wp[i]] * fprop.component_molar_fractions[:,:,wp[i]], axis = 1).ravel()

                if mob_k[0] == 0:
                    ref = 1
                    n_components -= 1

                if n_components > 1:
                    C = np.diag(np.ones(n_components))
                    C[1:,0] = - mob_k[(ref+1):]/ mob_k[ref]
                    C[0,:] = self.dVtk[ref:,wp[i]].T
                    self.q[ref:,wp[i]] = linalg.solve(C,well_term[:,0,i])

                else: self.q[ref,wp[i]] = well_term[0,0,i] / self.dVtk[ref,wp[i]]

            #(fprop.component_molar_fractions * self.phase_existance * fprop.phase_molar_densities).sum(axis=1)[:,wells['ws_inj']]
