from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from packs.compositional.properties_calculation import PropertiesCalc
from ..utils import relative_permeability2, phase_viscosity, capillary_pressure
from .partial_derivatives_new import PartialDerivatives
from packs.compositional.stability_check import StabilityCheck
from ..solvers.solvers_scipy.solver_sp import SolverSp
from .. import directories as direc
import scipy.sparse as sp
from scipy import linalg
from .update_time import delta_time
import numpy as np


class CompositionalFVM:
    def __init__(self, M, n_volumes, wells):
        self.n_volumes = n_volumes
        self.all_wells = wells['all_wells'].astype(int)
        self.v0 = M.faces.bridge_adjacencies(M.faces.internal,2,3)
        self.internal_faces = M.faces.internal
        self.n_internal_faces = len(self.v0[:,0])
        self.g = data_loaded['compositional_data']['g']
        self.z = -M.data['centroid_volumes'][:,2]

    def runIMPEC(self, M, data_loaded, data_impress, wells, fprop, kprop, delta_t):

        self.dVt_derivatives(fprop, kprop)
        self.update_gravity_term(M, kprop, fprop)
        self.get_faces_properties_upwind(fprop, kprop)
        self.get_phase_densities_internal_faces(fprop)

        r = 1/2# enter the while loop
        while (r!=1.):
            T = self.update_transmissibility( M, wells, fprop, kprop, delta_t)
            D = self.update_independent_terms(M, fprop, kprop, wells, delta_t)
            self.update_pressure(T, D, data_impress, fprop)
            self.update_flux_internal_faces(M, data_impress, wells, fprop, kprop, delta_t)
            self.update_flux_volumes(fprop, kprop)
            # For the composition calculation the time step might be different because it treats
            #composition explicitly and this explicit models are conditionally stable - wich can
            #be based on the CFL parameter.
            delta_t_new = delta_time.update_CFL(delta_t, wells, fprop)
            r = delta_t_new/delta_t
            delta_t = delta_t_new

        self.update_composition(fprop, kprop, wells, delta_t)
        return delta_t

    def update_gravity_term(self, M, kprop, fprop):
        self.G = self.g * fprop.phase_densities * self.z

    def dVt_derivatives(self, fprop, kprop):
        self.dVtk = np.zeros([kprop.n_components, self.n_volumes])

        if kprop.load_k:

            if not kprop.compressible_k:
                dVtP = np.zeros(self.n_volumes)
                self.dVtk[0:kprop.Nc,:] = 1 / fprop.ksi_L
            else: self.dVtk[0:kprop.Nc,:], dVtP = PartialDerivatives(fprop).get_all_derivatives(kprop, fprop)

        else: dVtP = np.zeros(self.n_volumes)

        if kprop.load_w:

            self.dVtk[kprop.n_components-1,:] = 1 / fprop.ksi_W
            dVwP = -fprop.component_mole_numbers[kprop.Nc,:] * fprop.ksi_W0 * fprop.Cw / fprop.ksi_W**2

        else: dVwP = np.zeros(self.n_volumes)

        self.dVtP = dVtP + dVwP

    def get_faces_properties_upwind(self, fprop, kprop):
        ''' Using one-point upwind approximation '''
        Pot_hid = fprop.P + fprop.Pcap - self.G[0,:,:]
        Pot_hidj = Pot_hid[:,self.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,self.v0[:,1]]

        self.mobilities_internal_faces = np.zeros([1, kprop.n_phases, self.n_internal_faces])
        mobilities_vols = fprop.mobilities[:,:,self.v0[:,0]]
        mobilities_vols_up = fprop.mobilities[:,:,self.v0[:,1]]
        self.mobilities_internal_faces[0,Pot_hidj_up <= Pot_hidj] = mobilities_vols[0,Pot_hidj_up <= Pot_hidj]
        self.mobilities_internal_faces[0,Pot_hidj_up > Pot_hidj] = mobilities_vols_up[0,Pot_hidj_up > Pot_hidj]

        self.phase_molar_densities_internal_faces = np.zeros([1, kprop.n_phases, self.n_internal_faces])
        phase_molar_densities_vols = fprop.phase_molar_densities[:,:,self.v0[:,0]]
        phase_molar_densities_vols_up = fprop.phase_molar_densities[:,:,self.v0[:,1]]
        self.phase_molar_densities_internal_faces[0,Pot_hidj_up <= Pot_hidj] = phase_molar_densities_vols[0,Pot_hidj_up <= Pot_hidj]
        self.phase_molar_densities_internal_faces[0,Pot_hidj_up > Pot_hidj] = phase_molar_densities_vols_up[0,Pot_hidj_up > Pot_hidj]

        self.component_molar_fractions_internal_faces = np.zeros([kprop.n_components, kprop.n_phases, self.n_internal_faces])
        component_molar_fractions_vols = fprop.component_molar_fractions[:,:,self.v0[:,0]]
        component_molar_fractions_vols_up = fprop.component_molar_fractions[:,:,self.v0[:,1]]
        self.component_molar_fractions_internal_faces[:,Pot_hidj_up <= Pot_hidj] = component_molar_fractions_vols[:,Pot_hidj_up <= Pot_hidj]
        self.component_molar_fractions_internal_faces[:,Pot_hidj_up > Pot_hidj] = component_molar_fractions_vols_up[:,Pot_hidj_up > Pot_hidj]

        self.t0_internal_faces_prod = self.component_molar_fractions_internal_faces * self.phase_molar_densities_internal_faces \
                                * self.mobilities_internal_faces

        #self.pretransmissibility_internal_faces
        # return t0_internal_faces_prod

    def get_phase_densities_internal_faces(self, fprop):

        self.phase_densities_internal_faces = (fprop.Vp[self.v0[:,0]] * fprop.phase_densities[:,:,self.v0[:,0]] +
                                                fprop.Vp[self.v0[:,1]] * fprop.phase_densities[:,:,self.v0[:,1]]) /  \
                                                (fprop.Vp[self.v0[:,0]] + fprop.Vp[self.v0[:,1]])

    def update_transmissibility(self, M, wells, fprop, kprop, delta_t):
        pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
        self.pretransmissibility_internal_faces = pretransmissibility_faces[self.internal_faces]#[100]*np.ones(len(self.internal_faces))

        ''' Transmissibility '''
        t0 = (self.t0_internal_faces_prod).sum(axis = 1)
        t0 = t0 * self.pretransmissibility_internal_faces
        T = np.zeros([self.n_volumes, self.n_volumes])

        # Look for a way of doing this not using a loop
        for i in range(kprop.n_components):
            lines = np.array([self.v0[:, 0], self.v0[:, 1], self.v0[:, 0], self.v0[:, 1]]).flatten()
            cols = np.array([self.v0[:, 1], self.v0[:, 0], self.v0[:, 0], self.v0[:, 1]]).flatten()
            data = np.array([-t0[i,:], -t0[i,:], +t0[i,:], +t0[i,:]]).flatten()

            Ta = (sp.csc_matrix((data, (lines, cols)), shape = (self.n_volumes, self.n_volumes))).toarray()
            T += Ta*self.dVtk[i,:]
        #T = (T * self.dVtk.T[:,np.newaxis, :]).sum(axis = 2)

        T = T * delta_t
        ''' Transmissibility diagonal term '''
        diag = np.diag((fprop.Vbulk * fprop.porosity * fprop.cf - self.dVtP))
        T += diag
        self.T_noCC = np.copy(T)

        ''' Includding contour conditions '''
        T[wells['ws_p'],:] = 0
        T[wells['ws_p'], wells['ws_p']] = 1
        return T

    def pressure_independent_term(self, fprop):
        vector = fprop.Vbulk * fprop.porosity * fprop.cf - self.dVtP
        pressure_term = vector * fprop.P
        return pressure_term

    def capillary_and_gravity_independent_term(self, wells, fprop, kprop):

        t0_j = self.t0_internal_faces_prod * self.pretransmissibility_internal_faces
        t0_k = self.g*np.sum(self.phase_densities_internal_faces * t0_j, axis=1)

        # Look for a better way to do this
        cap = np.zeros([self.n_volumes])
        grav = np.zeros([self.n_volumes,self.n_volumes])
        if any((self.z-self.z[0]) != 0):
            for i in range(kprop.n_components):
                lines = np.array([self.v0[:, 0], self.v0[:, 1], self.v0[:, 0], self.v0[:, 1]]).flatten()
                cols = np.array([self.v0[:, 1], self.v0[:, 0], self.v0[:, 0], self.v0[:, 1]]).flatten()
                data = np.array([t0_k[i,:], t0_k[i,:], -t0_k[i,:], -t0_k[i,:]]).flatten()
                t0_rho = sp.csc_matrix((data, (lines, cols)), shape = (self.n_volumes, self.n_volumes)).toarray()
                grav += t0_rho*self.dVtk[i,:]

                for j in range(kprop.n_phases):
                    lines = np.array([self.v0[:, 0], self.v0[:, 1], self.v0[:, 0], self.v0[:, 1]]).flatten()
                    cols = np.array([self.v0[:, 1], self.v0[:, 0], self.v0[:, 0], self.v0[:, 1]]).flatten()
                    data = np.array([t0_j[i,j,:], t0_j[i,j,:], -t0_j[i,j,:], -t0_j[i,j,:]]).flatten()
                    t0 = sp.csc_matrix((data, (lines, cols)), shape = (self.n_volumes, self.n_volumes))*self.dVtk[i,:]
                    cap += t0 @ fprop.Pcap[j,:]

        gravity_term = grav @ self.z

        # capillary_term = np.sum(self.dVtk * np.sum (fprop.component_molar_fractions *
        #         fprop.phase_molar_densities * fprop.mobilities * fprop.Pcap, axis = 1), axis = 0)
        return cap, gravity_term

    def volume_discrepancy_independent_term(self, fprop):
        volume_discrepancy_term = fprop.Vp - fprop.Vt
        if np.max(volume_discrepancy_term/fprop.Vp) > 5e-4:
            raise ValueError('diminuir delta_t')
        return volume_discrepancy_term

    def well_term(self, kprop, fprop, wells):
        self.q = np.zeros([kprop.n_components, self.n_volumes]) #for now
        #self.phase_existance = np.zeros([1, kprop.n_phases, self.n_volumes])

        #if kprop.load_k:
        #    self.phase_existance[0,0,:] = np.sign(fprop.L)
        #    self.phase_existance[0,1,:] = np.sign(fprop.V)**2
        #if kprop.load_w:
        #    self.phase_existance[0,kprop.n_phases-1,:] = 1

        if len(wells['ws_q']) > 0:
            self.injected_fluid_molar_density = data_loaded['Wells']['P1']['ksi_inj']
            self.q[:,wells['ws_q']] = (wells['values_q'] * self.injected_fluid_molar_density).T #[:,np.newaxis]
            #(fprop.component_molar_fractions * self.phase_existance * fprop.phase_molar_densities).sum(axis=1)[:,wells['ws_inj']]
             #for now - its going to change
        well_term = np.sum(self.dVtk * self.q, axis = 0)
        return well_term

    def update_independent_terms(self, M, fprop, kprop, wells, delta_t):
        self.pressure_term = self.pressure_independent_term(fprop)
        self.capillary_term, self.gravity_term = self.capillary_and_gravity_independent_term(wells, fprop, kprop)
        self.volume_term = self.volume_discrepancy_independent_term(fprop)
        well_term = self.well_term(kprop, fprop, wells)
        independent_terms = self.pressure_term  - self.volume_term  +  delta_t * well_term - delta_t * (self.capillary_term + self.gravity_term)
        if len(wells['ws_p'])>1:
            bhp_ind = np.argwhere(M.volumes.center[wells['ws_p']][:,2] == min(M.volumes.center[wells['ws_p']][:,2])).ravel()
        else: bhp_ind = wells['ws_p']
        independent_terms[wells['ws_p']] = wells['values_p'] + self.g * fprop.phase_densities[0,0,wells['ws_p']]*(self.z[wells['ws_p']] - self.z[bhp_ind])
        return independent_terms

    def update_pressure(self, T, D, data_impress, fprop):
        fprop.P = linalg.solve(T,D)
        #fprop.P = np.linalg.solve(T,D)
        data_impress['pressure'] = fprop.P

    def update_flux_internal_faces(self, M, data_impress, wells, fprop, kprop, delta_t):
        Pot_hid = fprop.P + fprop.Pcap
        Pot_hidj = Pot_hid[:,self.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,self.v0[:,1]]
        z = self.z[self.v0[:,0]]
        z_up = self.z[self.v0[:,1]]

        #total_flux_internal_faces = - np.sum(self.mobilities_internal_faces * self.pretransmissibility_internal_faces
        #                                 * ((Pot_hidj_up - Pot_hidj) - self.g * self.phase_densities_internal_faces
        #                                 * (z_up - z)), axis = 1)

        self.phase_flux_internal_faces = - (self.mobilities_internal_faces * self.pretransmissibility_internal_faces
                                     * (Pot_hidj_up - Pot_hidj - self.g * self.phase_densities_internal_faces
                                     * (z_up - z)))

        #total_flux_internal_faces1 = - np.sum(self.mobilities_internal_faces * self.pretransmissibility_internal_faces
        #                                 * (Pot_hidj_up - Pot_hidj), axis=1)
        #total_flux_internal_faces2 = total_flux_internal_faces1 + np.sum(self.mobilities_internal_faces *
        #                                self.pretransmissibility_internal_faces * self.phase_densities_internal_faces
        #                                    * self.g *(z_up - z), axis = 1)
        #frj = self.mobilities_internal_faces[0,:,:] / np.sum(self.mobilities_internal_faces[0,:,:], axis = 0)
        #phase_flux_internal_faces = frj * total_flux_internal_faces

        # M.flux_faces[M.faces.internal] = total_flux_internal_faces * M.faces.normal[M.faces.internal].T

    def update_flux_volumes(self, fprop, kprop):

        component_flux_internal_faces = np.sum(self.component_molar_fractions_internal_faces * self.phase_molar_densities_internal_faces *
                                self.phase_flux_internal_faces, axis = 1)
        cx = np.arange(kprop.n_components)
        lines = np.array([np.repeat(cx,len(self.v0[:,0])), np.repeat(cx,len(self.v0[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(self.v0[:,0],kprop.n_components), np.tile(self.v0[:,1], kprop.n_components)]).flatten()
        data = np.array([-component_flux_internal_faces, component_flux_internal_faces]).flatten()
        fprop.component_flux_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (kprop.n_components, self.n_volumes)).toarray()

    def update_flux_wells(self, fprop, kprop, wells, delta_t):
        """Tentar melhorar isso daqui"""
        wsp = set.difference(set(wells['ws_prod']), set(wells['ws_q']))
        wp = np.array(list(wsp))
        well_term = np.zeros([kprop.n_components,1,len(wp)])
        if len(wp)>=1:
            well_term[0,0,0:len(wp)] = (self.T_noCC[wp,:] @ fprop.P - self.pressure_term[wp] + self.volume_term[wp] ) / delta_t \
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

        if len(wells['ws_q']) > 0:
            self.q[:,wells['ws_q']] = (wells['values_q'] * self.injected_fluid_molar_density).T

            #(fprop.component_molar_fractions * self.phase_existance * fprop.phase_molar_densities).sum(axis=1)[:,wells['ws_inj']]

    def update_composition(self, fprop, kprop, wells, delta_t):
        self.update_flux_wells(fprop, kprop, wells, delta_t)
        import pdb; pdb.set_trace()
        fprop.component_mole_numbers = fprop.component_mole_numbers + delta_t * (self.q + fprop.component_flux_vols_total)
        fprop.z = fprop.component_mole_numbers[0:kprop.Nc,:] / np.sum(fprop.component_mole_numbers[0:kprop.Nc,:], axis = 0)
