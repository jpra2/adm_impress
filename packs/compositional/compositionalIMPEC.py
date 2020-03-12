from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from packs.compositional.properties_calculation import PropertiesCalc
from ..utils import relative_permeability2, phase_viscosity, capillary_pressure
from .partial_derivatives import PartialDerivatives
from .. import directories as direc
import scipy.sparse as sp
import numpy as np


class CompositionalFVM:
    def __init__(self, M, data_impress, wells, fprop, fprop_block, kprop, load, deltaT):
        self.n_phases = 3 #includding water
        self.n_volumes = data_impress.len_entities['volumes']
        self.n_components = fprop.Nc + 1
        self.all_wells = wells['all_wells'].astype(int)
        fprop.Vbulk = data_impress['volume']
        self.porosity = data_impress['poro']
        self.cf = np.array(data_loaded['compositional_data']['rock_compressibility']).astype(float)
        self.relative_permeability = getattr(relative_permeability2, data_loaded['compositional_data']['relative_permeability'])
        self.relative_permeability = self.relative_permeability()
        self.phase_viscosity = getattr(phase_viscosity, data_loaded['compositional_data']['phase_viscosity'])
        self.phase_viscosity = self.phase_viscosity(self.n_volumes, fprop, kprop)
        self.v0 = M.faces.bridge_adjacencies(M.faces.internal,2,3)
        self.runIMPEC(M, data_loaded, data_impress, wells, fprop, fprop_block, kprop, deltaT)

    def runIMPEC(self, M, data_loaded, data_impress, wells, fprop, fprop_block, kprop, deltaT):
        self.update_relative_permeabilities(fprop)
        self.update_phase_viscosities(data_loaded, fprop, kprop)
        self.update_mobilities()
        self.dVt_derivatives(fprop, fprop_block, kprop)
        Pcap = self.update_capillary_pressure(data_loaded)
        self.get_faces_properties_upwind(M, fprop)
        T = self.update_transmissibility(M, data_impress, wells, data_loaded, fprop)
        D = self.update_independent_terms(fprop, data_loaded, wells, deltaT)
        self.update_pressure(T, D, data_impress, fprop)
        self.update_flux_internal_faces(M, fprop)
        self.update_composition(M, fprop, deltaT)

    def update_relative_permeabilities(self, fprop):
        # So, Sw, Sg = self.update_saturations_without_contours()
        # saturations = np.array([So, Sg, Sw])
        saturations = np.array([fprop.So, fprop.Sg, fprop.Sw])
        kro,krg,krw = self.relative_permeability(saturations)
        self.relative_permeabilities = np.zeros([1, self.n_phases, self.n_volumes])
        self.relative_permeabilities[0,0,:] = kro
        self.relative_permeabilities[0,1,:] = krg
        self.relative_permeabilities[0,2,:] = krw

    def update_phase_viscosities(self, data_loaded, fprop, kprop):
        mi_W = data_loaded['compositional_data']['water_data']['mi_W']
        self.phase_viscosities = np.zeros(self.relative_permeabilities.shape)
        self.phase_viscosities_oil_and_gas = self.phase_viscosity(fprop, kprop)
        self.phase_viscosities[0,0,:] = self.phase_viscosities_oil_and_gas[0,0,:]
        self.phase_viscosities[0,1,:] = self.phase_viscosities_oil_and_gas[0,1,:]
        self.phase_viscosities[0,2,:] = mi_W

    def update_mobilities(self):
        self.mobilities = self.relative_permeabilities / self.phase_viscosities

    def dVt_derivatives(self, fprop, fprop_block, kprop):
        self.dVtk = np.zeros([fprop.Nc + 1, self.n_volumes])
        self.dVtk[0:fprop.Nc,:], dVtP = PartialDerivatives().dVt_derivatives(fprop, fprop_block, kprop)
        self.dVtk[fprop.Nc,:] = 1 / fprop.ksi_W
        dVwP = np.zeros(self.n_volumes)
        self.dVtP = dVtP + dVwP

    def get_mobilities_upwind(self, M, fprop):
        ''' Using one-point upwind approximation '''
        P_Pcap = fprop.P + self.Pcap
        Pj = P_Pcap[:,self.v0[:,0]]
        Pj_up = P_Pcap[:,self.v0[:,1]]

        self.mobilities_internal_faces = np.zeros([1, self.n_phases, len(M.faces.internal)])
        mobilities_vols = self.mobilities[:,:,self.v0[:,0]]
        mobilities_vols_up = self.mobilities[:,:,self.v0[:,1]]
        self.mobilities_internal_faces[0,Pj_up <= Pj] = mobilities_vols[0,Pj_up <= Pj]
        self.mobilities_internal_faces[0,Pj_up > Pj] = mobilities_vols_up[0,Pj_up > Pj]

    def get_faces_properties_upwind(self, M, fprop):
        self.get_mobilities_upwind(M, fprop)
        P_Pcap = fprop.P + self.Pcap
        Pj = P_Pcap[:,self.v0[:,0]]
        Pj_up = P_Pcap[:,self.v0[:,1]]

        self.phase_molar_densities_internal_faces = np.zeros([1, self.n_phases, len(M.faces.internal)])
        phase_molar_densities_vols = fprop.phase_molar_densities[:,:,self.v0[:,0]]
        phase_molar_densities_vols_up = fprop.phase_molar_densities[:,:,self.v0[:,1]]
        self.phase_molar_densities_internal_faces[0,Pj_up <= Pj] = phase_molar_densities_vols[0,Pj_up <= Pj]
        self.phase_molar_densities_internal_faces[0,Pj_up > Pj] = phase_molar_densities_vols_up[0,Pj_up > Pj]

        self.component_molar_fractions_internal_faces = np.zeros([self.n_components, self.n_phases, len(M.faces.internal)])
        component_molar_fractions_vols = fprop.component_molar_fractions[:,:,self.v0[:,0]]
        component_molar_fractions_vols_up = fprop.component_molar_fractions[:,:,self.v0[:,1]]
        self.component_molar_fractions_internal_faces[:,Pj_up <= Pj] = component_molar_fractions_vols[:,Pj_up <= Pj]
        self.component_molar_fractions_internal_faces[:,Pj_up > Pj] = component_molar_fractions_vols_up[:,Pj_up > Pj]

        self.t0_internal_faces_prod = self.component_molar_fractions_internal_faces * self.phase_molar_densities_internal_faces \
                                * self.mobilities_internal_faces

        # return t0_internal_faces_prod

    def update_transmissibility(self, M, data_impress, wells, data_loaded, fprop):
        pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
        self.pretransmissibility_internal_faces = pretransmissibility_faces[M.faces.internal]

        ''' Transmissibility '''
        t0 = (self.t0_internal_faces_prod).sum(axis = 1)
        t0 = t0 * self.pretransmissibility_internal_faces
        T = np.zeros([self.n_volumes, self.n_volumes, self.n_components])

        # Look for a way of doing this not using a loop
        for i in range(self.n_components):
            lines = np.array([self.v0[:, 0], self.v0[:, 1], self.v0[:, 0], self.v0[:, 1]]).flatten()
            cols = np.array([self.v0[:, 1], self.v0[:, 0], self.v0[:, 0], self.v0[:, 1]]).flatten()
            data = np.array([t0[i,:], t0[i,:], -t0[i,:], -t0[i,:]]).flatten()

            Ta = sp.csc_matrix((data, (lines, cols)), shape = (self.n_volumes, self.n_volumes))
            T[:,:,i] = Ta.toarray()

        T = (T * self.dVtk.T[:,np.newaxis, :]).sum(axis = 2)

        # Teste
        t0 = np.sum(self.dVtk * t0, axis = 0)
        lines = np.array([self.v0[:, 0], self.v0[:, 1], self.v0[:, 0], self.v0[:, 1]]).flatten()
        cols = np.array([self.v0[:, 1], self.v0[:, 0], self.v0[:, 0], self.v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        Ttest = sp.csc_matrix((data, (lines, cols)), shape = (self.n_volumes, self.n_volumes))
        import pdb; pdb.set_trace()

        ''' Transmissibility diagonal term '''
        diag = np.diag(fprop.Vbulk * self.porosity * self.cf - self.dVtP)
        T += diag

        ''' Includding contour conditions '''
        T[wells['ws_p'],:] = 0
        T[wells['ws_p'], wells['ws_p']] = 1
        return T


    def update_capillary_pressure(self, data_loaded):
        """ not working yet"""
        # get_capillary_pressure = getattr(capillary_pressure, data_loaded['compositional_data']['capillary_pressure'])
        # get_capillary_pressure = get_capillary_pressure(data_loaded, self.phase_molar_densities, self.component_molar_fractions,
        #                         self.pretransmissibility_internal_faces)
        # So, Sw, Sg = self.update_saturations_without_contours()
        # Pcow, Pcog = get_capillary_pressure(Sw, So, Sg)

        self.Pcap = np.zeros([self.n_phases,self.n_volumes])
        # Pcap[0,0,:] = Pcog
        # Pcap[0,1,:] = Pcow

    def pressure_independent_term(self, fprop):
        vector = fprop.Vbulk * self.porosity * self.cf - self.dVtP
        pressure_term = vector * fprop.P
        return pressure_term

    def capillary_independent_term(self, fprop, data_loaded):
        t0 = self.t0_internal_faces_prod

        # Look for a better way to do this
        cap = np.zeros([self.n_components, self.n_phases, self.n_volumes])
        for i in range(self.n_components):
            for j in range(self.n_phases):
                lines = np.array([self.v0[:, 0], self.v0[:, 1], self.v0[:, 0], self.v0[:, 1]]).flatten()
                cols = np.array([self.v0[:, 1], self.v0[:, 0], self.v0[:, 0], self.v0[:, 1]]).flatten()
                data = np.array([t0[i,j,:], t0[i,j,:], -t0[i,j,:], -t0[i,j,:]]).flatten()
                T = sp.csc_matrix((data, (lines, cols)), shape = (self.n_volumes, self.n_volumes))
                cap[i,j,:] = T.toarray() @ self.Pcap[j,:]

        cap = cap.sum(axis = 1)
        capillary_term = (cap * self.dVtk).sum(axis = 0)

        # capillary_term = np.sum(self.dVtk * np.sum (fprop.component_molar_fractions *
        #         fprop.phase_molar_densities * self.mobilities * self.Pcap, axis = 1), axis = 0)
        return capillary_term

    def volume_discrepancy_independent_term(self, fprop):
        volume_discrepancy_term = fprop.Vp - fprop.Vt
        return volume_discrepancy_term

    def well_term(self, wells):
        self.q = np.zeros([self.n_components, self.n_volumes]) #for now
        self.q[:,wells['ws_q']] = wells['values_q'] #for now - its going to change
        well_term = np.sum(self.dVtk * self.q, axis = 0)
        return well_term

    def update_independent_terms(self, fprop, data_loaded, wells, deltaT):
        pressure_term = self.pressure_independent_term(fprop)
        capillary_term = self.capillary_independent_term(fprop, data_loaded)
        volume_term = self.volume_discrepancy_independent_term(fprop)
        well_term = self.well_term(wells)
        independent_terms = pressure_term - volume_term + deltaT * well_term #+ deltaT * capillary_term
        independent_terms[wells['ws_p']] = wells['values_p']
        return independent_terms

    def update_pressure(self, T, D, data_impress, fprop):
        fprop.P = np.array(np.linalg.inv(T)@(D[:,np.newaxis])).ravel()
        data_impress['pressure'] = fprop.P

    def update_flux_internal_faces(self, M, fprop):
        P_Pcap = fprop.P + self.Pcap
        Pj = P_Pcap[:,self.v0[:,0]]
        Pj_up = P_Pcap[:,self.v0[:,1]]
        total_flux_internal_faces = - np.sum(self.mobilities_internal_faces * self.pretransmissibility_internal_faces * (Pj_up - Pj) ,axis = 1)
        self.get_mobilities_upwind(M, fprop) # a mobilidade aqui é em n + 1
        phase_flux_internal_faces = self.mobilities_internal_faces / np.sum(self.mobilities_internal_faces, axis = 1) * total_flux_internal_faces
        # M.flux_faces[M.faces.internal] = total_flux_internal_faces * M.faces.normal[M.faces.internal].T
        return phase_flux_internal_faces

    def update_composition(self, M, fprop, deltaT):
        phase_flux_internal_faces = self.update_flux_internal_faces(M, fprop)
        component_flux = np.sum(self.component_molar_fractions_internal_faces * self.phase_molar_densities_internal_faces * phase_flux_internal_faces, axis = 1)

        cx = np.arange(self.n_components)
        lines = np.array([np.repeat(cx,len(self.v0[:,0])), np.repeat(cx,len(self.v0[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(self.v0[:,0],self.n_components), np.tile(self.v0[:,1], self.n_components)]).flatten()
        data = np.array([component_flux, -component_flux]).flatten()
        fprop.component_flux_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (self.n_components, self.n_volumes)).toarray()

        # só ta funcionando pra 1d:
        #flux_vols = np.zeros([self.n_components,2,self.n_volumes])
        #flux_vols[:,0,self.v0[:,0]] = flux
        #flux_vols[:,1,self.v0[:,1]] = -flux
        #flux_vols_total = np.sum(flux_vols,axis = 1)

        fprop.component_mole_numbers = fprop.component_mole_numbers + deltaT * (self.q + fprop.component_flux_vols_total)
        fprop.z = fprop.component_mole_numbers[0:fprop.Nc,:] / np.sum(fprop.component_mole_numbers[0:fprop.Nc,:], axis = 0)
