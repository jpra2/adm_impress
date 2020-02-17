from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from packs.compositional.properties_calculation import PropertiesCalc
from ..utils import relative_permeability2, phase_viscosity, capillary_pressure
from .partial_derivatives import PartialDerivatives
from .. import directories as direc
import scipy.sparse as sp
import numpy as np


class IMPEC(DataManager):
    def __init__(self, M, data_impress, wells, fluid_properties, elements_lv0, load, data_name: str='CompositionalTPFA.npz'):
        super().__init__(data_name, load=load)
        self.n_phases = 3 #includding water
        self.n_volumes = data_impress.len_entities['volumes']
        self.n_components = fluid_properties.Nc + 1
        self.internal_faces = elements_lv0['internal_faces']
        self.n_faces = len(self.internal_faces)
        self.all_wells = wells['all_wells'].astype(int)
        self.Vbulk = data_impress['volume']
        self.porosity = data_impress['poro']
        self.cf = np.array(data_loaded['compositional_data']['rock_compressibility']).astype(float)
        self.relative_permeability = getattr(relative_permeability2, data_loaded['compositional_data']['relative_permeability'])
        self.relative_permeability = self.relative_permeability()
        self.phase_viscosity = getattr(phase_viscosity, data_loaded['compositional_data']['phase_viscosity'])
        self.phase_viscosity = self.phase_viscosity(self.n_volumes, fluid_properties)

        if not load:
            self.loop = 0
            self.vpi = 0.0
            self.t = 0.0
            self.contador_vtk = 0
            self.run(M, data_loaded, data_impress, wells, fluid_properties, elements_lv0)
        else: self.load_infos()

    def run(self, M, data_loaded, data_impress, wells, fluid_properties, elements_lv0):
        self.update_relative_permeabilities(fluid_properties)
        self.update_phase_viscosities(data_loaded, fluid_properties)
        self.update_mobilities()
        self.dVt_derivatives(data_impress, fluid_properties)
        Pcap = self.update_capillary_pressure(data_loaded)
        T = self.update_transmissibility(M, data_impress, wells, data_loaded, elements_lv0, fluid_properties)
        D = self.update_independent_terms(fluid_properties, data_loaded, wells)
        IMPEC.update_pressure(T, D, data_impress, fluid_properties)
        self.update_composition(fluid_properties, elements_lv0)

    def solve_composition(self):
        # self.component_mole_numbers = self.component_molar_fractions * self.phase_mole_numbers # em t = 0
        fluid_properties.component_mole_numbers = fluid_properties.component_mole_numbers + deltaT * (self.q +
        np.sum(fluid_properties.component_molar_fractions * fluid_properties.phase_molar_densities * phase_flux, axis = 1))

    def update_relative_permeabilities(self, fluid_properties):
        # So, Sw, Sg = self.update_saturations_without_contours()
        # saturations = np.array([So, Sg, Sw])
        saturations = np.array([fluid_properties.So, fluid_properties.Sg, fluid_properties.Sw])
        kro,krg,krw = self.relative_permeability(saturations)
        self.relative_permeabilities = np.zeros([1, self.n_phases, self.n_volumes])
        self.relative_permeabilities[0,0,:] = kro
        self.relative_permeabilities[0,1,:] = krg
        self.relative_permeabilities[0,2,:] = krw

    def update_phase_viscosities(self, data_loaded, fluid_properties):
        mi_W = data_loaded['compositional_data']['water_data']['mi_W']
        self.phase_viscosities = np.zeros(self.relative_permeabilities.shape)
        self.phase_viscosities_oil_and_gas = self.phase_viscosity(fluid_properties)
        self.phase_viscosities[0,0,:] = self.phase_viscosities_oil_and_gas[0,0,:]
        self.phase_viscosities[0,1,:] = self.phase_viscosities_oil_and_gas[0,1,:]
        self.phase_viscosities[0,2,:] = mi_W

    def update_mobilities(self):
        self.mobilities = self.relative_permeabilities / self.phase_viscosities

    def dVt_derivatives(self, data_impress, fluid_properties):
        self.dVtk = np.zeros([fluid_properties.Nc + 1, self.n_volumes])
        self.dVtk[0:fluid_properties.Nc,:], dVtP = PartialDerivatives().dVt_derivatives(fluid_properties)
        self.dVtk[fluid_properties.Nc,:] = 1 / fluid_properties.ksi_W
        dVwP = np.zeros(self.n_volumes)
        self.dVtP = dVtP + dVwP

    def get_faces_properties_upwind(self, data_loaded, fluid_properties, v0):
        P_Pcap = fluid_properties.P + self.Pcap
        Pj = P_Pcap[:,:,v0[:,0]]
        Pj_up = P_Pcap[:,:,v0[:,1]]

        self.mobilities_internal_faces = np.zeros([1,self.n_phases,self.n_faces])
        mobilities_vols = self.mobilities[:,:,v0[:,0]]
        mobilities_vols_up = self.mobilities[:,:,v0[:,1]]
        self.mobilities_internal_faces[0,Pj_up[0,:,:] <= Pj[0,:,:]] = mobilities_vols[0,Pj_up[0,:,:] <= Pj[0,:,:]]
        self.mobilities_internal_faces[0,Pj_up[0,:,:] > Pj[0,:,:]] = mobilities_vols_up[0,Pj_up[0,:,:] > Pj[0,:,:]]

        self.phase_molar_densities_internal_faces = np.zeros([1,self.n_phases,self.n_faces])
        phase_molar_densities_vols = fluid_properties.phase_molar_densities[:,:,v0[:,0]]
        phase_molar_densities_vols_up = fluid_properties.phase_molar_densities[:,:,v0[:,1]]
        self.phase_molar_densities_internal_faces[0,Pj_up[0,:,:] <= Pj[0,:,:]] = phase_molar_densities_vols[0,Pj_up[0,:,:] <= Pj[0,:,:]]
        self.phase_molar_densities_internal_faces[0,Pj_up[0,:,:] > Pj[0,:,:]] = phase_molar_densities_vols_up[0,Pj_up[0,:,:] > Pj[0,:,:]]

        self.component_molar_fractions_internal_faces = np.zeros([self.n_components,self.n_phases,self.n_faces])
        component_molar_fractions_vols = fluid_properties.component_molar_fractions[:,:,v0[:,0]]
        component_molar_fractions_vols_up = fluid_properties.component_molar_fractions[:,:,v0[:,1]]
        self.component_molar_fractions_internal_faces[Pj_up <= Pj] = component_molar_fractions_vols[Pj_up <= Pj]
        self.component_molar_fractions_internal_faces[Pj_up > Pj] = component_molar_fractions_vols_up[Pj_up > Pj]

        t0_internal_faces_prod = self.component_molar_fractions_faces * self.phase_molar_densities_faces * self.mobilities_faces

        v = M.faces.bridge_adjacencies(self.all_wells,2,3).flatten()
        self.mobilities_faces = np.copy(self.mobilities)
        self.mobilities_faces[:,:,self.all_wells] = self.mobilities[self.all_wells]
        self.mobilities_faces[:,:,] = self.mobilities_internal_faces
        # t0_vols_prod_all = fluid_properties.component_molar_fractions * fluid_properties.phase_molar_densities * self.mobilities
        # t0_vols_prod = t0_vols_prod_all[:,:,v0[:,0]]
        # t0_vols_prod_up = t0_vols_prod_all[:,:,v0[:,1]]
        # t0_face_prod = np.zeros([self.n_components, self.n_phases, self.n_faces])
        # t0_face_prod[Pj_up <= Pj] = t0_vols_prod[Pj_up <= Pj]
        # t0_face_prod[Pj_up > Pj] = t0_vols_prod_up[Pj_up > Pj]

        return t0_internal_faces_prod

    def update_transmissibility(self, M, data_impress, wells, data_loaded, elements_lv0, fluid_properties):
        v0 = elements_lv0['neig_internal_faces']
        pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
        self.pretransmissibility_internal_faces = pretransmissibility_faces[self.internal_faces]

        ''' Using one-point upwind approximation '''
        dVtk = np.delete(self.dVtk, self.all_wells, axis = 1)
        t0_face_prod = self.get_faces_properties_upwind(data_loaded, fluid_properties, v0)

        ''' Transmissibility '''
        t0 = (dVtk * (t0_face_prod).sum(axis=1)).sum(axis=0)
        t0 = t0 * self.pretransmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape = (self.n_volumes, self.n_volumes))

        ''' Transmissibility diagonal term '''
        diag = np.diag(self.Vbulk * self.porosity * self.cf - self.dVtP)
        T += diag
        self['Tini'] = T

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

        self.Pcap = np.zeros([self.n_components,self.n_phases,self.n_volumes])
        # Pcap[0,0,:] = Pcog
        # Pcap[0,1,:] = Pcow

    def update_deltaT(self):
        """ need to do that """
        deltaT = 2
        return deltaT

    def pressure_independent_term(self, fluid_properties):
        vector = self.Vbulk * self.porosity * self.cf - self.dVtP
        pressure_term = vector * fluid_properties.P
        return pressure_term

    def capillary_independent_term(self, fluid_properties, data_loaded):
        capillary_term = np.sum(self.dVtk * np.sum (fluid_properties.component_molar_fractions *
                fluid_properties.phase_molar_densities * self.mobilities * self.Pcap, axis = 1), axis = 0)
        return capillary_term

    def volume_discrepancy_independent_term(self, fluid_properties):
        volume_discrepancy_term = fluid_properties.Vp - np.sum(np.sum(fluid_properties.phase_mole_numbers
                                / fluid_properties.phase_molar_densities, axis = 0), axis = 0)
        return volume_discrepancy_term

    def well_term(self, wells):
        self.q = np.zeros([self.n_components,self.n_volumes]) #for now
        self.q[:,wells['ws_q']] = wells['values_q'] #for now - its going to change
        well_term = np.sum(self.dVtk * self.q, axis = 0)
        return well_term

    def update_independent_terms(self, fluid_properties, data_loaded, wells):
        deltaT = self.update_deltaT()
        pressure_term = self.pressure_independent_term(fluid_properties)
        capillary_term = self.capillary_independent_term(fluid_properties, data_loaded)
        volume_term = self.volume_discrepancy_independent_term(fluid_properties)
        well_term = self.well_term(wells)
        independent_terms = pressure_term - volume_term + deltaT * well_term #+ deltaT * capillary_term
        independent_terms[wells['ws_p']] = wells['values_p']
        return independent_terms

    def update_pressure(T, D, data_impress, fluid_properties):
        fluid_properties.P = np.array(np.linalg.inv(T)*(D[:,np.newaxis])).ravel()
        data_impress['pressure'] = fluid_properties.P

    def update_composition(self, fluid_properties, elements_lv0):
        v0 = elements_lv0['neig_internal_faces']
        deltaT = self.update_deltaT()
        P_Pcap = fluid_properties.P + self.Pcap
        Pj = P_Pcap[:,:,v0[:,0]]
        Pj_up = P_Pcap[:,:,v0[:,1]]

        total_flux = - np.sum(self.mobilities_internal_faces * self.pretransmissibility_internal_faces ,axis = 1) * (Pj_up[0,:,:] - Pj[0,:,:])
        phase_flux = self.mobilities_internal_faces / np.sum(self.mobilities_internal_faces, axis = 1) * total_flux

        fluid_properties.component_mole_numbers = fluid_properties.component_mole_numbers + deltaT * (self.q +
        np.sum(self.component_molar_fractions_faces * self.phase_molar_densities_faces * phase_flux, axis = 1))
