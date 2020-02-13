from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from ..utils import relative_permeability2, phase_viscosity, capillary_pressure
from .. import directories as direc
from .partial_derivatives import PartialDerivatives
import scipy.sparse as sp
import numpy as np


class CompositionalIMPEC(DataManager):
    def __init__(self, M, data_impress, wells, fluid_properties, elements_lv0, load, data_name: str='CompositionalTPFA.npz'):
        super().__init__(data_name, load=load)
        self.n_phases = 3 #includding water
        self.n_volumes = data_impress.len_entities['volumes']
        self.n_components = fluid_properties.Nc + 1 #1 = water
        self.internal_faces = elements_lv0['internal_faces']
        self.relative_permeability = getattr(relative_permeability2, data_loaded['compositional_data']['relative_permeability'])
        self.relative_permeability = self.relative_permeability()
        self.phase_viscosity = getattr(phase_viscosity, data_loaded['compositional_data']['phase_viscosity'])
        self.phase_viscosity = self.phase_viscosity(len(self.internal_faces), fluid_properties)
        self.all_wells = wells['all_wells'].astype(int)
        self.Vbulk = data_impress['volume']
        self.porosity = data_impress['poro']
        self.cf = np.array(data_loaded['compositional_data']['rock_compressibility']).astype(float)

        if not load:
            self.loop = 0
            self.vpi = 0.0
            self.t = 0.0
            self.contador_vtk = 0
            self.run(M, data_loaded, data_impress, wells, fluid_properties, elements_lv0)
        else: self.load_infos()

    def run(self, M, data_loaded, data_impress, wells, fluid_properties, elements_lv0):
        self.get_water_data(data_loaded['compositional_data']['water_data'],fluid_properties)
        self.set_properties(fluid_properties, elements_lv0)
        self.update_saturations(data_impress, wells, fluid_properties)
        self.update_relative_permeabilities(fluid_properties)
        self.update_phase_viscosities(fluid_properties)
        T = self.update_transmissibility(M, data_impress, data_loaded, elements_lv0, fluid_properties)
        D = self.update_independent_terms(fluid_properties, data_loaded)

        # self.mesh_name = os.path.join(direc.flying, 'compositional_')
        # self.all_compositional_results = self.get_empty_current_compositional_results()
        # self.solver = SolverSp()

    def get_water_data(self, water_data, fluid_properties):
        """ Attention: Water is assumed to be incompressible """
        fluid_properties.rho_W = water_data['rho_W']
        self.mi_W = water_data['mi_W']
        fluid_properties.Mw_w = water_data['Mw_w']
        fluid_properties.eta_W = fluid_properties.rho_W/fluid_properties.Mw_w

    def set_properties(self, fluid_properties, elements_lv0):
        self.component_molar_fractions = np.zeros([self.n_components, self.n_phases, self.n_volumes])
        self.phase_molar_densities = np.zeros([1, self.n_phases, self.n_volumes])

        self.component_molar_fractions[0:fluid_properties.Nc,0,:] = fluid_properties.x
        self.component_molar_fractions[0:fluid_properties.Nc,1,:] = fluid_properties.y
        self.component_molar_fractions[fluid_properties.Nc,2,:] = 1 #water molar fraction in water component

        self.phase_molar_densities[0,0,:] = fluid_properties.eta_L
        self.phase_molar_densities[0,1,:] = fluid_properties.eta_V
        self.phase_molar_densities[0,2,:] = fluid_properties.eta_W

    def update_saturations(self, data_impress, wells, fluid_properties):
        self.Sw = data_impress['saturation']

        if fluid_properties.V != 0:
            self.Sg = (1 - self.Sw) * (fluid_properties.V / fluid_properties.rho_V) / \
                (fluid_properties.V / fluid_properties.rho_V +
                fluid_properties.L / fluid_properties.rho_L )
        else: self.Sg = np.zeros(self.n_volumes)
        self.So = 1 - self.Sw - self.Sg

    def update_saturations_without_contours(self):
        So = np.delete(self.So, self.all_wells)
        Sw = np.delete(self.Sw, self.all_wells)
        Sg = np.delete(self.Sg, self.all_wells)
        return So, Sw, Sg

    def update_phase_volumes(self, fluid_properties, data_impress):
        Pf = np.array(data_loaded['compositional_data']['Pf']).astype(float)
        # Pf - reference pressure at witch porosity is obtained
        self.Vp = self.porosity * self.Vbulk * (1 + self.cf*(fluid_properties.P - Pf))
        self.Vo = self.Vp * self.So
        self.Vg = self.Vp * self.Sg
        self.Vw = self.Vp * self.Sw

    def update_phase_mole_numbers(self, fluid_properties, data_impress):
        self.update_phase_volumes(fluid_properties, data_impress)
        self.phase_mole_numbers = np.zeros([1, self.n_phases, self.n_volumes])
        eta = np.ones([1, 2, self.n_volumes])
        V = np.ones([1, 2, self.n_volumes])
        eta[0,0,:] = fluid_properties.eta_V
        eta[0,1,:] = fluid_properties.eta_L
        V[0,0,:] = self.Vg
        V[0,1,:] = self.Vo
        moles_oil_and_gas = eta * V
        moles_w = fluid_properties.eta_W * self.Vw
        self.phase_mole_numbers[0,0:2,:] = moles_oil_and_gas
        self.phase_mole_numbers[0,2,:] = moles_w
        return eta, moles_oil_and_gas, moles_w
        #saem como um vetor em  função dos blocos

    def update_relative_permeabilities(self, fluid_properties):
        So, Sw, Sg = self.update_saturations_without_contours()
        saturations = np.array([So, Sg, Sw])
        kro,krg,krw = self.relative_permeability(saturations)
        self.relative_permeabilities = np.zeros([1,self.n_phases,len(self.internal_faces)])
        self.relative_permeabilities[0,0,:] = kro
        self.relative_permeabilities[0,1,:] = krg
        self.relative_permeabilities[0,2,:] = krw

    def update_phase_viscosities(self,fluid_properties):
        self.phase_viscosities = np.zeros(self.relative_permeabilities.shape)
        self.phase_viscosities_oil_and_gas = self.phase_viscosity(fluid_properties)
        self.phase_viscosities[0,0,:] = self.phase_viscosities_oil_and_gas[0,0,:]
        self.phase_viscosities[0,1,:] = self.phase_viscosities_oil_and_gas[0,1,:]
        self.phase_viscosities[0,2,:] = self.mi_W

    def dVt_derivatives(self, data_impress, fluid_properties):
        eta, No_g, Nw = self.update_phase_mole_numbers(fluid_properties, data_impress)
        dVtk = np.zeros([fluid_properties.Nc + 1, self.n_volumes])
        dVtk[0:fluid_properties.Nc,:], dVtP = PartialDerivatives().dVt_derivatives(
                fluid_properties, No_g, self.component_molar_fractions, eta)
        dVtk[fluid_properties.Nc,:] = 1 / fluid_properties.eta_W
        dVwP = np.zeros(self.n_volumes)
        dVtP = dVtP + dVwP
        return dVtk, dVtP

    def update_transmissibility(self, M, data_impress, data_loaded, elements_lv0, fluid_properties):
        v0 = elements_lv0['neig_internal_faces']
        pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
        self.pretransmissibility_internal_faces = pretransmissibility_faces[self.internal_faces]

        self.mobilities = self.relative_permeabilities / self.phase_viscosities
        self.dVtk, self.dVtP = self.dVt_derivatives(data_impress, fluid_properties)
        dVtk = np.delete(self.dVtk, self.all_wells, axis = 1)

        ''' Retracting the well volume (transmissibility without contours)'''
        self.component_molar_fractions = np.delete(self.component_molar_fractions, self.all_wells, axis=2)
        self.phase_molar_densities = np.delete(self.phase_molar_densities, self.all_wells, axis=2)

        ''' Transmissibility '''
        t0 = (dVtk * (self.component_molar_fractions * self.phase_molar_densities *\
              self.mobilities).sum(axis=1)).sum(axis=0)

        t0 = t0 * self.pretransmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape = (self.n_volumes, self.n_volumes))

        diag = np.diag(self.Vbulk * self.porosity * self.cf - self.dVtP)
        T += diag
        self['Tini'] = T
        return T


    def update_capillary_pressure(self, data_loaded):
        """ not working yet"""
        # get_capillary_pressure = getattr(capillary_pressure, data_loaded['compositional_data']['capillary_pressure'])
        # get_capillary_pressure = get_capillary_pressure(data_loaded, self.phase_molar_densities, self.component_molar_fractions,
        #                         self.pretransmissibility_internal_faces)
        # So, Sw, Sg = self.update_saturations_without_contours()
        # Pcow, Pcog = get_capillary_pressure(Sw, So, Sg)

        Pcap = np.zeros([2,3,self.n_volumes - len(self.all_wells)])
        # Pcap[0,0,:] = Pcog
        # Pcap[0,1,:] = Pcow
        return Pcap

    def update_deltaT(self):
        """ need to do that """
        deltaT=2
        return deltaT

    def pressure_independent_term(self, fluid_properties):
        vector = self.Vbulk * self.porosity * self.cf - self.dVtP
        pressure_term = vector * fluid_properties.P
        pressure_term = np.delete(pressure_term, self.all_wells)
        return pressure_term

    def capillary_independent_term(self, data_loaded):
        Pcap = self.update_capillary_pressure(data_loaded)
        self.dVtk = np.delete(self.dVtk,self.all_wells,axis = 1)
        capillary_term = np.sum(self.dVtk * np.sum (self.component_molar_fractions * self.phase_molar_densities *\
              self.mobilities * Pcap,axis =1), axis = 0)
        return capillary_term

    def volume_discrepancy_independent_term(self):
        self.Vp = np.delete(self.Vp, self.all_wells)
        self.phase_mole_numbers = np.delete(self.phase_mole_numbers, self.all_wells, axis= 2)
        volume_discrepancy_term = self.Vp - np.sum(np.sum(self.phase_mole_numbers / self.phase_molar_densities, axis = 0), axis = 0)
        return volume_discrepancy_term

    def well_term(self):
        q = np.zeros(self.n_volumes - len(self.all_wells))
        well_term = np.sum(self.dVtk,axis=0) * q
        return well_term

    def update_independent_terms(self, fluid_properties, data_loaded):
        deltaT = self.update_deltaT()
        pressure_term = self.pressure_independent_term(fluid_properties)
        capillary_term = self.capillary_independent_term(data_loaded)
        volume_term = self.volume_discrepancy_independent_term()
        well_term = self.well_term()
        well_pressures = 13.1 #ver como obter isso
        independent_terms = np.zeros(self.n_volumes)
        independent_terms[self.all_wells] = well_pressures #?
        ind = np.arange(self.n_volumes)
        independent_terms[ind!=self.all_wells] = pressure_term - volume_term + deltaT * well_term #+ deltaT * capillary_term
        return independent_terms
