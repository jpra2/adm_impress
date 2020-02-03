from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from ..utils import relative_permeability2, phase_viscosity
from .. import directories as direc
from .partial_derivatives import ParcialDerivatives
import scipy.sparse as sp
import numpy as np

#Next step: por os parâmetros de entrada e entender como a class fluid_properties vai entrar
# e tentar rodar
# MUDAR AS COISAS PARA TER COMPRIMENTO DO NÚMERO DE FACES INTERNAS - PELO MENOS PRA RESOLVER A transmissibilidade

class CompositionalTPFA(DataManager):
    def __init__(self, M, data_impress, wells, fluid_properties, elements_lv0, load, data_name: str='CompositionalTPFA.npz'):
        super().__init__(data_name, load=load)
        self.internal_faces = elements_lv0['internal_faces']
        self.relative_permeability = getattr(relative_permeability2, data_loaded['compositional_data']['relative_permeability'])
        self.relative_permeability = self.relative_permeability()
        self.phase_viscosity = getattr(phase_viscosity, data_loaded['compositional_data']['phase_viscosity'])
        self.phase_viscosity = self.phase_viscosity(len(self.internal_faces), fluid_properties)
        self.n_phases = 3 #includding water
        self.n_volumes = data_impress.len_entities['volumes']
        self.all_wells = wells['all_wells']
        self.Vbulk = data_impress['volume']
        self.porosity = data_loaded['compositional_data']['porosity']
        self.cf = data_loaded['compositional_data']['rock_compressibility']

        if not load:
            self.loop = 0
            self.vpi = 0.0
            self.t = 0.0
            self.contador_vtk = 0
            self.run(M, data_loaded, data_impress, wells, fluid_properties, elements_lv0)
        else: self.load_infos()

    def run(self, M, data_loaded, data_impress, wells, fluid_properties, elements_lv0):
        self.get_water_data(data_loaded,fluid_properties)
        self.set_properties(fluid_properties, elements_lv0)
        self.update_saturations(data_impress, wells, fluid_properties)
        self.update_relative_permeabilities(fluid_properties)
        self.update_phase_viscosities(fluid_properties)
        T = self.update_transmissibility(M, data_impress, data_loaded, elements_lv0, fluid_properties)
        D = self.update_independent_terms(fluid_properties)
        # self.M = M
        # self.elements_lv0 = elements_lv0
        # self.relative_permeability = getattr(relative_permeability, self.compositional_data['relative_permeability'])
        # load = data_loaded['load_compositional_data']
        # self.mesh_name = os.path.join(direc.flying, 'compositional_')
        # self.all_compositional_results = self.get_empty_current_compositional_results()
        # self.solver = SolverSp()

    def get_water_data(self, data_loaded, fluid_properties):
        ''' Attention: Water assumed to be incompressible '''
        fluid_properties.rho_W = data_loaded['compositional_data']['rho_W']
        self.mi_W = data_loaded['compositional_data']['mi_W']
        fluid_properties.Mw_w = data_loaded['compositional_data']['Mw_w']
        fluid_properties.eta_W = fluid_properties.rho_W/fluid_properties.Mw_w

    def set_properties(self, fluid_properties, elements_lv0): #provavelmente n vai estar aqui
        self.component_molar_fractions = np.zeros([fluid_properties.Nc+1, self.n_phases, self.n_volumes])
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
        else: self.Sg = np.zeros(len(self.Sw))
        self.So = 1 - self.Sw - self.Sg

    def update_phase_volumes(self, fluid_properties, data_impress):
        Pf = data_loaded['compositional_data']['Pf']
        Pf = np.array(Pf).astype(float)
        self.Vp = self.porosity * self.Vbulk * (1 + self.cf*(fluid_properties.P - Pf))
        self.Vo = self.Vp * self.So
        self.Vg = self.Vp * self.Sg
        self.Vw = self.Vp * self.Sw

    def update_phase_mole_numbers(self, fluid_properties, data_impress):
        self.update_phase_volumes(fluid_properties, data_impress)
        self.phase_mole_numbers = np.zeros([1, self.n_phases, self.n_volumes])
        eta = np.ones([1,2,len(fluid_properties.eta_V)])
        V = np.ones([1,2,len(self.Vo)])
        eta[0,0,:] = fluid_properties.eta_V
        eta[0,1,:] = fluid_properties.eta_L
        V[0,0,:] = self.Vg
        V[0,1,:] = self.Vo
        No_g = eta * V
        Nw = fluid_properties.eta_W * self.Vw
        self.phase_mole_numbers[0,0:2,:] = No_g
        self.phase_mole_numbers[0,2,:] = Nw
        return eta, No_g, Nw
        #saem como um vetor em  função dos blocos

    def update_relative_permeabilities(self, fluid_properties):
        So = np.delete(self.So, self.all_wells)
        Sw = np.delete(self.Sw, self.all_wells)
        Sg = np.delete(self.Sg, self.all_wells)
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
        dVtk[0:fluid_properties.Nc,:], dVtP = ParcialDerivatives().dVt_derivatives(
                fluid_properties, No_g, self.component_molar_fractions, eta)
        dVtk[fluid_properties.Nc,:] = 1 / fluid_properties.eta_W
        dVwP = np.zeros(self.n_volumes)
        dVtP = dVtP + dVwP
        return dVtk, dVtP

    def update_flux_volumes(self):
        "need to do that"

    def update_deltaT(self):
        ###
        ## de acordo com o fluxo nos volumes
        ###

        flux_volumes = np.absolute(self.data_impress['flux_volumes'])
        phis = self.data_impress['poro']
        volume = self.data_impress['volume']
        self.delta_t = (self.biphasic_data['cfl']*(volume*phis)/flux_volumes).min()

    def update_transmissibility(self, M, data_impress, data_loaded, elements_lv0, fluid_properties):
        v0 = elements_lv0['neig_internal_faces']
        pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
        pretransmissibility_internal_faces = pretransmissibility_faces[self.internal_faces]

        mobilities = self.relative_permeabilities / self.phase_viscosities
        self.dVtk, self.dVtP = self.dVt_derivatives(data_impress, fluid_properties)
        dVtk = np.delete(self.dVtk, self.all_wells, axis = 1)
        ''' Retracting the well volume (transmissibility without contours)'''
        component_molar_fractions = np.delete(self.component_molar_fractions, self.all_wells, axis=2)
        phase_molar_densities = np.delete(self.phase_molar_densities, self.all_wells, axis=2)

        t0 = (dVtk * (component_molar_fractions * phase_molar_densities *\
              mobilities).sum(axis=1)).sum(axis=0)

        t0 = t0 * pretransmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape = (self.n_volumes, self.n_volumes))

        # diagonal matrix:
        # see from where I get the bulk volume (sum of the fluid and rock and volumes
        # i guess the impress already has that (it can be calculated with the mesh
        # dimensions))
        diag = np.diag(self.Vbulk * self.porosity * self.cf - self.dVtP)
        T += diag
        self['Tini'] = T
        return T

    def update_independent_terms(self, fluid_properties):
        vec_Pn = self.Vbulk * self.porosity * self.cf - self.dVtP
        deltaT = 2
        q = np.zeros(self.n_volumes)
        deltaV = self.Vp - np.sum(np.sum(self.phase_mole_numbers / self.phase_molar_densities, axis = 0),axis=0)
        import pdb; pdb.set_trace()
        independent_terms = vec_Pn * fluid_properties.P - deltaV + deltaT * np.sum(self.dVtk,axis=0) * q
        return independent_terms
