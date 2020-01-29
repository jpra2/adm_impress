from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from ..utils import relative_permeability2, phase_viscosity
from .. import directories as direc
from .partial_derivatives import FugacityParcialDerivative
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
        self.n_phases = 3 #len(relative_permeabilities[:,0,0])

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
        self.update_transmissibility(M, data_impress, data_loaded, elements_lv0, fluid_properties)

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

        self.component_molar_fractions = np.zeros([fluid_properties.Nc+1, self.n_phases, len(self.internal_faces)])
        self.phase_mass_densities = np.zeros([1, self.n_phases, len(self.internal_faces)])
        self.phase_molar_densities = np.copy(self.phase_mass_densities)

        self.component_molar_fractions[:,0,:] = fluid_properties.x
        self.component_molar_fractions[:,1,:] = fluid_properties.y
        self.component_molar_fractions[fluid_properties.Nc,2,:] = 1 #water molar fraction in water component

        self.phase_mass_densities[0,0,:] = fluid_properties.rho_L
        self.phase_mass_densities[0,1,:] = fluid_properties.rho_V
        self.phase_mass_densities[0,2,:] = fluid_properties.rho_W

        self.phase_molar_densities[0,0,:] = fluid_properties.eta_L
        self.phase_molar_densities[0,1,:] = fluid_properties.eta_V
        self.phase_molar_densities[0,2,:] = fluid_properties.eta_W


    def update_saturations(self, data_impress, wells, fluid_properties):
        Sw = data_impress['saturation']
        # ind = np.arange(len(Sw))
        self.all_wells = wells['all_wells']
        # index_wells = np.argwhere(Sw == Sw[all_wells])
        self.Sw = np.delete(Sw, self.all_wells)

        if fluid_properties.V != 0:
            self.Sg = (1 - self.Sw) * (fluid_properties.V / fluid_properties.rho_V) / \
                (fluid_properties.V / fluid_properties.rho_V +
                fluid_properties.L / fluid_properties.rho_L )
        else: self.Sg = 0
        self.So = 1 - self.Sw - self.Sg


    def update_phase_volumes(self, fluid_properties, data_impress):
        cf = data_loaded['compositional_data']['rock_compressibility']
        Pf = data_loaded['compositional_data']['Pf']
        Pf = np.array(Pf).astype(float)
        Vbulk = data_impress['volume']
        Vbulk = np.delete(Vbulk, self.all_wells)
        porosities = data_loaded['compositional_data']['porosity']
        self.Vp = porosities * Vbulk * (1 + cf*(fluid_properties.P - Pf))
        self.Vo = self.Vp * self.So
        self.Vg = self.Vp * self.Sg
        self.Vw = self.Vp * self.Sw
        import pdb; pdb.set_trace()

    def update_phase_mole_numbers(self, fluid_properties, data_impress):
        self.update_phase_volumes(fluid_properties, data_impress)
        eta = np.ones([1,2,len(fluid_properties.eta_V)])
        V = np.ones([1,2,len(self.Vo)])
        eta[0,0,:] = fluid_properties.eta_V
        eta[0,1,:] = fluid_properties.eta_L
        V[0,0,:] = self.Vg
        V[0,1,:] = self.Vo
        self.Nphase = eta * V
        self.Nw = fluid_properties.eta_W * self.Vw
        #saem como um vetor em  função dos blocos

    def update_relative_permeabilities(self, fluid_properties):
        saturations = np.array([self.So, self.Sg, self.Sw])
        # self.data_impress['saturation'] = saturations
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

    def dVtdNk(self, data_impress, fluid_properties):
        self.update_phase_mole_numbers(fluid_properties, data_impress)

        for b in range(len(data_impress['volume'])):
            Nphase = self.Nphase[0,:,b]
            matrix = FugacityParcialDerivative().dlnphi_dn_matrix_numerically(fluid_properties, Nphase)
            for k in range(0, self.Nc):
                dlnphi_dnk = FugacityParcialDerivative.get_dlnphi_dnk_numerically(fluid_properties, fluid_properties.x[k], Nphase, 1)
                dnkx_dNk = np.linalg.inv(matrix)@dlnphi_dnk
                dnky_dNk = - dnkx_dNk
                dnky_dNk[k] = 1 - dnkx_dNk
            import pdb; pdb.set_trace()
        # dvl_dn #falta
        # dVxdNk = sum(dnkx_dNk*(1/fluid_properties.eta_L + self.Nphase[1]*dvl_dn))
        return dVtdNk

    def dVtdNw(self):
        return dVtdNw

    def dVtdP(self):
        """é o mesmo fumo do de cima"""
        return dVtdP

    def update_deltaT(self):
        """include CFL condition and flux calculations"""

    def update_transmissibility(self, M, data_impress, data_loaded, elements_lv0, fluid_properties):
        v0 = elements_lv0['neig_internal_faces']
        pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
        pretransmissibility_internal_faces = pretransmissibility_faces[self.internal_faces]
        n_volumes = data_impress.len_entities['volumes']
        mobilities = self.relative_permeabilities / self.phase_viscosities
        dVtdNk = self.dVtdNk(data_impress, fluid_properties)
        import pdb; pdb.set_trace()
        dVtdP = self.dVtdP()
        porosity = data_loaded['compositional_data']['Porosity']
        cf = data_loaded['compositional_data']['rock_compressibility']

        t0 = (dVtdNk*(self.component_molar_fractions * self.phase_molar_densities *\
              mobilities).sum(axis=1)).sum(axis=0)

        t0 = t0 * pretransmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape = (n_volumes, n_volumes))

        # diagonal matrix:
        # see from where I get the bulk volume (sum of the fluid and rock and volumes
        # i guess the impress already has that (it can be calculated with the mesh
        # dimensions))
        # diag = np.diag(V_bulk * porosity * cf - dVtdP)
        # T += diag
        self['Tini'] = T

    # def independent_terms(self):
    #     Vp = V_bulk * porosity * (1 + cf * (P - self.P))
    #     Q = Vp - sum(phase_moles_number/self.phase_molar_densities)
    #
    # def solve_pressure(self):
    #     """pressure solver"""






        # Falta incluir o termo da derivada de V por P que soma na diagonal principal
        #e o termo da derivada de V por Nk que multiplica a transmissibilidade t0 (ele
        #entra no segundo somatório)
