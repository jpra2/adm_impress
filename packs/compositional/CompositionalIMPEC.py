from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from ..utils import relative_permeability2, phase_viscosity
from .. import directories as direc
import numpy as np

#Next step: por os parâmetros de entrada e entender como a class fluid_properties vai entrar
# e tentar rodar

class CompositionalTPFA(DataManager):
    def __init__(self, M, fluid_properties, elements_lv0, load, data_name: str='CompositionalTPFA.npz'):
        super().__init__(data_impress, load=load)
        self.compositional_data = data_loaded['compositional_data']
        self.n_blocks = len(elements_lv0['volumes']) #número de blocos da malha
        self.relative_permeability = getattr(relative_permeability2, self.compositional_data['relative_permeability'])
        self.relative_permeability = self.relative_permeability2()
        self.phase_viscosity = getattr(phase_viscosity, self.compositional_data['phase_viscosity'])
        self.phase_viscosity = self.phase_viscosity(self.n_blocks, fluid_properties)
        self.n_phases = 3 #len(relative_permeabilities[:,0,0])
        self.Nc = self.compositional_data['Nc']
        if not load:
            self.loop = 0
            self.vpi = 0.0
            self.t = 0.0
            self.contador_vtk = 0
            self.run(M, fluid_properties, elements_lv0)
        else: self.load_infos()

    def run(self, M, fluid_properties, elements_lv0):
        self.get_water_data()
        self.set_properties(fluid_properties, elements_lv0)
        self.update_saturations(fluid_properties)
        self.update_relative_permeability(fluid_properties)
        self.update_transmissibility(M, elements_lv0, fluid_properties)

        # self.M = M
        # self.elements_lv0 = elements_lv0
        # self.relative_permeability = getattr(relative_permeability, self.compositional_data['relative_permeability'])
        # load = data_loaded['load_compositional_data']
        # self.mesh_name = os.path.join(direc.flying, 'compositional_')
        # self.all_compositional_results = self.get_empty_current_compositional_results()
        # self.solver = SolverSp()

    def get_water_data():
        self.rho_W = self.compositional_data['rho_W']
        self.mi_W = self.compositional_data['mi_W']
        self.Mw_w = self.compositional_data['Mw_w']
        self.eta_W = self.rho_W/self.Mw_w

    def set_properties(self, fluid_properties, elements_lv0): #provavelmente n vai estar aqui

        component_molar_fractions = np.zeros([self.Nc+1, self.n_phases, self.n_blocks])
        phase_mass_densities = np.zeros([1, self.n_phases, self.n_blocks])
        phase_molar_densities = np.copy(phase_mass_densities)

        self.component_molar_fractions[:,0,:] = fluid_properties.x
        self.component_molar_fractions[:,1,:] = fluid_properties.y
        self.component_molar_fractions[Nc,2,:] = 1 #water molar fraction in water component

        self.phase_mass_densities[0,0,:] = fluid_properties.rho_L
        self.phase_mass_densities[0,1,:] = fluid_properties.rho_V
        self.phase_mass_densities[0,2,:] = fluid_properties.rho_W

        self.phase_molar_densities[0,0,:] = fluid_properties.eta_L
        self.phase_molar_densities[0,1,:] = fluid_properties.eta_V
        self.phase_molar_densities[0,2,:] = fluid_properties.eta_W

        self.phase_viscosities = phase_viscosity(fluid_properties)

    def update_saturations(self,fluid_properties):
        self.Sw = self.data_impress['saturation']
        self.Sg = (1 - Sw) * (fluid_properties.V / fluid_properties.rho_V) / \
            (fluid_properties.V / fluid_properties.rho_V +
            fluid_properties.L / fluid_properties.rho_L )
        self.So = 1 - Sw - Sg
        # saturations = np.zeros([self.Nc,self.n_phases,self.n_blocks])
        # saturations[0,0,:] = So
        # saturations[0,1,:] = Sg
        # saturations[0,2,:] = Sw
        # self.data_impress['saturation'] = saturations

    def update_relative_permeability(self, fluid_properties):
        kro,krg,krw = self.relative_permeability(self.So, self.Sg, self.Sw)
        self.relative_permeabilities = np.zeros([1,self.n_phases,self.n_blocks])
        self.relative_permeabilities[0,0,:] = kro
        self.relative_permeabilities[0,1,:] = krg
        self.relative_permeabilities[0,2,:] = krw

    def update_transmissibility(self, M, elements_lv0, fluid_properties):
        v0 = elements_lv0['neig_internal_faces']
        internal_faces = elements_lv0['internal_faces']
        pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
        pretransmissibility_internal_faces = transmissibility_faces[internal_faces]

        mobilities = self.relative_permeabilities / self.phase_viscosities
        print('Yay')
        t0 = (self.component_molar_fractions * self.phase_molar_densities *
              mobilities).sum(axis=1).sum(axis=1)

        t0 = t0 * pretransmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape = (self.n_volumes, self.n_volumes))
        self['Tini'] = T
        # Falta incluir o termo da derivada de V por P que soma na diagonal principal
        #e o termo da derivada de V por Nk que multiplica a transmissibilidade t0 (ele
        #entra no segundo somatório)
