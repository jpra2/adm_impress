from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from ..utils import relative_permeability2, phase_viscosity
from .. import directories as direc
import numpy as np

#Next step: por os parâmetros de entrada e entender como a class fluid_properties vai entrar
# e tentar rodar

class CompositionalTPFA(DataManager):
    def __init__(self,data_name: str='CompositionalTPFA.npz'):
        super().__init__(data_impress, load=load)
        self.biphasic_data = data_loaded['biphasic_data']
        self.relative_permeability = getattr(relative_permeability2, self.compositional_data['relative_permeability'])
        self.relative_permeability = self.relative_permeability2()
        self.n_blocks = len(elements_lv0['volumes']) #número de blocos da malha
        self.phase_viscosity = phase_viscosity(self.n_blocks)
        self.n_phases = 3 #len(relative_permeabilities[:,0,0])
        self.Nc = fluid_properties.Nc
        if not load:
            self.loop = 0
            self.vpi = 0.0
            self.t = 0.0
            self.contador_vtk = 0
            self.update_viscosities(T,fluid_properties)
            self.update_saturations()
            self.update_relative_permeability()
            #self.update_mobilities()
            self.update_transmissibility_ini()
        else:
            self.load_infos()
        # self.M = M
        # self.elements_lv0 = elements_lv0
        # self.relative_permeability = getattr(relative_permeability, self.compositional_data['relative_permeability'])
        # load = data_loaded['load_compositional_data']
        # self.mesh_name = os.path.join(direc.flying, 'compositional_')
        # self.all_compositional_results = self.get_empty_current_compositional_results()
        # self.solver = SolverSp()

    def set_properties(fluid_properties,elements_lv0):

        component_molar_fractions = np.zeros([self.Nc,self.n_phases,self.n_blocks])
        phase_mass_densities = np.zeros([1,self.n_phases,self.n_blocks])
        phase_molar_densities = np.copy(phase_mass_densities)
        phase_molar_fractions = np.copy(phase_mass_densities)
        phase_viscosities = np.copy(phase_mass_densities)

        component_molar_fractions[:,0,:] = fluid_properties.x
        component_molar_fractions[:,1,:] = fluid_properties.y
        component_molar_fractions[:,2,:] = fluid_properties.w

        phase_mass_densities[0,0,:] = fluid_properties.rho_L
        phase_mass_densities[0,1,:] = fluid_properties.rho_V
        phase_mass_densities[0,2,:] = fluid_properties.rho_W

        phase_molar_densities[0,0,:] = fluid_properties.neta_L
        phase_molar_densities[0,1,:] = fluid_properties.neta_V
        phase_molar_densities[0,2,:] = fluid_properties.neta_W

        phase_molar_fractions[0,0,:] = fluid_properties.L
        phase_molar_fractions[0,1,:] = fluid_properties.V
        phase_molar_fractions[0,2,:] = fluid_properties.W

        phase_viscosities[0,0,:] = fluid_properties.mi_L
        phase_viscosities[0,1,:] = fluid_properties.mi_V
        phase_viscosities[0,2,:] = fluid_properties.mi_W


        fluid_properties.component_molar_fractions = component_molar_fractions
        fluid_properties.phase_mass_densities = phase_mass_densities
        fluid_properties.phase_molar_fractions = phase_molar_fractions
        fluid_properties.phase_molar_densities = phase_molar_densities
        fluid_properties.phase_viscosities = phase_viscosities

    def update_viscosities(self,T,fluid_properties):
        phase_viscosities = self.phase_viscosity(T,fluid_properties)

    def update_saturations(self,fluid_properties):
        Sw = self.data_impress['saturation']
        Sg = (1 - Sw) * (fluid_properties.V / fluid_properties.rho_V) / \
            (fluid_properties.V / fluid_properties.rho_V +
            fluid_properties.L / fluid_properties.rho_L )
        So = 1 - Sw - Sg
        saturations = np.zeros([self.Nc,self.n_phases,self.n_blocks])
        saturations[0,0,:] = So
        saturations[0,1,:] = Sg
        saturations[0,2,:] = Sw
        self.data_impress['saturation'] = saturations

    def update_relative_permeability(self, fluid_properties):
        kro,krg,krw = self.relative_permeability(self.data_impress['saturation'])
        relative_permeabilities = np.zeros([1,self.n_phases,self.n_blocks])
        relative_permeabilities[0,0,:] = kro
        relative_permeabilities[0,1,:] = krg
        relative_permeabilities[0,2,:] = krw
        fluid_properties.relative_permeabilities = relative_permeabilities

    def update_transmissibility(self,M,data_impress,elements_lv0,fluid_properties):
        v0 = elements_lv0['neig_internal_faces']
        internal_faces = elements_lv0['internal_faces']
        pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
        pretransmissibility_internal_faces = transmissibility_faces[internal_faces]

        mobilities = fluid_properties.relative_permeabilities / fluid_properties.phase_viscosities

        t0 = (fluid_properties.component_molar_fractions *
                fluid_properties.phase_molar_densities * mobilities).sum(axis=1).sum(axis=1)

        t0 = t0 * pretransmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, self.n_volumes))
        self['Tini'] = T
        # Falta incluir o termo da derivada de V por P que soma na diagonal principal
        #e o termo da derivada de V por Nk que multiplica a transmissibilidade t0 (ele
        #entra no segundo somatório)
