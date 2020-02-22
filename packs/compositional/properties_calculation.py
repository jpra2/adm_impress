from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from .. import directories as direc
import numpy as np


class PropertiesCalc(DataManager):
    def __init__(self, M, data_impress, wells, fluid_properties, load, data_name: str='Compositional_data.npz'):
        super().__init__(data_name, load=load)
        self.n_phases = 3 #includding water
        self.n_volumes = data_impress.len_entities['volumes']
        self.n_components = fluid_properties.Nc + 1
        self.run(data_loaded, data_impress, wells, fluid_properties )

    def run(self, data_loaded, data_impress, wells, fluid_properties):
        self.set_properties(fluid_properties)
        self.update_saturations(data_impress, wells, fluid_properties)
        self.update_mole_numbers(data_impress, fluid_properties)

    def set_properties(self, fluid_properties):
        fluid_properties.component_molar_fractions = np.zeros([self.n_components, self.n_phases, self.n_volumes])
        fluid_properties.phase_molar_densities = np.zeros([1, self.n_phases, self.n_volumes])


        fluid_properties.component_molar_fractions[0:fluid_properties.Nc,0,:] = fluid_properties.x
        fluid_properties.component_molar_fractions[0:fluid_properties.Nc,1,:] = fluid_properties.y
        fluid_properties.component_molar_fractions[fluid_properties.Nc,2,:] = 1 #water molar fraction in water component

        fluid_properties.phase_molar_densities[0,0,:] = fluid_properties.ksi_L
        fluid_properties.phase_molar_densities[0,1,:] = fluid_properties.ksi_V
        fluid_properties.phase_molar_densities[0,2,:] = fluid_properties.ksi_W


    def update_saturations(self, data_impress, wells, fluid_properties):
        fluid_properties.Sw = data_impress['saturation']
        fluid_properties.Sg = np.zeros(fluid_properties.Sw.shape)
        fluid_properties.Sg[fluid_properties.V!=0] = (1 - fluid_properties.Sw[fluid_properties.V!=0]) * \
            (fluid_properties.V[fluid_properties.V!=0] / fluid_properties.rho_V[fluid_properties.V!=0]) / \
            (fluid_properties.V[fluid_properties.V!=0] / fluid_properties.rho_V[fluid_properties.V!=0] +
            fluid_properties.L[fluid_properties.V!=0] / fluid_properties.rho_L[fluid_properties.V!=0] )
        fluid_properties.Sg[fluid_properties.V==0] = 0
        fluid_properties.So = 1 - fluid_properties.Sw - fluid_properties.Sg

    def update_phase_volumes(self, data_impress, fluid_properties):
        Pf = np.array(data_loaded['compositional_data']['Pf']).astype(float)
        Vbulk = data_impress['volume']
        porosity = data_impress['poro']
        cf = np.array(data_loaded['compositional_data']['rock_compressibility']).astype(float)
        # Pf - reference pressure at witch porosity is obtained
        fluid_properties.Vp = porosity * Vbulk * (1 + cf*(fluid_properties.P - Pf))
        fluid_properties.Vo = fluid_properties.Vp * fluid_properties.So
        fluid_properties.Vg = fluid_properties.Vp * fluid_properties.Sg
        fluid_properties.Vw = fluid_properties.Vp * fluid_properties.Sw

    def update_mole_numbers(self, data_impress, fluid_properties):
        self.update_phase_volumes(data_impress, fluid_properties)
        fluid_properties.phase_mole_numbers = np.zeros([1, self.n_phases, self.n_volumes])
        fluid_properties.ksi_o_and_g = np.ones([1, 2, self.n_volumes])
        V = np.ones([1, 2, self.n_volumes])
        fluid_properties.ksi_o_and_g[0,0,:] = fluid_properties.ksi_V
        fluid_properties.ksi_o_and_g[0,1,:] = fluid_properties.ksi_L
        V[0,0,:] = fluid_properties.Vg
        V[0,1,:] = fluid_properties.Vo
        fluid_properties.mole_numbers_o_and_g = fluid_properties.ksi_o_and_g * V
        fluid_properties.mole_number_w = fluid_properties.ksi_W * fluid_properties.Vw
        fluid_properties.phase_mole_numbers[0,0:2,:] = fluid_properties.mole_numbers_o_and_g
        fluid_properties.phase_mole_numbers[0,2,:] = fluid_properties.mole_number_w

        #saem como um vetor em  função dos blocos - rever isso daqui. Se eu atualizo o número de moles do componente,
        # não faz sentido eu calcular ele sempre desse jeito, ele teria que vim como um dado de entrada em t=0 e depois
        # essa conta some
