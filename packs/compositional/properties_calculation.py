from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from .. import directories as direc
import numpy as np


class PropertiesCalc:
    def __init__(self, data_impress, wells, fprop):
        self.n_phases = 3 #includding water
        self.n_volumes = data_impress.len_entities['volumes']
        self.n_components = fprop.Nc + 1
        self.set_properties(fprop)
        self. update_porous_volume(data_impress, fprop)

    def run_outside_loop(self, data_impress, wells, fprop):
        self.update_saturations(data_impress, fprop)
        self.update_mole_numbers(fprop)

    def run_inside_loop(self, data_impress, wells, fprop):
        self.update_water_saturation(data_impress, fprop)
        self.update_saturations(data_impress, fprop)
        self.update_mole_numbers(fprop)

    def update_porous_volume(self, data_impress, fprop):
        Pf = np.array(data_loaded['compositional_data']['Pf']).astype(float)
        Vbulk = data_impress['volume']
        porosity = data_impress['poro']
        cf = np.array(data_loaded['compositional_data']['rock_compressibility']).astype(float)
        # Pf - reference pressure at witch porosity is obtained
        fprop.Vp = porosity * Vbulk * (1 + cf*(fprop.P - Pf))

    def set_properties(self, fprop):
        fprop.component_molar_fractions = np.zeros([self.n_components, self.n_phases, self.n_volumes])
        fprop.phase_molar_densities = np.zeros([1, self.n_phases, self.n_volumes])

        fprop.component_molar_fractions[0:fprop.Nc,0,:] = fprop.x
        fprop.component_molar_fractions[0:fprop.Nc,1,:] = fprop.y
        fprop.component_molar_fractions[fprop.Nc,2,:] = 1 #water molar fraction in water component

        fprop.phase_molar_densities[0,0,:] = fprop.ksi_L
        fprop.phase_molar_densities[0,1,:] = fprop.ksi_V
        fprop.phase_molar_densities[0,2,:] = fprop.ksi_W


    def update_saturations(self, data_impress, fprop):
        fprop.Sw = data_impress['saturation']
        fprop.Sg = np.zeros(fprop.Sw.shape)
        fprop.Sg[fprop.V!=0] = (1 - fprop.Sw[fprop.V!=0]) * \
            (fprop.V[fprop.V!=0] / fprop.rho_V[fprop.V!=0]) / \
            (fprop.V[fprop.V!=0] / fprop.rho_V[fprop.V!=0] +
            fprop.L[fprop.V!=0] / fprop.rho_L[fprop.V!=0] )
        fprop.Sg[fprop.V==0] = 0
        fprop.So = 1 - fprop.Sw - fprop.Sg

    def update_phase_volumes(self, fprop):
        fprop.Vo = fprop.Vp * fprop.So
        fprop.Vg = fprop.Vp * fprop.Sg
        fprop.Vw = fprop.Vp * fprop.Sw
        fprop.Vt = fprop.Vo + fprop.Vg + fprop.Vw  #np.sum(fprop.phase_mole_numbers / fprop.phase_molar_densities, axis = 1).ravel()

    def update_mole_numbers(self, fprop):
        self.update_phase_volumes(fprop)
        fprop.phase_mole_numbers = np.zeros([1, self.n_phases, self.n_volumes])
        fprop.ksi_o_and_g = np.ones([1, 2, self.n_volumes])
        V = np.ones([1, 2, self.n_volumes])
        fprop.ksi_o_and_g[0,0,:] = fprop.ksi_V
        fprop.ksi_o_and_g[0,1,:] = fprop.ksi_L
        V[0,0,:] = fprop.Vg
        V[0,1,:] = fprop.Vo
        fprop.mole_numbers_o_and_g = fprop.ksi_o_and_g * V
        fprop.mole_number_w = fprop.ksi_W * fprop.Vw
        fprop.phase_mole_numbers[0,0:2,:] = fprop.mole_numbers_o_and_g
        fprop.phase_mole_numbers[0,2,:] = fprop.mole_number_w

        #saem como um vetor em  função dos blocos - rever isso daqui. Se eu atualizo o número de moles do componente,
        # não faz sentido eu calcular ele sempre desse jeito, ele teria que vim como um dado de entrada em t=0 e depois
        # essa conta some
    def update_water_saturation(self, data_impress, fprop):
        data_impress['saturation'] = fprop.component_mole_numbers[self.n_components-1,:] * (1 / fprop.ksi_W) / fprop.Vp
