from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from ..utils import relative_permeability2, phase_viscosity, capillary_pressure
from .. import directories as direc
import numpy as np


class PropertiesCalc:
    def __init__(self, data_impress, wells, fprop):
        self.n_phases = 3 #includding water
        self.n_volumes = data_impress.len_entities['volumes']
        self.n_components = fprop.Nc + 1
        self.set_properties(fprop)
        self. update_porous_volume(data_impress, fprop)

    def run_outside_loop(self, data_impress, wells, fprop, kprop):
        self.update_saturations(data_impress, fprop)
        self.update_mole_numbers(fprop)

        self.relative_permeability = getattr(relative_permeability2, data_loaded['compositional_data']['relative_permeability'])
        self.relative_permeability = self.relative_permeability()
        self.phase_viscosity = getattr(phase_viscosity, data_loaded['compositional_data']['phase_viscosity'])
        self.phase_viscosity = self.phase_viscosity(self.n_volumes, fprop, kprop)

        self.update_relative_permeabilities(fprop)
        self.update_phase_viscosities(data_loaded, fprop, kprop)

    def run_inside_loop(self, data_impress, wells, fprop, kprop):
        self.update_water_saturation(data_impress, fprop)
        self.update_saturations(data_impress, fprop)
        self.update_mole_numbers(fprop)

        self.relative_permeability = getattr(relative_permeability2, data_loaded['compositional_data']['relative_permeability'])
        self.relative_permeability = self.relative_permeability()
        self.phase_viscosity = getattr(phase_viscosity, data_loaded['compositional_data']['phase_viscosity'])
        self.phase_viscosity = self.phase_viscosity(self.n_volumes, fprop, kprop)

        self.update_relative_permeabilities(fprop)
        self.update_phase_viscosities(data_loaded, fprop, kprop)

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
            (fprop.V[fprop.V!=0] / fprop.ksi_V[fprop.V!=0]) / \
            (fprop.V[fprop.V!=0] / fprop.ksi_V[fprop.V!=0] +
            fprop.L[fprop.V!=0] / fprop.ksi_L[fprop.V!=0] )
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
        fprop.ksi_o_and_g[0,0,:] = fprop.ksi_L
        fprop.ksi_o_and_g[0,1,:] = fprop.ksi_V
        V[0,0,:] = fprop.Vo
        V[0,1,:] = fprop.Vg
        fprop.mole_numbers_o_and_g = fprop.ksi_o_and_g * V
        fprop.mole_number_w = fprop.ksi_W * fprop.Vw
        fprop.phase_mole_numbers[0,0:2,:] = fprop.mole_numbers_o_and_g
        fprop.phase_mole_numbers[0,2,:] = fprop.mole_number_w

        #saem como um vetor em  função dos blocos - rever isso daqui. Se eu atualizo o número de moles do componente,
        # não faz sentido eu calcular ele sempre desse jeito, ele teria que vim como um dado de entrada em t=0 e depois
        # essa conta some


    def update_relative_permeabilities(self, fprop):
        # So, Sw, Sg = self.update_saturations_without_contours()
        # saturations = np.array([So, Sg, Sw])
        saturations = np.array([fprop.So, fprop.Sg, fprop.Sw])
        kro,krg,krw = self.relative_permeability(saturations)
        fprop.relative_permeabilities = np.zeros([1, self.n_phases, self.n_volumes])
        fprop.relative_permeabilities[0,0,:] = kro
        fprop.relative_permeabilities[0,1,:] = krg
        fprop.relative_permeabilities[0,2,:] = krw

    def update_phase_viscosities(self, data_loaded, fprop, kprop):
        mi_W = data_loaded['compositional_data']['water_data']['mi_W']
        fprop.phase_viscosities = np.zeros(fprop.relative_permeabilities.shape)
        self.phase_viscosities_oil_and_gas = self.phase_viscosity(fprop, kprop)
        fprop.phase_viscosities[0,0,:] = self.phase_viscosities_oil_and_gas[0,0,:]
        fprop.phase_viscosities[0,1,:] = self.phase_viscosities_oil_and_gas[0,1,:]
        fprop.phase_viscosities[0,2,:] = mi_W

    def update_water_saturation(self, data_impress, fprop):
        Pw = np.array(data_loaded['compositional_data']['water_data']['Pw']).astype(float)
        Cw = np.array(data_loaded['compositional_data']['water_data']['Cw']).astype(float)
        fprop.ksi_W = fprop.ksi_W0*(1 + Cw * (fprop.P - Pw))
        fprop.rho_W = fprop.ksi_W * fprop.Mw_w
        #self.Sw = fprop.component_mole_numbers[self.n_components-1,:] * (1 / fprop.ksi_W) / fprop.Vp
        data_impress['saturation'] = fprop.component_mole_numbers[self.n_components-1,:] * (1 / fprop.ksi_W) / fprop.Vp
