from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from ..utils import relative_permeability2, phase_viscosity, capillary_pressure
from .. import directories as direc
import numpy as np


class PropertiesCalc:
    def __init__(self, n_volumes):
        self.n_volumes = n_volumes
        self.relative_permeability = getattr(relative_permeability2, data_loaded['compositional_data']['relative_permeability'])
        self.relative_permeability = self.relative_permeability()
        self.phase_viscosity = getattr(phase_viscosity, data_loaded['compositional_data']['phase_viscosity'])

    def run_outside_loop(self, data_impress, wells, fprop, kprop):
        self.set_properties(fprop, kprop)
        self.update_porous_volume(data_impress, fprop)
        self.update_saturations(data_impress, fprop, kprop)
        self.set_initial_volume(fprop)
        self.set_initial_mole_numbers(fprop, kprop)
        self.update_relative_permeabilities(fprop, kprop)
        self.update_phase_viscosities(data_loaded, fprop, kprop)
        self.update_mobilities(fprop)

    def run_inside_loop(self, data_impress, wells, fprop, kprop):
        self.set_properties(fprop, kprop)
        self.update_porous_volume(data_impress, fprop)
        if kprop.load_w: self.update_water_saturation(data_impress, wells, fprop, kprop)
        self.update_saturations(data_impress, fprop, kprop)
        self.update_mole_numbers(fprop, kprop)
        self.update_total_volume(fprop)
        self.update_relative_permeabilities(fprop, kprop)
        self.update_phase_viscosities(data_loaded, fprop, kprop)
        self.update_mobilities(fprop)

    def set_properties(self, fprop, kprop):
        fprop.component_molar_fractions = np.zeros([kprop.n_components, kprop.n_phases, self.n_volumes])
        fprop.phase_molar_densities = np.zeros([1, kprop.n_phases, self.n_volumes])

        if kprop.load_k:
            fprop.component_molar_fractions[0:kprop.Nc,0,:] = fprop.x
            fprop.component_molar_fractions[0:kprop.Nc,1,:] = fprop.y

            fprop.phase_molar_densities[0,0,:] = fprop.ksi_L
            fprop.phase_molar_densities[0,1,:] = fprop.ksi_V

        fprop.component_molar_fractions[kprop.n_components-1, kprop.n_phases-1,:] = 1 #water molar fraction in water component
        if kprop.load_w: fprop.phase_molar_densities[0, kprop.n_phases-1,:] = fprop.ksi_W

    def set_initial_volume(self, fprop):
        self.Vo = fprop.Vp * fprop.So
        self.Vg = fprop.Vp * fprop.Sg
        self.Vw = fprop.Vp * fprop.Sw
        fprop.Vt = self.Vo +self.Vg + self.Vw

    def set_initial_mole_numbers(self, fprop, kprop):
        fprop.phase_mole_numbers = np.zeros([1, kprop.n_phases, self.n_volumes])

        if kprop.load_k:
            fprop.phase_mole_numbers[0,0,:] = fprop.ksi_L * self.Vo
            fprop.phase_mole_numbers[0,1,:] = fprop.ksi_V * self.Vg
        if kprop.load_w:
            fprop.phase_mole_numbers[0,kprop.n_phases-1,:] = fprop.ksi_W * self.Vw

        component_phase_mole_numbers = fprop.component_molar_fractions * fprop.phase_mole_numbers
        fprop.component_mole_numbers = np.sum(component_phase_mole_numbers, axis = 1)

    def update_porous_volume(self, data_impress, fprop):
        fprop.Vp = fprop.porosity * fprop.Vbulk * (1 + fprop.cf*(fprop.P - fprop.Pf))

    def update_saturations(self, data_impress, fprop, kprop):
        fprop.Sw = data_impress['saturation']
        if kprop.load_k:
            fprop.Sg = np.zeros(fprop.Sw.shape)
            fprop.Sg[fprop.V!=0] = (1 - fprop.Sw[fprop.V!=0]) * \
                (fprop.V[fprop.V!=0] / fprop.ksi_V[fprop.V!=0]) / \
                (fprop.V[fprop.V!=0] / fprop.ksi_V[fprop.V!=0] +
                fprop.L[fprop.V!=0] / fprop.ksi_L[fprop.V!=0] )
            fprop.Sg[fprop.V==0] = 0
            fprop.So = 1 - fprop.Sw - fprop.Sg
        else: fprop.So = np.zeros(len(fprop.Sw)); fprop.Sg = np.zeros(len(fprop.Sw))

    def update_mole_numbers(self, fprop, kprop):
        # este daqui foi criado separado pois, quando compressivel, o volume poroso pode
        #diferir do volume total, servindo como um termo de correção de erro na equação da pressão,
        #como se sabe. A equação do outside loop ela relaciona o volume total com o poroso, de modo
        #que eles nunca vão ser diferentes e um erro vai ser propagado por toda a simulação. Todavia,
        #ele funciona para o primeiro passo de tempo uma vez que a pressão não mudou e Vp = Vt ainda.

        fprop.phase_mole_numbers = np.zeros([1, kprop.n_phases, self.n_volumes])

        if kprop.load_k:
            fprop.phase_mole_numbers[0,0,:] = fprop.component_mole_numbers[kprop.Nc-1,:]/fprop.x[kprop.Nc-1,:]
            fprop.phase_mole_numbers[0,1,:] = fprop.phase_mole_numbers[0,0,:]*fprop.V/fprop.L
            #fprop.phase_mole_numbers[0,1,:] = fprop.component_mole_numbers[fprop.y!=0,:][0,:]/fprop.y[fprop.y!=0][0]

        if kprop.load_w:
            fprop.phase_mole_numbers[0,kprop.n_phases-1,:] = fprop.component_mole_numbers[kprop.n_components-1,:]
        else: fprop.mole_number_w = np.zeros(self.n_volumes)

    def update_total_volume(self, fprop):
        fprop.Vt = np.sum(fprop.phase_mole_numbers / fprop.phase_molar_densities, axis = 1).ravel()

    def update_relative_permeabilities(self, fprop, kprop):
        Sgr = float(direc.data_loaded['compositional_data']['residual_saturations']['Sgr'])
        Swr = float(direc.data_loaded['compositional_data']['residual_saturations']['Swr'])

        saturations = np.array([fprop.So, fprop.Sg, fprop.Sw])
        kro,krg,krw, Sor = self.relative_permeability(saturations)
        self.relative_permeabilities = np.zeros([1, kprop.n_phases, self.n_volumes])
        if kprop.load_k:
            self.relative_permeabilities[0,0,:] = kro
            self.relative_permeabilities[0,1,:] = krg
        if kprop.load_w:
            self.relative_permeabilities[0, kprop.n_phases-1,:] = krw
            #if kprop.load_k:
                #if any(fprop.Sw > (1 - Sor - Sgr)) or any(fprop.So > (1 - Swr - Sgr)):
                #    raise ValueError('valor errado da saturacao - mudar delta_t_ini')

    def update_phase_viscosities(self, data_loaded, fprop, kprop):
        self.phase_viscosities = np.zeros(self.relative_permeabilities.shape)
        if kprop.load_k:
            self.phase_viscosity = self.phase_viscosity(self.n_volumes, fprop, kprop)
            #self.phase_viscosities[0,0:2,:] = 0.02*np.ones([2,self.n_volumes]) #only for BL test
            self.phase_viscosities_oil_and_gas = self.phase_viscosity(fprop, kprop)
            self.phase_viscosities[0,0:kprop.n_phases-1*kprop.load_w,:] = self.phase_viscosities_oil_and_gas
        if kprop.load_w:
            self.phase_viscosities[0,kprop.n_phases-1,:] = data_loaded['compositional_data']['water_data']['mi_W']

    def update_mobilities(self, fprop):
        fprop.mobilities = self.relative_permeabilities / self.phase_viscosities

    def update_water_saturation(self, data_impress, wells, fprop, kprop):
        fprop.ksi_W = fprop.ksi_W0 * (1 + fprop.Cw * (fprop.P - fprop.Pw))
        fprop.rho_W = fprop.ksi_W * fprop.Mw_w
        data_impress['saturation'] = fprop.component_mole_numbers[kprop.n_components-1,:] * (1 / fprop.ksi_W) / fprop.Vp #or Vt ?
