from ..directories import data_loaded
from ..data_class.data_manager import DataManager
from ..utils import relative_permeability2, phase_viscosity, capillary_pressure
from .. import directories as direc
import numpy as np


class PropertiesCalc:
    def __init__(self, data_impress, n_volumes):
        self.n_volumes = n_volumes
        self.Pf = np.array(data_loaded['compositional_data']['Pf']).astype(float)
        self.Vbulk = data_impress['volume']
        self.porosity = data_impress['poro']
        self.cf = np.array(data_loaded['compositional_data']['rock_compressibility']).astype(float)

    def run_outside_loop(self, data_impress, wells, fprop, kprop):
        self.set_properties(fprop, kprop)
        self.update_porous_volume(data_impress, fprop)
        self.update_saturations(data_impress, fprop, kprop)
        self.update_mole_numbers_outside_loop(fprop, kprop)
        self.update_relative_permeabilities(fprop, kprop)
        fprop.phase_viscosities = np.zeros(fprop.relative_permeabilities.shape)
        self.update_phase_viscosities(data_loaded, fprop, kprop)
        self.set_missing_properties_outside_loop(kprop, fprop)


    def run_inside_loop(self, data_impress, wells, fprop, kprop):
        self.set_properties(fprop, kprop)
        self.update_porous_volume(data_impress, fprop)
        if kprop.load_w: self.update_water_saturation(data_impress, fprop, kprop)
        self.update_saturations(data_impress, fprop, kprop)
        self.update_mole_numbers_inside_loop(fprop, kprop)
        self.update_relative_permeabilities(fprop, kprop)
        self.update_phase_viscosities(data_loaded, fprop, kprop)

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

    def update_porous_volume(self, data_impress, fprop):
        fprop.Vp = self.porosity * self.Vbulk * (1 + self.cf*(fprop.P - self.Pf))

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

    def update_phase_volumes(self, fprop):
        fprop.Vo = fprop.Vp * fprop.So
        fprop.Vg = fprop.Vp * fprop.Sg
        fprop.Vw = fprop.Vp * fprop.Sw

    def update_mole_numbers_outside_loop(self, fprop, kprop):
        self.update_phase_volumes(fprop)
        fprop.phase_mole_numbers = np.zeros([1, kprop.n_phases, self.n_volumes])

        if kprop.load_k:
            fprop.phase_mole_numbers[0,0,:] = fprop.ksi_L*fprop.Vo
            fprop.phase_mole_numbers[0,1,:] = fprop.ksi_V*fprop.Vg
        if kprop.load_w:
            fprop.phase_mole_numbers[0,kprop.n_phases-1,:] = fprop.ksi_W * fprop.Vw
        fprop.Vt = fprop.Vo + fprop.Vg + fprop.Vw

    def update_mole_numbers_inside_loop(self, fprop, kprop):
        # este daqui foi criado separado pois, quando compressivel, o volume poroso pode
        #diferir do volume total, servindo como um termo de correção de erro na equação da pressão,
        #como se sabe. A equação do outside loop ela relaciona o volume total com o poroso, de modo
        #que eles nunca vão ser diferentes e um erro vai ser propagado por toda a simulação. Todavia,
        #ele funciona para o primeiro passo de tempo uma vez que a pressão não mudou e Vp = Vt ainda.
        self.update_phase_volumes(fprop)
        fprop.phase_mole_numbers = np.zeros([1, kprop.n_phases, self.n_volumes])

        if kprop.load_k:
            fprop.phase_mole_numbers[0,0,:] = fprop.component_mole_numbers[kprop.Nc-1,:]/fprop.x[kprop.Nc-1,:]
            fprop.phase_mole_numbers[0,1,:] = fprop.phase_mole_numbers[0,0,:]*fprop.V/fprop.L
            #fprop.phase_mole_numbers[0,1,:] = fprop.component_mole_numbers[fprop.y!=0,:][0,:]/fprop.y[fprop.y!=0][0]

        if kprop.load_w:
            fprop.phase_mole_numbers[0,kprop.n_phases-1,:] = fprop.component_mole_numbers[kprop.n_components-1,:]
        else: fprop.mole_number_w = np.zeros(self.n_volumes)

        fprop.Vt = np.sum(fprop.phase_mole_numbers / fprop.phase_molar_densities, axis = 1).ravel()

    def set_missing_properties_outside_loop(self, kprop, fprop):
        if kprop.load_w:
            mi_W = data_loaded['compositional_data']['water_data']['mi_W']
            fprop.phase_viscosities[0,kprop.n_phases-1,:] = mi_W

        fprop.component_phase_mole_numbers = fprop.component_molar_fractions * fprop.phase_mole_numbers
        fprop.component_mole_numbers = np.sum(fprop.component_phase_mole_numbers, axis = 1)


    def update_relative_permeabilities(self, fprop, kprop):
        self.relative_permeability = getattr(relative_permeability2, data_loaded['compositional_data']['relative_permeability'])
        self.relative_permeability = self.relative_permeability()

        saturations = np.array([fprop.So, fprop.Sg, fprop.Sw])
        kro,krg,krw = self.relative_permeability(saturations)
        fprop.relative_permeabilities = np.zeros([1, kprop.n_phases, self.n_volumes])
        if kprop.load_k:
            fprop.relative_permeabilities[0,0,:] = kro
            fprop.relative_permeabilities[0,1,:] = krg
        if kprop.load_w: fprop.relative_permeabilities[0, kprop.n_phases-1,:] = krw

    def update_phase_viscosities(self, data_loaded, fprop, kprop):
        if kprop.load_k:
            self.phase_viscosity = getattr(phase_viscosity, data_loaded['compositional_data']['phase_viscosity'])
            self.phase_viscosity = self.phase_viscosity(self.n_volumes, fprop, kprop)

            #fprop.phase_viscosities[0,0:2,:] = 0.02*np.ones([2,self.n_volumes]) #only for BL test
            self.phase_viscosities_oil_and_gas = self.phase_viscosity(fprop, kprop)
            fprop.phase_viscosities[0,0:kprop.n_phases-1*kprop.load_w,:] = self.phase_viscosities_oil_and_gas

    def update_water_saturation(self, data_impress, fprop, kprop):
        Pw = np.array(data_loaded['compositional_data']['water_data']['Pw']).astype(float)
        fprop.Cw = np.array(data_loaded['compositional_data']['water_data']['Cw']).astype(float)
        fprop.ksi_W = fprop.ksi_W0 * (1 + fprop.Cw * (fprop.P - Pw))
        fprop.rho_W = fprop.ksi_W * fprop.Mw_w
        data_impress['saturation'] = fprop.component_mole_numbers[kprop.n_components-1,:] * (1 / fprop.ksi_W) / fprop.Vp
