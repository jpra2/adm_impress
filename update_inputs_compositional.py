from packs.utils.info_manager import InfoManager
from packs.compositional.properties_calculation import PropertiesCalc
from packs.compositional.stability_check import StabilityCheck
import numpy as np
import os

dd = InfoManager('input_cards/inputs_compositional.yml', 'input_cards/inputs0.yml')
dd2 = InfoManager('input_cards/variable_inputs_compositional.yml','input_cards/variable_input.yml')
dd['load_data'] = True
dd.save_obj()
dd2.save_obj()

if dd['deletar_results']:

    results = 'results'
    ff = os.listdir(results)

    for f in ff:
        if f[-4:] == '.vtk':
            os.remove(os.path.join(results, f))

def inputs_overall_properties(data_loaded):
    z = np.array(data_loaded['compositional_data']['component_data']['z']).astype(float)
    T = np.array(data_loaded['Temperature']['r1']['value']).astype(float)
    P = np.array(data_loaded['Pressure']['r1']['value']).astype(float)
    Nc = len(z)
    R = np.array(data_loaded['compositional_data']['component_data']['R']).astype(float)
    return z, P, T, R, Nc

class ComponentProperties:
    def __init__(self, data_loaded):
        self.w = np.array(data_loaded['compositional_data']['component_data']['w']).astype(float)
        self.Bin = np.array(data_loaded['compositional_data']['component_data']['Bin']).astype(float)
        self.Tc = np.array(data_loaded['compositional_data']['component_data']['Tc']).astype(float)
        self.Pc = np.array(data_loaded['compositional_data']['component_data']['Pc']).astype(float)
        self.vc = np.array(data_loaded['compositional_data']['component_data']['vc']).astype(float)
        self.Mw = np.array(data_loaded['compositional_data']['component_data']['Mw']).astype(float)
        self.C7 = np.array(data_loaded['compositional_data']['component_data']['C7']).astype(float)
        self.SG = np.array(data_loaded['compositional_data']['component_data']['SG']).astype(float)

class FluidProperties:
    def __init__(self, fprop_block, data_loaded, n_volumes):
        self.run_inputs(fprop_block, data_loaded, n_volumes)

    def run_inputs(self, fprop_block, data_loaded, n_volumes):
        self.inputs_water_properties(data_loaded, n_volumes)
        self.inputs_all_volumes(fprop_block, n_volumes)

    def inputs_water_properties(self, data_loaded, n_volumes):
        self.rho_W = data_loaded['compositional_data']['water_data']['rho_W'] * np.ones(n_volumes)
        self.Mw_w = data_loaded['compositional_data']['water_data']['Mw_w'] * np.ones(n_volumes)
        self.ksi_W0 = self.rho_W/self.Mw_w
        self.ksi_W = self.ksi_W0
        #self.Sw = data_loaded['Saturation']['r1']['value'] * np.ones(n_volumes)

    def inputs_all_volumes(self, fprop_block, n_volumes):
        self.T = fprop_block.T
        self.R = fprop_block.R
        self.P = fprop_block.P * np.ones(n_volumes)
        self.Nc = len(fprop_block.z)
        self.z = fprop_block.z * np.ones([fprop_block.Nc, n_volumes])
        self.x = fprop_block.x * np.ones([fprop_block.Nc, n_volumes])
        self.y = fprop_block.y * np.ones([fprop_block.Nc, n_volumes])
        self.L = fprop_block.L * np.ones(n_volumes)
        self.V = fprop_block.V * np.ones(n_volumes)
        self.Mw_L = fprop_block.Mw_L * np.ones(n_volumes)
        self.Mw_V = fprop_block.Mw_V * np.ones(n_volumes)
        self.ksi_L = fprop_block.ksi_L * np.ones(n_volumes)
        self.ksi_V = fprop_block.ksi_V * np.ones(n_volumes)
        self.rho_L = fprop_block.rho_L * np.ones(n_volumes)
        self.rho_V = fprop_block.rho_V * np.ones(n_volumes)


    def inputs_missing_properties(self, kprop):
        #coef_vc7 = np.array([21.573, 0.015122, -27.6563, 0.070615])
        #kprop.vc[kprop.C7 == 1] =  coef_vc7[0] + coef_vc7[1] * np.mean(kprop.Mw[kprop.C7 == 1]) + \
                            #coef_vc7[2] * np.mean([kprop.SG[kprop.C7 == 1]]) + coef_vc7[3] \
                            #* np.mean(kprop.Mw[kprop.C7 == 1]) * np.mean([kprop.SG[kprop.C7 == 1]])

        self.component_phase_mole_numbers = self.component_molar_fractions * self.phase_mole_numbers
        self.component_mole_numbers = np.sum(self.component_phase_mole_numbers, axis = 1)

    def update_all_volumes(self, fprop_block, i):
        self.z[:,i] = fprop_block.z
        self.x[:,i] = fprop_block.x
        self.y[:,i] = fprop_block.y
        self.L[i] = fprop_block.L
        self.V[i] = fprop_block.V
        self.Mw_L[i] = fprop_block.Mw_L
        self.Mw_V[i] = fprop_block.Mw_V
        self.ksi_L[i] = fprop_block.ksi_L
        self.ksi_V[i] = fprop_block.ksi_V
        self.rho_L[i] = fprop_block.rho_L
        self.rho_V[i] = fprop_block.rho_V
