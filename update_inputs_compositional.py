from packs.utils.info_manager import InfoManager
from packs.compositional.properties_calculation import PropertiesCalc
from packs.compositional.stability_check import StabilityCheck
from packs.directories import data_loaded
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


class ComponentProperties:
    def __init__(self):
        self.load_k = data_loaded['hidrocarbon_components']
        self.load_w = data_loaded['water_component']
        self.compressible_k = data_loaded['compressible_fluid']
        self.n_phases = 2 * self.load_k + 1 * self.load_w

        if self.load_k:
            self.w = np.array(data_loaded['compositional_data']['component_data']['w']).astype(float)
            self.Bin = np.array(data_loaded['compositional_data']['component_data']['Bin']).astype(float)
            self.Tc = np.array(data_loaded['compositional_data']['component_data']['Tc']).astype(float)
            self.Pc = np.array(data_loaded['compositional_data']['component_data']['Pc']).astype(float)
            self.vc = np.array(data_loaded['compositional_data']['component_data']['vc']).astype(float)
            self.Mw = np.array(data_loaded['compositional_data']['component_data']['Mw']).astype(float)
            self.s = np.array(data_loaded['compositional_data']['component_data']['vshift_parameter']).astype(float)
            self.z = np.array(data_loaded['compositional_data']['component_data']['z']).astype(float)
            self.Nc = len(self.z)
        else: self.Nc = 0; self.z = []
        self.n_components = self.Nc + 1 * self.load_w

class FluidProperties:
    def __init__(self, kprop, n_volumes):
        P = np.array(data_loaded['Pressure']['r1']['value']).astype(float)
        self.P = P*np.ones(n_volumes)
        self.T = np.array(data_loaded['Temperature']['r1']['value'])
        self.component_molar_fractions = np.zeros([kprop.n_components, kprop.n_phases, n_volumes])

    def run_inputs_k(self, fprop_block, kprop, n_volumes):
        self.inputs_all_volumes_k(fprop_block, kprop, n_volumes)
        if not kprop.load_w: self.Cw = 0

    def inputs_all_volumes_k(self, fprop_block, kprop, n_volumes):
        self.T = fprop_block.T
        self.z = kprop.z * np.ones(n_volumes)
        self.x = fprop_block.x[:,np.newaxis] * np.ones([kprop.Nc, n_volumes])
        self.y = fprop_block.y[:,np.newaxis] * np.ones([kprop.Nc, n_volumes])
        self.L = fprop_block.L * np.ones(n_volumes)
        self.V = fprop_block.V * np.ones(n_volumes)
        '''self.Mw_L = fprop_block.Mw_L * np.ones(n_volumes)
        self.Mw_V = fprop_block.Mw_V * np.ones(n_volumes)
        self.ksi_L = fprop_block.ksi_L * np.ones(n_volumes)
        self.ksi_V = fprop_block.ksi_V * np.ones(n_volumes)
        self.rho_L = fprop_block.rho_L * np.ones(n_volumes)
        self.rho_V = fprop_block.rho_V * np.ones(n_volumes)'''

        #coef_vc7 = np.array([21.573, 0.015122, -27.6563, 0.070615])
        #kprop.vc[kprop.C7 == 1] =  coef_vc7[0] + coef_vc7[1] * np.mean(kprop.Mw[kprop.C7 == 1]) + \
                            #coef_vc7[2] * np.mean([kprop.SG[kprop.C7 == 1]]) + coef_vc7[3] \
                            #* np.mean(kprop.Mw[kprop.C7 == 1]) * np.mean([kprop.SG[kprop.C7 == 1]])

    def update_all_volumes(self, fprop_block, i):
        self.component_molar_fractions[:,0,i] = fprop_block.x
        self.component_molar_fractions[:,1,i] = fprop_block.y
        self.L[i] = fprop_block.L
        self.V[i] = fprop_block.V
        '''self.Mw_L[i] = fprop_block.Mw_L
        self.Mw_V[i] = fprop_block.Mw_V
        self.ksi_L[i] = fprop_block.ksi_L
        self.ksi_V[i] = fprop_block.ksi_V
        self.rho_L[i] = fprop_block.rho_L
        self.rho_V[i] = fprop_block.rho_V'''
