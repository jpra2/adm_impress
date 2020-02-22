from packs.utils.info_manager import InfoManager
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
        self.Vc = np.array(data_loaded['compositional_data']['component_data']['Vc']).astype(float)
        self.Mw = np.array(data_loaded['compositional_data']['component_data']['Mw']).astype(float)
        self.C7 = np.array(data_loaded['compositional_data']['component_data']['C7']).astype(float)

class FluidProperties:
    def __init__(self, fp, data_loaded, n_volumes):
        self.run_inputs(fp, data_loaded, n_volumes)

    def run_inputs(self, fp, data_loaded, n_volumes):
        self.inputs_water_properties(data_loaded, n_volumes)
        self.inputs_all_volumes(fp,n_volumes)

    def inputs_water_properties(self, data_loaded, n_volumes):
        self.rho_W = data_loaded['compositional_data']['water_data']['rho_W'] * np.ones(n_volumes)
        self.Mw_w = data_loaded['compositional_data']['water_data']['Mw_w'] * np.ones(n_volumes)
        self.ksi_W = self.rho_W/self.Mw_w

    def inputs_all_volumes(self, fp, n_volumes):
        self.T = fp.T
        self.R = fp.R
        self.P = fp.P * np.ones(n_volumes)
        self.Nc = len(fp.z)
        self.z = fp.z * np.ones([fp.Nc, n_volumes])
        self.x = fp.x * np.ones([fp.Nc, n_volumes])
        self.y = fp.y * np.ones([fp.Nc, n_volumes])
        self.L = fp.L * np.ones(n_volumes)
        self.V = fp.V * np.ones(n_volumes)
        self.Mw_L = fp.Mw_L * np.ones(n_volumes)
        self.Mw_V = fp.Mw_V * np.ones(n_volumes)
        self.ksi_L = fp.ksi_L * np.ones(n_volumes)
        self.ksi_V = fp.ksi_V * np.ones(n_volumes)
        self.rho_L = fp.rho_L * np.ones(n_volumes)
        self.rho_V = fp.rho_V * np.ones(n_volumes)

    def inputs_missing_properties(self):
        self.component_phase_mole_numbers = self.component_molar_fractions * self.phase_mole_numbers
        self.component_mole_numbers = np.sum(self.component_phase_mole_numbers, axis = 1)

    def update_all_volumes(self, fp, i):
        self.z[:,i] = fp.z
        self.x[:,i] = fp.x
        self.y[:,i] = fp.y
        self.L[i] = fp.L
        self.V[i] = fp.V
        self.Mw_L[i] = fp.Mw_L
        self.Mw_V[i] = fp.Mw_V
        self.ksi_L[i] = fp.ksi_L
        self.ksi_V[i] = fp.ksi_V
        self.rho_L[i] = fp.rho_L
        self.rho_V[i] = fp.rho_V
