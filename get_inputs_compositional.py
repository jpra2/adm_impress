from packs.utils.info_manager import InfoManager
from packs.directories import data_loaded
from packs.utils import constants as ctes
import numpy as np
import os
import yaml

with open ('input_cards/input_file_name.yml','r') as f:
    names_files_load = yaml.safe_load(f)
    name_input_file_load = names_files_load['name_file']

dd = InfoManager(name_input_file_load, 'input_cards/inputs0.yml')
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

class FluidProperties:
    def __init__(self, M, wells):
        self.P = np.array([data_loaded['Pressure']['r1']['value']]).astype(float)
        self.T = data_loaded['Temperature']['r1']['value']
        self.xkj = np.ones([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        self.Csi_j = np.empty([1, ctes.n_phases, ctes.n_volumes])
        self.rho_j = np.empty_like(self.Csi_j)
        self.Nk = np.empty([ctes.n_components, ctes.n_volumes])
        self.Nj = np.zeros([1, ctes.n_phases, ctes.n_volumes])
        self.update_initial_porous_volume()
        self.P = self.P * np.ones(ctes.n_volumes)
        self.P[wells['ws_p']] = wells['values_p']
        self.Sw = M.data['saturation']
        if ctes.load_k:
            self.z = np.array([data_loaded['compositional_data']['component_data']['z']]).astype(float).T
            self.z = self.z * np.ones(ctes.n_volumes)
            self.L = np.empty(len(self.P))
            self.V = np.empty(len(self.P))
        else: z = []

    def update_initial_porous_volume(self):
        self.Vp = ctes.porosity * ctes.Vbulk * (1 + ctes.Cf*(self.P - ctes.Pf))

    '''
    def inputs_fluid_properties(self):
        #coef_vc7 = np.array([21.573, 0.015122, -27.6563, 0.070615])
        #ctes.vc[ctes.C7 == 1] =  coef_vc7[0] + coef_vc7[1] * np.mean(ctes.Mw[ctes.C7 == 1]) + \
                            #coef_vc7[2] * np.mean([ctes.SG[ctes.C7 == 1]]) + coef_vc7[3] \
                            #* np.mean(ctes.Mw[ctes.C7 == 1]) * np.mean([ctes.SG[ctes.C7 == 1]])
    '''

    def inputs_water_properties(self, M):
        self.rho_j[0,ctes.n_phases-1,:] = data_loaded['compositional_data']['water_data']['rho_W']
        self.Csi_W0 = self.rho_j[0,ctes.n_phases-1,:] / ctes.Mw_w
        self.Csi_W = self.Csi_W0
        self.rho_W = self.Csi_W * ctes.Mw_w
        #self.Csi_j[0,ctes.n_phases-1,:] = self.Csi_W0
        self.xkj[ctes.n_components-1,ctes.n_phases-1,:] = 1
        self.xkj[ctes.n_components-1,0:ctes.n_phases-1,:] = 0
        self.xkj[0:ctes.n_components-1,ctes.n_phases-1,:] = 0
        self.Nk[-1,:] = self.Vp * self.Csi_W * self.Sw
