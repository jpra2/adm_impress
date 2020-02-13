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

def inputs_components_properties(data_loaded, n_volumes):
    w = np.array(data_loaded['compositional_data']['component_data']['w']).astype(float)
    Bin = np.array(data_loaded['compositional_data']['component_data']['Bin']).astype(float)
    R = np.array(data_loaded['compositional_data']['component_data']['R']).astype(float)
    Tc = np.array(data_loaded['compositional_data']['component_data']['Tc']).astype(float)
    Pc = np.array(data_loaded['compositional_data']['component_data']['Pc']).astype(float)
    Vc = np.array(data_loaded['compositional_data']['component_data']['Vc']).astype(float)
    T = np.array(data_loaded['Temperature']['r1']['value']).astype(float)
    P = np.array(data_loaded['Pressure']['r1']['value']).astype(float)
    C7 = np.array(data_loaded['compositional_data']['component_data']['C7']).astype(float)
    Mw = np.array(data_loaded['compositional_data']['component_data']['Mw']).astype(float)
    z = np.array(data_loaded['compositional_data']['component_data']['z']).astype(float)
    P = P * np.ones(n_volumes)
    return w, Bin, R, Tc, Pc, Vc, T, P, Mw, C7, z
