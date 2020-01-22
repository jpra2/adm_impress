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

def inputs_components_properties(data_loaded_compositional):
    w = data_loaded_compositional['w']
    Bin = data_loaded_compositional['Bin']
    R = data_loaded_compositional['R']
    Tc = data_loaded_compositional['Tc']
    Pc = data_loaded_compositional['Pc']
    Vc = data_loaded_compositional['Vc']
    T = data_loaded_compositional['T']
    P = data_loaded_compositional['P']
    C7 = data_loaded_compositional['C7']
    Mw = data_loaded_compositional['Mw']
    z = data_loaded_compositional['z']
    z = np.array(z).astype(float)
    return w, Bin, R, Tc, Pc, Vc, T, P, Mw, C7, z
