import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

descriptions = {
    'case59': 'Cr7',
    'case60': 'Cr5'
}

def get_file_by_name(all_names, file_name):
    arqs = pd.Series(all_names)
    test = arqs.str.contains(file_name)
    result = arqs.values[test.values]
    return result

def load_file_in_flying(file_name):
    flying = 'flying'
    all_files = os.listdir(flying)
    return get_file_by_name(all_files, file_name)
    

arqs = load_file_in_flying('erros')
names_prod = ['gas', 'oil']

def plot_prod(arqs, name_prod='', plot_name=''):
    
    names_by_prod = get_file_by_name(arqs, name_prod)
    plt.clf()
    
    for case in list(descriptions.keys()):
        name_by_case = get_file_by_name(names_by_prod, case)[0]
        file_name = os.path.join('flying', name_by_case)
        arrays = np.load(file_name)
        log10 = arrays['log10_erro_rel']
        erro_rel = arrays['erro_rel']
        erro_abs = arrays['erro_abs']
        x = arrays['x']
        plt.plot(x, log10, label=descriptions[case])
        
    
    plt.ylabel('Log10 do erro relativo')
    plt.xlabel('time [days]')
    plt.legend()
    plt.savefig(plot_name)
    plt.clf()

plot_prod(arqs, name_prod='oil', plot_name='LOG10_PROBLEMA1.png')