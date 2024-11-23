import os
from packs import directories as direc
import numpy as np
data_loaded = direc.data_loaded
__all__ = []

load = data_loaded['load_biphasic']
loop_maximo = data_loaded['biphasic_data']['loop_maximo']
tempo_maximo = data_loaded['biphasic_data']['tempo_maximo']

if not load:
    resposta = input('\nVoce deseja perder os dados da simulacao s/n\n')
    if resposta == 's':
        pass
    else:
        import sys
        print('\nAltere o valor de load_biphasic para True\n')
        sys.exit(0)

np.save(direc.name_load, np.array([load]))
aqui = os.path.dirname(os.path.abspath(__file__))
run_biphasic = os.path.join(aqui, 'run_biphasic.py')
run_biphasic = 'python3 ' + run_biphasic


os.system(run_biphasic)

import pdb; pdb.set_trace()
load = True
np.save(direc.name_load, np.array([load]))

while verif:
    os.system(run_biphasic)
    if loop_maximo or tempo_maximo:
        hist = np.load(direc.name_hist)

        if loop_maximo:
            if hist[0] > loop_maximo:
                verif = False
        if tempo_maximo:
            if hist[5] > tempo_maximo:
                verif = False
    import pdb; pdb.set_trace()
