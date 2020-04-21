import os
import yaml
import numpy as np


netas_limite=np.array([1000000000.0, 10.0, 5.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,0.05,0.025,0.0])
for neta in netas_limite:
    np.save('flying/neta_lim.npy',np.array([neta]))
    os.system('python ADM.py')
