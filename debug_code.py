import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy import interpolate

flying = 'flying'
name = 'results_'
arquivos = os.listdir(flying)


for  arq in arquivos:
    if  arq.startswith(name):
        datas = np.load('flying/results_Firoo_2k_ex1_IMPEC_FOU_2545.npy', allow_pickle=True)
        for data in datas[0:]:
            import pdb; pdb.set_trace()
            zC1 = data[10][0]
            vpi = data[1]
            print('vpi: ', vpi)
            
            
