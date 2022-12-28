import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
n = 20

for arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FI_e4 = data[5]
            So_FI_e4 = data[6]
            Sg_FI_e4 = data[7]
            Oil_p_FI_e4 = data[8]
            Gas_p_FI_e4 = data[9]
            pressure_FI_e4 = data[4]/1e3
            time_FI_e4 = data[3]
            zC1_FI_e4 = data[10][0]
            x1 = np.linspace(0, 170.69, n)

        import pdb; pdb.set_trace()
