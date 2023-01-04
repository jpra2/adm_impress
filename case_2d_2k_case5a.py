import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
n = 40

for arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_Hoteit_Firoo_2k_ex5a_40x20_FI_teste_VPI0_4_dtmin_864_tempo_1025.npy', allow_pickle=True)
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

        datas2 = np.load('flying/results_Hoteit_Firoo_2k_ex5a_80x40_FI_VPI0_45_dtmin_864_1331.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data_2 in datas2[1:]:
            Sw_FI_e4 = data_2[5]
            So_FI_e4 = data_2[6]
            Sg_FI_e4 = data_2[7]
            Oil_p_FI_e4 = data_2[8]
            Gas_p_FI_e4 = data_2[9]
            pressure_FI_e4 = data_2[4]/1e3
            time_FI_e4 = data_2[3]
            zC1_FI_e4 = data_2[10][0]
            x1 = np.linspace(0, 170.69, n)

        import pdb; pdb.set_trace()
