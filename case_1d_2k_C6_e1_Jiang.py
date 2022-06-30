import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
L = 400
for  arq in arquivos:
    if  arq.startswith(name):
        '''------------------------- RESULTS------------------------'''

        datas = np.load('flying/results_case1_Jiang_2k_400_FOU_491.npy', allow_pickle=True)
        datas = np.load('flying/results_case3_Jiang_3k_400_FOU_1704.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FOU_400 = data[7]
            z_FOU_400 = data[10]
            x400 = np.linspace(0+L/800,L*(1-1/800),400)


        plt.figure(1)
        plt.plot(x400, Sg_FOU_400, 'r')
        plt.grid()
        plt.ylabel('Gas Saturation')
        plt.xlabel('Cells')
        plt.savefig('results/compositional/2k_C6_e1_Jiang_Sg.png')

        plt.figure(2)
        plt.plot(x400, z_FOU_400[0,:], 'r')
        plt.grid()
        plt.ylabel('C1')
        plt.xlabel('Cells')
        plt.savefig('results/compositional/2k_C6_e1_Jiang_zC1.png')

        import pdb; pdb.set_trace()