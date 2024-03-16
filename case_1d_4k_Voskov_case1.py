import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)

L = 500
for  arq in arquivos:
    if  arq.startswith(name):
        '''------------------------- RESULTS --------------------------'''



        '------------------------------- FOU ----------------------------------'

        datas = np.load('flying/results_4k_Voskov_case1_adapt_200x1x1_IMPEC_FOU_776.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_200_FOU = data[10]
            n=200
            x_200 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_200_FOU = data[7]
            t_200_FOU = data[2]

        datas = np.load('flying/results_4k_Voskov_case1_500x1x1_IMPEC_FOU_1396.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_500_FOU = data[10]
            n=500
            x_500 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_500_FOU = data[7]
            t_500_FOU = data[2]

        datas = np.load('flying/results_4k_Voskov_case1_1000x1x1_IMPEC_FOU_2758.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_1000_FOU = data[10]
            n=1000
            x_1000 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_1000_FOU = data[7]
            t_1000_FOU = data[2]

        plt.figure(1)
        plt.plot(x_200, Sg_200_FOU, 'r')
        plt.plot(x_500, Sg_500_FOU, 'g')
        plt.plot(x_1000, Sg_1000_FOU, 'k')
        plt.grid()
        plt.ylabel('Gas Saturation')
        plt.xlabel('Cells')
        plt.savefig('results/compositional/SPE_paper/5k_Sg_Voskov_case1_500.png')

        plt.figure(2)
        plt.plot(x_200, z_200_FOU[0,:], 'r')
        plt.grid()
        plt.ylabel('C1')
        plt.xlabel('Cells')
        plt.savefig('results/compositional/SPE_paper/5k_C1_Voskov_case1_500.png')

        import pdb; pdb.set_trace()
