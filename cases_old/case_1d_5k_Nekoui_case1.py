import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)

L = 1
for  arq in arquivos:
    if  arq.startswith(name):
        '''------------------------- RESULTS --------------------------'''



        '------------------------------- FOU ----------------------------------'

        datas = np.load('flying/results_5k_Nekoui_case2_200x1x1_IMPEC_FOU_11066.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_200_FOU = data[10]
            n=200
            x_200 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_200_FOU = data[7]
            t_200_FOU = data[2]

        plt.figure(1)
        plt.plot(x_200, Sg_200_FOU, 'r')
        plt.grid()
        plt.ylabel('Gas Saturation')
        plt.xlabel('Cells')
        plt.savefig('results/compositional/SPE_paper/5k_Sg_Nekoui_case2_200.png')

        plt.figure(2)
        plt.plot(x_200, z_200_FOU[0,:], 'r')
        plt.grid()
        plt.ylabel('C1')
        plt.xlabel('Cells')
        plt.savefig('results/compositional/SPE_paper/5k_N2_Nekoui_case2_200.png')

        import pdb; pdb.set_trace()
