import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
x_CMG = np.loadtxt('xCH4_case2MM_CMG.txt')
xCH4_CMG = np.loadtxt('xxCH4_case2MM_CMG.txt')

for  arq in arquivos:
    if  arq.startswith(name):
        '''-------------------------MUSCL LLF RESULTS------------------------'''

        zCH4_x = np.array([0, 0.241998, 0.24356, 0.24356, 0.529274, 0.67031, \
                           0.959147, 1.13089, 1.13297, 1.48374, 1.48374, 1.99948])
        zCH4 = np.array([0.602529, 0.600723, 0.523035, 0.490515, 0.456188, \
                         0.439928, 0.414634, 0.411924, 0.592593, 0.592593, \
                         0.303523, 0.301716])



        datas = np.load('flying/results_Orr_2k_C4_IMPEC_FOU_1475.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCO2_200_FOU = data[10][0]
            n=200
            x_200 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)


        plt.figure(1)
        plt.plot(x_200, zCO2_200_FOU, '-bo', mfc='none')

        #plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        plt.ylabel('CO$_2$ global molar fraction ')
        plt.title('HIgh volatile intermediate component case - Orr')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/2k_CO2_Orr_NVCM.png')
        import pdb; pdb.set_trace()
