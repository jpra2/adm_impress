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
        L=2
        '''-------------------------MUSCL LLF RESULTS------------------------'''

        zCH4_x = np.array([0, 0.243506, 0.243506, 0.296898, 0.384921, 0.516955, \
                           0.631313, 0.645382, 0.655123, 0.665224, 0.677489, 0.689755, \
                           0.702742, 0.720779, 0.73557, 0.755411, 0.774531, \
                           0.790765, 0.814214, 0.840548, 1.49856, 1.49928, \
                           1.99956])*1
        zCH4 = np.array([0.799874, 0.799874, 0.652146, 0.640152, 0.624369, \
                         0.604798, 0.589015, 0.570076, 0.561869, 0.558712, 0.557449, \
                         0.55556, 0.554924, 0.556818, 0.560606, 0.564394, \
                         0.569444, 0.575126, 0.587121, 0.600379, 0.600379, \
                         0.300505, 0.300505])

        datas = np.load('flying/results_Orr_3k_C518_100x1x1_IMPEC_FR2_402.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_100_FR2 = data[10][1]
            n=100
            x_100 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C518_100x1x1_IMPEC_FOU_127.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_100_FOU = data[10][1]

        datas = np.load('flying/results_Orr_3k_C518_100x1x1_IMPEC_FR3_657.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_100_FR3 = data[10][1]

        datas = np.load('flying/results_Orr_3k_C518_100x1x1_IMPEC_FR4_910.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_100_FR4 = data[10][1]
            t100_FR4 = data[2]


        datas = np.load('flying/results_Orr_3k_C518_100x1x1_IMPEC_FR4_RK3_374.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_100_FR4_RK3 = data[10][1] #CFL2
            t100_FR4_RK3_CFL2 = data[2]

        datas = np.load('flying/results_Orr_3k_C518_200x1x1_IMPEC_FR2_932.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_200_FR2 = data[10][0] #CFL2


        datas = np.load('flying/results_Orr_3k_C518_200x1x1_IMPEC_FOU_253.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_200_FOU = data[10][1]
            n=200
            x_200 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)

        plt.figure(1)
        plt.plot(x_100, zCH4_100_FR2, '-bo', mfc='none')
        plt.plot(x_100, zCH4_100_FR3, '-gs', mfc='none')
        plt.plot(x_100, zCH4_100_FR4, '-mv', mfc='none')
        plt.plot(x_100, zCH4_100_FOU, '-r*', mfc='none')
        #plt.plot(x_200, zCH4_200_FOU, '-r*', mfc='none')
        plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('FR P2 - 100 CV','FR P3 - 100 CV', 'FR P4 - 100 CV', \
         'FOU - 100 CV' ,'MOC (Orr, 2005)'))
        plt.ylabel('Methane global molar fraction ')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_LVI_methane_Orr_NVCM_100.png')

        plt.figure(2)
        plt.plot(x_100, zCH4_100_FR4_RK3, '-mv', mfc='none')
        plt.plot(x_200, zCH4_200_FR2, '-co', mfc='none')
        plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        plt.legend(('FR P4 - 100 CV','FR P2 RK3 - 200 CV','MOC (Orr, 2005)'))
        plt.ylabel('Methane global molar fraction ')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_LVI_methane_Orr_NVCM_FR4_100.png')

        import pdb; pdb.set_trace()
