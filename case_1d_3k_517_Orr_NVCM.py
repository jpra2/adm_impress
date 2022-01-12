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

        zCH4_x = np.array([0, 0.24356, 0.24356, 0.39084, 0.531356, \
                           0.710383, 0.81707, 0.934166, 1.06583, 1.13193, \
                           1.13193, 1.4827, 1.48322, 1.99532])
        zCH4 = np.array([0.600723, 0.600723, 0.48785, 0.466125, 0.450768, \
                         0.435411, 0.427281, 0.419151, 0.411021, 0.409214, \
                         0.589883, 0.589883, 0.300813, 0.300813])



        datas = np.load('flying/results_Orr_3k_C5_IMPEC_FOU_617.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_200_FOU = data[10][1]
            n=200
            x_200 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C5_100x1x1_IMPEC_FR2_398.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_100_FR2 = data[10][1]
            n=100
            x_100 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C5_100x1x1_IMPEC_FR3_653.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_100_FR3 = data[10][1]

        datas = np.load('flying/results_Orr_3k_C5_IMPEC_FR3_1337.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_200_FR3 = data[10][1]

        datas = np.load('flying/results_Orr_3k_C5_100x1x1_IMPEC_FR4_937.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_100_FR4 = data[10][1]

        '''Convergence test'''
        'FR2'
        datas = np.load('flying/results_Orr_3k_C5_8x1x1_IMPEC_FR2_31.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_8_FR2 = data[10][1]
            n=8
            x_8 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C5_16x1x1_IMPEC_FR2_63.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_16_FR2 = data[10][1]
            n=16
            x_16 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C5_32x1x1_IMPEC_FR2_127.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_32_FR2 = data[10][1]
            n=32
            x_32 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C5_64x1x1_IMPEC_FR2_257.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_64_FR2 = data[10][1]
            n=64
            x_64 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C5_64x1x1_IMPEC_FR2_753.npy', allow_pickle=True) #CFL 0,3
        for data in datas[datas.shape[0]-1:]:
            zCH4_64_FR2_2 = data[10][1]
            n=64
            x_64 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C5_128x1x1_IMPEC_FR2_518.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_128_FR2 = data[10][1]
            n=128
            x_128 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C5_256x1x1_IMPEC_FR2_1890.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCH4_256_FR2 = data[10][1]
            n=256
            x_256 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)


        plt.figure(1)
        plt.plot(x_8, zCH4_8_FR2, '-g', mfc='none')
        plt.plot(x_16, zCH4_16_FR2, '-r', mfc='none')
        plt.plot(x_32, zCH4_32_FR2, '-m', mfc='none')
        plt.plot(x_64, zCH4_64_FR2, '-b', mfc='none')
        plt.plot(x_128, zCH4_128_FR2, '-y', mfc='none')
        plt.plot(x_256, zCH4_256_FR2, '-c', mfc='none')
        plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        plt.legend(('8 CV', '16 CV', '32 CV', '64 CV', '128 CV', '256 CV', 'MOC (Orr, 2005)'))
        plt.ylabel('Methane global molar fraction ')
        plt.title('FR P1')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_methane_Orr_NVCM_FR2.png')

        plt.figure(2)
        plt.plot(x_200, zCH4_200_FOU, '-bo')
        plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        
        plt.legend(('FOU 200 CV', 'MOC (Orr, 2005)'))
        plt.ylabel('Methane global molar fraction ')
        plt.title('FR P1 64 CV')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_methane_Orr_NVCM_FOU.png')
        import pdb; pdb.set_trace()
