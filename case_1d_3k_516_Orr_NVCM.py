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
        '''------------------------- RESULTS --------------------------'''

        datas = np.load('flying/results_Orr_3k_C516_8000x1x1_IMPEC_FOU_10980.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_8000_FOU = data[10]
            n=8000
            x_8000 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        'FR3'
        datas = np.load('flying/results_Orr_3k_C516_16x1x1_IMPEC_FR3_LLF_75.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_16_FR3 = data[10]
            n=16
            x_16 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C516_32x1x1_IMPEC_FR3_LLF_166.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_32_FR3 = data[10]
            n=32
            x_32 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C516_64x1x1_IMPEC_FR3_LLF_355.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_64_FR3 = data[10]
            n=64
            x_64 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C516_128x1x1_IMPEC_FR3_LLF_738.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_128_FR3 = data[10]
            n=128
            x_128 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        'FR4'
        datas = np.load('flying/results_Orr_3k_C516_16x1x1_IMPEC_FR4_LLF_103.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_16_FR4 = data[10]

        datas = np.load('flying/results_Orr_3k_C516_32x1x1_IMPEC_FR4_LLF_219.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_32_FR4 = data[10]

        datas = np.load('flying/results_Orr_3k_C516_64x1x1_IMPEC_FR4_LLF_471.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_64_FR4 = data[10]

        datas = np.load('flying/results_Orr_3k_C516_128x1x1_IMPEC_FR4_LLF_966.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_128_FR4 = data[10]


        plt.figure(1)
        s=3
        #plt.plot(x_128, z_128_FOU[0,:], '-..b', mfc='none', markersize=s)
        #plt.plot(x_128, z_128_MUSCL[0,:], ':m', mfc='none', markersize=s)
        #plt.plot(x_128, z_128_FR2[0,:], '-g', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR3[0,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR4[0,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k')
        plt.plot(x_8000, z_8000_FOU[0,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'Reference'))
        plt.title('Results for 128x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{CO2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_CO2_orr_C516_NVCM_128.png')
        plt.close()

        plt.figure(2)
        s=3
        #plt.plot(x_64, z_64_FOU[0,:], '-..b', mfc='none', markersize=s)
        #plt.plot(x_64, z_64_MUSCL[0,:], ':m', mfc='none', markersize=s)
        #plt.plot(x_64, z_64_FR2[0,:], '-g', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR3[0,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR4[0,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k')
        plt.plot(x_8000, z_8000_FOU[0,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'Reference'))
        plt.title('Results for 64x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{CO2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_CO2_orr_C516_NVCM_64.png')
        plt.close()

        plt.figure(3)
        s=3
        #plt.plot(x_32, z_32_FOU[0,:], '-..b', mfc='none', markersize=s)
        #plt.plot(x_32, z_32_MUSCL[0,:], ':m', mfc='none', markersize=s)
        #plt.plot(x_32, z_32_FR2[0,:], '-g', mfc='none', markersize=s)
        plt.plot(x_32, z_32_FR3[0,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_32, z_32_FR4[0,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k')
        plt.plot(x_8000, z_8000_FOU[0,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'Reference'))
        plt.title('Results for 32x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{CO2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_CO2_orr_C516_NVCM_32.png')
        plt.close()

        plt.figure(1)
        s=3
        #plt.plot(x_128, z_128_FOU[1,:], '-..b', mfc='none', markersize=s)
        #plt.plot(x_128, z_128_MUSCL[1,:], ':m', mfc='none', markersize=s)
        #plt.plot(x_128, z_128_FR2[1,:], '-g', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR3[1,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR4[1,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k')
        plt.plot(x_8000, z_8000_FOU[1,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'Reference'))
        plt.title('Results for 128x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C1_orr_C516_NVCM_128.png')
        plt.close()

        plt.figure(2)
        s=3
        #plt.plot(x_64, z_64_FOU[1,:], '-..b', mfc='none', markersize=s)
        #plt.plot(x_64, z_64_MUSCL[1,:], ':m', mfc='none', markersize=s)
        #plt.plot(x_64, z_64_FR2[1,:], '-g', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR3[1,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR4[1,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k')
        plt.plot(x_8000, z_8000_FOU[1,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'Reference'))
        plt.title('Results for 64x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C1_orr_C516_NVCM_64.png')
        plt.close()

        plt.figure(3)
        s=3
        #plt.plot(x_32, z_32_FOU[1,:], '-..b', mfc='none', markersize=s)
        #plt.plot(x_32, z_32_MUSCL[1,:], ':m', mfc='none', markersize=s)
        #plt.plot(x_32, z_32_FR2[1,:], '-g', mfc='none', markersize=s)
        plt.plot(x_32, z_32_FR3[1,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_32, z_32_FR4[1,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k')
        plt.plot(x_8000, z_8000_FOU[1,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'Reference'))
        plt.title('Results for 32x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C1_orr_C516_NVCM_32.png')
        plt.close()

        plt.figure(1)
        s=3
        #plt.plot(x_128, z_128_FOU[2,:], '-..b', mfc='none', markersize=s)
        #plt.plot(x_128, z_128_MUSCL[2,:], ':m', mfc='none', markersize=s)
        #plt.plot(x_128, z_128_FR2[2,:], '-g', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR3[2,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR4[2,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k')
        plt.plot(x_8000, z_8000_FOU[2,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'Reference'))
        plt.title('Results for 128x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C10}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C10_orr_C516_NVCM_128.png')
        plt.close()

        plt.figure(2)
        s=3
        #plt.plot(x_64, z_64_FOU[2,:], '-..b', mfc='none', markersize=s)
        #plt.plot(x_64, z_64_MUSCL[2,:], ':m', mfc='none', markersize=s)
        #plt.plot(x_64, z_64_FR2[2,:], '-g', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR3[2,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR4[2,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k')
        plt.plot(x_8000, z_8000_FOU[2,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'Reference'))
        plt.title('Results for 64x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C10}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C10_orr_C516_NVCM_64.png')
        plt.close()

        plt.figure(3)
        s=3
        #plt.plot(x_32, z_32_FOU[2,:], '-..b', mfc='none', markersize=s)
        #plt.plot(x_32, z_32_MUSCL[2,:], ':m', mfc='none', markersize=s)
        #plt.plot(x_32, z_32_FR2[2,:], '-g', mfc='none', markersize=s)
        plt.plot(x_32, z_32_FR3[2,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_32, z_32_FR4[2,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k')
        plt.plot(x_8000, z_8000_FOU[2,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'Reference'))
        plt.title('Results for 32x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C10}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C10_orr_C516_NVCM_32.png')
        plt.close()

        import pdb; pdb.set_trace()
