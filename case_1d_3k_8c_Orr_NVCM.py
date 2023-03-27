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

        x_zCH4 = np.array([0, 0.254, 0.254, 0.769, 0.769, 1.288, 1.288, 2])
        zCH4 = np.array([1, 1, 0.85, 0.847, 0.617,0.576, 0, 0])
        x_zCO2 = np.array([0, 0.254, 0.772, 0.772, 1.288, 1.288, 2])
        zCO2 = np.array([0, 0, 0, 0.203, 0.216, 0.377, 0.377])

        datas = np.load('flying/results_Orr_3k_C8_8000x1x1_IMPEC_FOU_12266.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_8000_FOU = data[10]
            n=8000
            x_8000 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)


        'FR2'

        datas = np.load('flying/results_Orr_3k_C8_128x1x1_IMPEC_FR2_LLF_481.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_128_FR2 = data[10]

        datas = np.load('flying/results_Orr_3k_C8_64x1x1_IMPEC_FR2_LLF_234.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_64_FR2 = data[10]

        'FR3'

        datas = np.load('flying/results_Orr_3k_C8_128x1x1_IMPEC_FR3_LLF_724.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_128_FR3 = data[10]
            n=128
            x_128 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C8_64x1x1_IMPEC_FR3_LLF_346.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_64_FR3 = data[10]
            n=64
            x_64 = np.linspace(0+2/(2*n),2*(1-1/(2*n)),n)

        'FR4'

        datas = np.load('flying/results_Orr_3k_C8_128x1x1_IMPEC_FR4_LLF_951.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_128_FR4 = data[10]

        datas = np.load('flying/results_Orr_3k_C8_64x1x1_IMPEC_FR4_LLF_454.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_64_FR4 = data[10]

        'MUSCL'

        datas = np.load('flying/results_Orr_3k_C8_128x1x1_IMPEC_MUSCL_LLF_492.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_128_MUSCL = data[10]

        datas = np.load('flying/results_Orr_3k_C8_64x1x1_IMPEC_MUSCL_LLF_240.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_64_MUSCL = data[10]

        'FOU'

        datas = np.load('flying/results_Orr_3k_C8_128x1x1_IMPEC_FOU_132.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_128_FOU = data[10]

        datas = np.load('flying/results_Orr_3k_C8_64x1x1_IMPEC_FOU_67.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_64_FOU = data[10]


        plt.figure(1)
        s=3
        plt.plot(x_64, z_64_FOU[0,:], '-..b', mfc='none', markersize=s)
        plt.plot(x_64, z_64_MUSCL[0,:], ':m', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR2[0,:], '-g', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR3[0,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR4[0,:], '--r', mfc='none', markersize=s)
        plt.plot(x_zCO2, zCO2, '-k',linewidth=2)
        #plt.plot(x_8000, z_8000_FOU[0,:], '--k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'MOC'))
        plt.title('Results for 64x1x1 mesh')
        plt.grid()

        #plt.xlim(1.1,1.35)
        plt.ylabel('$z_{CO2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_CO2_orr_C8_NVCM_64.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_64, z_64_FOU[0,:], '-..b', mfc='none', markersize=s)
        plt.plot(x_64, z_64_MUSCL[0,:], ':m', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR2[0,:], '-g', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR3[0,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR4[0,:], '--r', mfc='none', markersize=s)
        plt.plot(x_zCO2, zCO2, '-k',linewidth=2)
        #plt.plot(x_8000, z_8000_FOU[0,:], '--k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'MOC'))
        plt.title('Results for 64x1x1 mesh')
        plt.grid()
        plt.xlim(0.7,1.5)
        plt.ylabel('$z_{CO2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_CO2_orr_C8_NVCM_64_zoom.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_64, z_64_FOU[1,:], '-..b', mfc='none', markersize=s)
        plt.plot(x_64, z_64_MUSCL[1,:], ':m', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR2[1,:], '-g', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR3[1,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR4[1,:], '--r', mfc='none', markersize=s)
        plt.plot(x_zCH4, zCH4, '-k',linewidth=2)
        #plt.plot(x_8000, z_8000_FOU[1,:], '--k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'MOC'))
        plt.title('Results for 64x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C1_orr_C8_NVCM_64.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_64, z_64_FOU[1,:], '-..b', mfc='none', markersize=s)
        plt.plot(x_64, z_64_MUSCL[1,:], ':m', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR2[1,:], '-g', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR3[1,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR4[1,:], '--r', mfc='none', markersize=s)
        plt.plot(x_zCH4, zCH4, '-k', linewidth=2)
        plt.xlim(0.1,1.4)
        plt.ylim(0.4,1.01)
        #plt.plot(x_8000, z_8000_FOU[1,:], '--k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'MOC'))
        plt.title('Results for 64x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C1_orr_C8_NVCM_64_zoom.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_64, z_64_FOU[2,:], '-..b', mfc='none', markersize=s)
        plt.plot(x_64, z_64_MUSCL[2,:], ':m', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR2[2,:], '-g', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR3[2,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_64, z_64_FR4[2,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k',linewidth=2)
        #plt.plot(x_8000, z_8000_FOU[0,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'MOC'))
        plt.title('Results for 64x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C10}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C10_orr_C8_NVCM_64.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_128, z_128_FOU[0,:], '-..b', mfc='none', markersize=s)
        plt.plot(x_128, z_128_MUSCL[0,:], ':m', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR2[0,:], '-g', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR3[0,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR4[0,:], '--r', mfc='none', markersize=s)
        plt.plot(x_zCO2, zCO2, '-k',linewidth=2)
        #plt.plot(x_8000, z_8000_FOU[0,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'MOC'))
        plt.title('Results for 128x1x1 mesh')
        plt.grid()
        #plt.xlim(1.1,1.35)
        plt.ylabel('$z_{CO2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_CO2_orr_C8_NVCM_128.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_128, z_128_FOU[0,:], '-..b', mfc='none', markersize=s)
        plt.plot(x_128, z_128_MUSCL[0,:], ':m', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR2[0,:], '-g', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR3[0,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR4[0,:], '--r', mfc='none', markersize=s)
        plt.plot(x_zCO2, zCO2, '-k',linewidth=2)
        #plt.plot(x_8000, z_8000_FOU[0,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'MOC'))
        plt.title('Results for 128x1x1 mesh')
        plt.grid()
        plt.xlim(0.7,1.5)
        plt.ylabel('$z_{CO2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_CO2_orr_C8_NVCM_128_zoom.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_128, z_128_FOU[1,:], '-..b', mfc='none', markersize=s)
        plt.plot(x_128, z_128_MUSCL[1,:], ':m', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR2[1,:], '-g', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR3[1,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR4[1,:], '--r', mfc='none', markersize=s)
        plt.plot(x_zCH4, zCH4, '-k',linewidth=2)
        #plt.plot(x_8000, z_8000_FOU[0,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'MOC'))
        plt.title('Results for 128x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C1_orr_C8_NVCM_128.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_128, z_128_FOU[1,:], '-..b', mfc='none', markersize=s)
        plt.plot(x_128, z_128_MUSCL[1,:], ':m', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR2[1,:], '-g', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR3[1,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR4[1,:], '--r', mfc='none', markersize=s)
        plt.plot(x_zCH4, zCH4, '-k', linewidth=2)
        plt.xlim(0.1,1.4)
        plt.ylim(0.4,1.01)
        #plt.plot(x_8000, z_8000_FOU[0,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'MOC'))
        plt.title('Results for 128x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C1_orr_C8_NVCM_128_zoom.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_128, z_128_FOU[2,:], '-..b', mfc='none', markersize=s)
        plt.plot(x_128, z_128_MUSCL[2,:], ':m', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR2[2,:], '-g', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR3[2,:], '-.y', mfc='none', markersize=s)
        plt.plot(x_128, z_128_FR4[2,:], '--r', mfc='none', markersize=s)
        #plt.plot(zC10_x, zC10, '-k',linewidth=2)
        plt.plot(x_8000, z_8000_FOU[2,:], '-k', linewidth=2)
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'MOC'))
        plt.title('Results for 128x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C10}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_C10_orr_C8_NVCM_128.png')
        plt.close()

        import pdb; pdb.set_trace()
