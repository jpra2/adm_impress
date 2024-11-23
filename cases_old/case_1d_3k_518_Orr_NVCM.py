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

        datas = np.load('flying/results_Orr_3k_C518_50x1x1_IMPEC_FR2_RK3_199.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_50_FR2 = data[10]
            Sg_50_FR2 = data[7]
            t50_FR2 = data[2]

        datas = np.load('flying/results_Orr_3k_C518_50x1x1_IMPEC_FR4_RK3_422.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_50_FR4 = data[10]
            Sg_50_FR4 = data[7]
            t50_FR4 = data[2]
            n=50
            x_50 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C518_100x1x1_IMPEC_FR4_RK3_864.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_100_FR4 = data[10]
            Sg_100_FR4 = data[7]
            t100_FR4 = data[2]
            n=100
            x_100 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)

        'MUSCL'
        datas = np.load('flying/results_Orr_3k_C518_50x1x1_IMPEC_MUSCL_LLF_196.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_50_MUSCL = data[10]
            Sg_50_MUSCL = data[7]
            t50_MUSCL = data[2]

        datas = np.load('flying/results_Orr_3k_C518_100x1x1_IMPEC_MUSCL_LLF_398.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_100_MUSCL = data[10]
            Sg_100_MUSCL = data[7]
            t100_MUSCL = data[2]
            n=100
            x_100 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C518_200x1x1_IMPEC_MUSCL_LLF_801.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_200_MUSCL = data[10]
            Sg_200_MUSCL = data[7]
            t200_MUSCL = data[2]
            n=200
            x_200 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C518_300x1x1_IMPEC_MUSCL_LLF_1226.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_300_MUSCL = data[10]
            Sg_300_MUSCL = data[7]
            t300_MUSCL = data[2]
            n=300
            x_300 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C518_400x1x1_IMPEC_MUSCL_LLF_1746.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_400_MUSCL = data[10]
            Sg_400_MUSCL = data[7]
            t400_MUSCL = data[2]
            n=400
            x_400 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)


        'FOU'
        datas = np.load('flying/results_Orr_3k_C518_50x1x1_IMPEC_FOU_63.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_50_FOU = data[10]
            Sg_50_FOU = data[7]
            t50_FOU = data[2]

        datas = np.load('flying/results_Orr_3k_C518_100x1x1_IMPEC_FOU_127.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_100_FOU = data[10]
            Sg_100_FOU = data[7]
            t100_FOU = data[2]

        datas = np.load('flying/results_Orr_3k_C518_500x1x1_IMPEC_FOU_638.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_500_FOU = data[10]
            Sg_500_FOU = data[7]
            t500_FOU = data[2]
            n=500
            x_500 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)

        datas = np.load('flying/results_Orr_3k_C518_8000x1x1_IMPEC_FOU_12404.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_8000_FOU = data[10]
            Sg_8000_FOU = data[7]
            t8000_FOU = data[2]
            n=8000
            x_8000 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)



        plt.figure(1)
        plt.plot(x_8000, z_8000_FOU[0,:], '-k')
        plt.plot(x_50, z_50_FOU[0,:], '--b', mfc='none')
        plt.plot(x_50, z_50_FR4[0,:], '-r*', mfc='none')
        #plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (50)', 'FR-P3 (50)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C1}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_methane_Orr_NVCM_50.png')

        plt.figure(2)
        plt.plot(x_8000, z_8000_FOU[0,:], '-k')
        plt.plot(x_100, z_100_FOU[0,:], '--b', mfc='none')
        plt.plot(x_100, z_100_FR4[0,:], '-r*', mfc='none')
        #plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (100)', 'FR-P3 (100)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C1}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_methane_Orr_NVCM_100.png')

        plt.figure(3)
        plt.plot(x_8000, z_8000_FOU[0,:], '-k')
        plt.plot(x_500, z_500_FOU[0,:], '--b', mfc='none')
        plt.plot(x_100, z_100_FR4[0,:], '-r*', mfc='none')
        #plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (500)','FR-P3 (100)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C1}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_methane_Orr_NVCM_500.png')
        plt.close()

        plt.figure(4)
        plt.plot(x_8000, z_8000_FOU[1,:], '-k')
        plt.plot(x_50, z_50_FOU[1,:], '--b', mfc='none')
        plt.plot(x_50, z_50_FR4[1,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (50)', 'FR-P3 (50)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C4}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C4_Orr_NVCM_50.png')
        plt.close()

        plt.figure(5)
        plt.plot(x_8000, z_8000_FOU[1,:], '-k')
        plt.plot(x_100, z_100_FOU[1,:], '--b', mfc='none')
        plt.plot(x_100, z_100_FR4[1,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (100)', 'FR-P3 (100)'))
        plt.ylabel('$z_{C4}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C4_Orr_NVCM_100.png')
        plt.close()

        plt.figure(6)
        plt.plot(x_8000, z_8000_FOU[1,:], '-k')
        plt.plot(x_500, z_500_FOU[1,:], '--b', mfc='none')
        plt.plot(x_100, z_100_FR4[1,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (500)', 'FR-P3 (100)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C4}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C4_Orr_NVCM_500.png')
        plt.close()

        plt.figure(7)
        plt.plot(x_8000, z_8000_FOU[2,:], '-k')
        plt.plot(x_50, z_50_FOU[2,:], '--b', mfc='none')
        plt.plot(x_50, z_50_FR4[2,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference','FOU (50)', 'FR-P3 (50)'))# 'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C10}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C10_Orr_NVCM_50.png')
        plt.close()

        plt.figure(8)
        plt.plot(x_8000, z_8000_FOU[2,:], '-k')
        plt.plot(x_100, z_100_FOU[2,:], '--b', mfc='none')
        plt.plot(x_100, z_100_FR4[2,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (100)', 'FR-P3 (100)'))
        plt.ylabel('$z_{C10}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C10_Orr_NVCM_100.png')
        plt.close()

        plt.figure(9)
        plt.plot(x_8000, z_8000_FOU[2,:], '-k')
        plt.plot(x_500, z_500_FOU[2,:], '--b', mfc='none')
        plt.plot(x_100, z_100_FR4[2,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (500)','FR-P3 (100)'))# 'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C10}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C10_Orr_NVCM_500.png')

        plt.figure(10)
        plt.plot(x_8000, Sg_8000_FOU, '-k')
        plt.plot(x_50, Sg_50_FOU, '--b', mfc='none')
        plt.plot(x_50, Sg_50_FR4, '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (50)', 'FR-P3 (50)'))
        plt.ylabel('$Sg$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_Sg_Orr_NVCM_50.png')

        plt.figure(11)
        plt.plot(x_8000, Sg_8000_FOU, '-k')
        plt.plot(x_100, Sg_100_FOU, '--b', mfc='none')
        plt.plot(x_100, Sg_100_FR4, '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (100)', 'FR-P3 (100)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$Sg$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_Sg_Orr_NVCM_100.png')

        plt.figure(12)
        plt.plot(x_8000, Sg_8000_FOU, '-k')
        plt.plot(x_500, Sg_500_FOU, '--b', mfc='none')
        plt.plot(x_100, Sg_100_FR4, '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (500)', 'FR-P3 (100)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$Sg$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_Sg_Orr_NVCM_500.png')

        plt.figure(13)
        plt.plot(x_8000, z_8000_FOU[0,:], '-k')
        plt.plot(x_50, z_50_FOU[0,:], '--b', mfc='none')
        plt.plot(x_50, z_50_MUSCL[0,:], '-.m', mfc='none')
        plt.plot(x_50, z_50_FR4[0,:], '-r*', mfc='none')
        #plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (50)', 'MUSCL (50)', 'FR-P3 (50)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C1}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_methane_Orr_NVCM_50_and_MUSCL.png')

        plt.figure(14)
        plt.plot(x_8000, z_8000_FOU[0,:], '-k')
        plt.plot(x_100, z_100_FOU[0,:], '--b', mfc='none')
        plt.plot(x_100, z_100_MUSCL[0,:], '-.m', mfc='none')
        plt.plot(x_100, z_100_FR4[0,:], '-r*', mfc='none')
        #plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (100)', 'MUSCL (100)', 'FR-P3 (100)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C1}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_methane_Orr_NVCM_100_and_MUSCL.png')

        plt.figure(15)
        plt.plot(x_8000, z_8000_FOU[0,:], '-k')
        plt.plot(x_500, z_500_FOU[0,:], '--b', mfc='none')
        plt.plot(x_100, z_100_MUSCL[0,:], '-.m', mfc='none')
        plt.plot(x_100, z_100_FR4[0,:], '-r*', mfc='none')
        #plt.plot(zCH4_x, zCH4, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (500)', 'MUSCL (100)','FR-P3 (100)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C1}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_methane_Orr_NVCM_500_and_MUSCL.png')

        plt.figure(16)
        plt.plot(x_8000, z_8000_FOU[1,:], '-k')
        plt.plot(x_50, z_50_FOU[1,:], '--b', mfc='none')
        plt.plot(x_50, z_50_MUSCL[1,:], '-.m', mfc='none')
        plt.plot(x_50, z_50_FR4[1,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (50)', 'MUSCL (50)', 'FR-P3 (50)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C4}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C4_Orr_NVCM_50_and_MUSCL.png')

        plt.figure(17)
        plt.plot(x_8000, z_8000_FOU[1,:], '-k')
        plt.plot(x_100, z_100_FOU[1,:], '--b', mfc='none')
        plt.plot(x_100, z_100_MUSCL[1,:], '-.m', mfc='none')
        plt.plot(x_100, z_100_FR4[1,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (100)', 'MUSCL (100)', 'FR-P3 (100)'))
        plt.ylabel('$z_{C4}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C4_Orr_NVCM_100_and_MUSCL.png')

        plt.figure(18)
        plt.plot(x_8000, z_8000_FOU[1,:], '-k')
        plt.plot(x_500, z_500_FOU[1,:], '--b', mfc='none')
        plt.plot(x_100, z_100_MUSCL[1,:], '-.m', mfc='none')
        plt.plot(x_100, z_100_FR4[1,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (500)', 'MUSCL (100)', 'FR-P3 (100)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C4}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C4_Orr_NVCM_500_and_MUSCL.png')

        plt.figure(19)
        plt.plot(x_8000, z_8000_FOU[2,:], '-k')
        plt.plot(x_50, z_50_FOU[2,:], '--b', mfc='none')
        plt.plot(x_50, z_50_MUSCL[2,:], '-.m', mfc='none')
        plt.plot(x_50, z_50_FR4[2,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference','FOU (50)', 'MUSCL (50)', 'FR-P3 (50)'))# 'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C10}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C10_Orr_NVCM_50_and_MUSCL.png')

        plt.figure(20)
        plt.plot(x_8000, z_8000_FOU[2,:], '-k')
        plt.plot(x_100, z_100_FOU[2,:], '--b', mfc='none')
        plt.plot(x_100, z_100_MUSCL[2,:], '-.m', mfc='none')
        plt.plot(x_100, z_100_FR4[2,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (100)', 'MUSCL (100)', 'FR-P3 (100)'))
        plt.ylabel('$z_{C10}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C10_Orr_NVCM_100_and_MUSCL.png')

        plt.figure(21)
        plt.plot(x_8000, z_8000_FOU[2,:], '-k')
        plt.plot(x_500, z_500_FOU[2,:], '--b', mfc='none')
        plt.plot(x_100, z_100_MUSCL[2,:], '-.m', mfc='none')
        plt.plot(x_100, z_100_FR4[2,:], '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (500)', 'MUSCL (100)','FR-P3 (100)'))# 'MOC (Orr, 2005)'))
        plt.ylabel('$z_{C10}$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_C10_Orr_NVCM_500_and_MUSCL.png')

        plt.figure(22)
        plt.plot(x_8000, Sg_8000_FOU, '-k')
        plt.plot(x_50, Sg_50_FOU, '--b', mfc='none')
        #plt.plot(x_50, Sg_50_FR2, '-g', mfc='none')
        plt.plot(x_50, Sg_50_MUSCL, '-.m', mfc='none')
        plt.plot(x_50, Sg_50_FR4, '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (50)', 'MUSCL (50)', 'FR-P3 (50)'))
        plt.ylabel('$Sg$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_Sg_Orr_NVCM_50_and_MUSCL.png')

        plt.figure(23)
        plt.plot(x_8000, Sg_8000_FOU, '-k')
        plt.plot(x_100, Sg_100_FOU, '--b', mfc='none')
        plt.plot(x_100, Sg_100_MUSCL, '-.m', mfc='none')
        plt.plot(x_100, Sg_100_FR4, '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (100)', 'MUSCL (100)', 'FR-P3 (100)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$Sg$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_Sg_Orr_NVCM_100_and_MUSCL.png')

        plt.figure(24)
        plt.plot(x_8000, Sg_8000_FOU, '-k')
        plt.plot(x_500, Sg_500_FOU, '--b', mfc='none')
        plt.plot(x_100, Sg_100_MUSCL, '-.m', mfc='none')
        plt.plot(x_100, Sg_100_FR4, '-r*', mfc='none')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        #plt.ylim(0,1)
        plt.xlim(0,L)
        plt.legend(('Reference', 'FOU (500)', 'MUSCL (100)', 'FR-P3 (100)'))#'MOC (Orr, 2005)'))
        plt.ylabel('$Sg$')
        plt.title('Low Volatile Intermediate Component Case')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/3k_LVI_Sg_Orr_NVCM_500_and_MUSCL.png')


        import pdb; pdb.set_trace()
