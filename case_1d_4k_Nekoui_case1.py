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

        '''x_zC1 = np.array([0, 0.162, 0.162, 0.223, 0.237, 0.736, 0.745, 0.755, \
                          0.763, 0.77, 0.78, 0.79, 0.796, 0.86, 0.86, 1])*L
        zC1 = np.array([0.90, 0.90, 0.736, 0.71,0.71, 0.71, 0.63, 0.55, 0.48, \
                        0.44, 0.38, 0.32, 0.28, 0.28, 0., 0])
        x_zC2 = np.array([0, 0.16, 0.163, 0.168, 0.736, 0.74, 0.7407, 0.7433, \
                          0.7477, 0.751, 0.755, 0.76, 0.764, 0.766, 0.77, \
                          0.774, 0.778, 0.782, 0.786, 0.794, 0.833, 0.86, 0.86, 1])*L
        zC2 = np.array([0.1, 0.1, 0.092, 0.091, 0.09, 0.1, 0.13, 0.16, \
                        0.193, 0.225, 0.26, 0.282, 0.305, 0.33, 0.35, 0.37, 0.39, \
                        0.415, 0.435, 0.469, 0.469, 0.467, 0.25, 0.25])'''
        x_ans_03t = np.array([0, 0.0288, 0.0288, 0.0968, 0.0968, 0.1522, 0.1596, \
                              0.1616, 0.207, 0.207, 1])

        Sg_ans_03t = np.array([1, 1, 0.775, 0.775, 0.571, 0.571, 0.472, 0.406, \
                               0.361, 0., 0.])

        x_ans_06t = np.array([0, 0.0645, 0.0645, 0.215, 0.215, 0.337, 0.339, \
                              0.341, 0.343, 0.346, 0.348, 0.350, 0.352, 0.354, \
                              0.357, 0.357, 0.461, 0.461, 1])

        Sg_ans_06t = np.array([1, 1, 0.775, 0.775, 0.571, 0.571, 0.560, 0.545, \
                               0.537, 0.517, 0.505, 0.488, 0.479, 0.459, 0.433, \
                               0.405, 0.369, 0., 0])

        x_ans_09t = np.array([0, 0.113086, 0.113333, 0.381975, 0.381728, 0.381235,
                          0.380988, 0.382222, 0.581975, 0.597778, 0.60321, 0.606667, \
                          0.609136, 0.613333, 0.615802, 0.62, 0.624198, 0.626914, \
                          0.628889, 0.631111, 0.632, 0.633333, 0.633333, 0.63358, \
                          0.644444, 0.817778, 0.817778, 1])

        Sg_ans_09t = np.array([1, 1, 0.775536, 0.775216,  0.736471, 0.62472, 0.587576,   \
                           0.570285, 0.570605, 0.570925, 0.552994, 0.542107, 0.533141, \
                           0.519052, 0.509766, 0.495677, 0.478066, 0.463977, 0.456612, \
                           0.443484, 0.435479, 0.425232, 0.415626, 0.403778, 0.400576, \
                           0.360231, 0., 0])

        x_ans_Sg = np.array([0, 0.219, 0.219, 0.736, 0.737, 1.153, 1.156, 1.160, \
                             1.163, 1.165, 1.169, 1.172, 1.175, 1.178, 1.181, \
                             1.184, 1.187, 1.191, 1.194, 1.197, 1.203, 1.206, \
                             1.209, 1.212, 1.215, 1.219, 1.222, 1.222, 1.576, \
                             1.577, 2])
        Sg_ans = np.array([1, 1, 0.775, 0.775, 0.57, 0.57, 0.565, 0.56, 0.555, \
                           0.55, 0.545, 0.54, 0.535, 0.53, 0.522, 0.517, 0.512, \
                           0.507, 0.507, 0.495, 0.48, 0.472, 0.465, 0.457, 0.447, \
                           0.432, 0.412, 0.400, 0.36, 0, 0])

        x_ans_zCO2 = np.array([0, 0.22, 0.22, 0.737, 0.737, 1.156, 1.156, 1.161, \
                               1.165, 1.169, 1.178, 1.183, 1.191, 1.199, 1.204, \
                               1.209, 1.214, 1.219, 1.223, 1.566, 1.575, 1.576, \
                               1.999])

        zCO2_ans = np.array([0.1, 0.1, 0.0916, 0.0916, 0.076, 0.082, 0.086, 0.101, \
                             0.119, 0.139, 0.166, 0.192, 0.211, 0.235, 0.250, \
                             0.259, 0.269, 0.277, 0.287, 0.284, 0.277, 0.2, 0.2])


        '------------------------------- FOU ----------------------------------'
        datas = np.load('flying/results_4k_Nekoui_case1_100x1x1_IMPEC_FOU_137.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_100_FOU = data[10]
            n=100
            x_100 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_100_FOU = data[7]
            t_100_FOU = data[2]

        datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_FOU_279.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_200_FOU = data[10]
            n=200
            x_200 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_200_FOU = data[7]
            t_200_FOU = data[2]

        datas = np.load('flying/results_4k_Nekoui_case1_350x1x1_IMPEC_FOU_494.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_350_FOU = data[10]
            n=350
            x_350 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_350_FOU = data[7]
            t_350_FOU = data[2]


        datas = np.load('flying/results_4k_Nekoui_case1_1000x1x1_IMPEC_FOU_1414.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_1000_FOU = data[10]
            n=1000
            x_1000 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_1000_FOU = data[7]

        datas = np.load('flying/results_4k_Nekoui_case1_2000x1x1_IMPEC_FOU_2838.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_2000_FOU = data[10]
            n=2000
            x_2000 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_2000_FOU = data[7]
            t_2000_FOU = data[2]

        datas = np.load('flying/results_4k_Nekoui_case1_5000x1x1_IMPEC_FOU_23167.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_5000_FOU = data[10]
            n=5000
            x_5000 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_5000_FOU = data[7]

        '---------------------------------- MUSCL -----------------------------'
        #datas = np.load('flying/results_4k_Nekoui_case1_100x1x1_IMPEC_MUSCL_LLF_minmod_473.npy', allow_pickle=True)
        #datas = np.load('flying/results_4k_Nekoui_case1_100x1x1_IMPEC_MUSCL_LLF_VA_461.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_100x1x1_IMPEC_MUSCL_LLF_minmod_RK3t_499.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            #CFL=0.7
            z_100_MUSCL = data[10]
            Sg_100_MUSCL = data[7]
            t_100_MUSCL = data[2]

        #datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_MUSCL_LLF_minmod_956.npy', allow_pickle=True)
        #datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_MUSCL_LLF_VA_934.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_MUSCL_LLF_minmod_RK3t_999.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_200_MUSCL = data[10]
            Sg_200_MUSCL = data[7]
            t_200_MUSCL = data[2]


        #datas = np.load('flying/results_4k_Nekoui_case1_300x1x1_IMPEC_MUSCL_LLF_minmod_1442.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_300x1x1_IMPEC_MUSCL_LLF_minmod_RK3t_1500.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_300_MUSCL = data[10]
            Sg_300_MUSCL = data[7]
            t_300_MUSCL = data[2]
            n=300
            x_300 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)

        #datas = np.load('flying/results_4k_Nekoui_case1_350x1x1_IMPEC_MUSCL_LLF_VA_1649.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_350x1x1_IMPEC_MUSCL_LLF_VA_RK3t_1750.npy', allow_pickle=True)
        #datas = np.load('flying/results_4k_Nekoui_case1_350x1x1_IMPEC_MUSCL_LLF_minmod_RK3t_1750.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_350_MUSCL = data[10]
            Sg_350_MUSCL = data[7]
            t_350_MUSCL = data[2]
            n=350
            x_350 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)

        #datas = np.load('flying/results_4k_Nekoui_case1_400x1x1_IMPEC_MUSCL_LLF_minmod_1926.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_400x1x1_IMPEC_MUSCL_LLF_VA_1887.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_400_MUSCL = data[10]
            Sg_400_MUSCL = data[7]
            t_400_MUSCL = data[2]
            n=400
            x_400 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)


        '---------------------------------- FR-P1 -----------------------------'
        #datas = np.load('flying/results_4k_Nekoui_case1_50x1x1_IMPEC_FR2_LLF_RK3_222.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_50x1x1_IMPEC_FR2_LLF_RK3t_251.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_50_FR2 = data[10]
            Sg_50_FR2 = data[7]
            n=50
            x_50 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)


        #datas = np.load('flying/results_4k_Nekoui_case1_100x1x1_IMPEC_FR2_LLF_RK3_455.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_100x1x1_IMPEC_FR2_LLF_RK3t_502.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-2:-1]:
            z_100_FR2 = data[10]
            Sg_100_FR2 = data[7]
            t_100_FR2 = data[2]


        #datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_FR2_LLF_RK3_924.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_FR2_LLF_RK3t_1002.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_200_FR2 = data[10]
            Sg_200_FR2 = data[7]
            t_200_FR2 = data[2]


        #datas = np.load('flying/results_4k_Nekoui_case1_400x1x1_IMPEC_FR2_LLF_RK3_1864.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_400x1x1_IMPEC_FR2_LLF_RK3t_2005.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_400_FR2 = data[10]
            Sg_400_FR2 = data[7]
            n=400
            x_400 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)



        '---------------------------------- FR-P2 -----------------------------'
        #datas = np.load('flying/results_4k_Nekoui_case1_50x1x1_IMPEC_FR3_LLF_RK3_322.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_50x1x1_IMPEC_FR3_LLF_RK3t_332.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-2:-1]:
            z_50_FR3 = data[10]
            Sg_50_FR3 = data[7]

        #datas = np.load('flying/results_4k_Nekoui_case1_100x1x1_IMPEC_FR3_LLF_RK3_679.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_100x1x1_IMPEC_FR3_LLF_RK3t_687.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-2:-1]:
            z_100_FR3 = data[10]
            Sg_100_FR3 = data[7]
            t_100_FR3 = data[2]

        #datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_FR3_LLF_RK3_1386.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_FR3_LLF_RK3t_1398.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-2:-1]:
            z_200_FR3 = data[10]
            Sg_200_FR3 = data[7]
            t_200_FR3 = data[2]

        #datas = np.load('flying/results_4k_Nekoui_case1_400x1x1_IMPEC_FR3_LLF_RK3_2810.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_400x1x1_IMPEC_FR3_LLF_RK3t_2819.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-2:-1]:
            z_400_FR3 = data[10]
            Sg_400_FR3 = data[7]


        '---------------------------------- FR-P3 -----------------------------'
        #datas = np.load('flying/results_4k_Nekoui_case1_50x1x1_IMPEC_FR4_LLF_RK3_430.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_50x1x1_IMPEC_FR4_LLF_RK3t_434.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_50_FR4 = data[10]
            Sg_50_FR4 = data[7]

        #datas = np.load('flying/results_4k_Nekoui_case1_100x1x1_IMPEC_FR4_LLF_RK3_887.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_100x1x1_IMPEC_FR4_LLF_RK3t_897.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            z_100_FR4 = data[10]
            n=100
            x_100 = np.linspace(0+L/(2*n),L*(1-1/(2*n)),n)
            Sg_100_FR4 = data[7]
            t_100_FR4 = data[2]


        #datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_FR4_LLF_RK3_1827.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_FR4_LLF_RK3t_1834.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-2:-1]:
            z_200_FR4 = data[10]
            Sg_200_FR4 = data[7]
            t_200_FR4 = data[2]


        #datas = np.load('flying/results_4k_Nekoui_case1_400x1x1_IMPEC_FR4_LLF_RK3_3717.npy', allow_pickle=True)
        datas = np.load('flying/results_4k_Nekoui_case1_400x1x1_IMPEC_FR4_LLF_RK3t_3725.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-2:-1]:
            z_400_FR4 = data[10]
            Sg_400_FR4 = data[7]

        """ tests for 3t-steps"""

        datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_3ts_FR2_LLF_RK3t_1287.npy', allow_pickle=True)
        i=0
        z_3t_200_FR2 = np.empty((4,200,3))
        Sg_3t_200_FR2 = np.empty((200,3))
        P_3t_200_FR2 = np.empty((200,3))
        for data in datas[1:]:
            P_3t_200_FR2[...,i] = data[4]
            z_3t_200_FR2[...,i] = data[10]
            Sg_3t_200_FR2[...,i] = data[7]
            i+=1

        datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_3ts_FR3_LLF_RK3t_1717.npy', allow_pickle=True)
        i=0
        z_3t_200_FR3 = np.empty((4,200,3))
        Sg_3t_200_FR3 = np.empty((200,3))
        P_3t_200_FR3 = np.empty((200,3))
        for data in datas[1:]:
            P_3t_200_FR3[...,i] = data[4]
            z_3t_200_FR3[...,i] = data[10]
            Sg_3t_200_FR3[...,i] = data[7]
            i+=1

        datas = np.load('flying/results_4k_Nekoui_case1_200x1x1_IMPEC_3ts_FR4_LLF_RK3t_1938.npy', allow_pickle=True)
        i=0
        z_3t_200_FR4 = np.empty((4,200,3))
        Sg_3t_200_FR4 = np.empty((200,3))
        P_3t_200_FR4 = np.empty((200,3))
        for data in datas[1:]:
            P_3t_200_FR4[...,i] = data[4]
            z_3t_200_FR4[...,i] = data[10]
            Sg_3t_200_FR4[...,i] = data[7]
            i+=1


        plt.figure(1)
        s=3
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s)
        plt.plot(x_50/0.5, z_50_FR2[1,:], '-bo', mfc='none', markersize=s)
        plt.plot(x_100/0.5, z_100_FR2[1,:], '-g*', mfc='none', markersize=s)
        #plt.plot(x_150/0.5, z_150_FR2[1,:], '-.c', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR2[1,:], '-..r', mfc='none', markersize=s)
        plt.plot(x_400/0.5, z_400_FR2[1,:], '--k', mfc='none', markersize=s)
        #plt.plot(x_500/0.5, z_500_FR2[1,:], '--k', mfc='none', markersize=s)
        #plt.xlim(0.6,1.75)
        #plt.ylim(0., 0.4)
        plt.legend(('Reference', '50 CV', '100 CV', '200 CV', '400 CV'))
        plt.title('Results for FR-P1')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_FR2.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s)
        plt.plot(x_50/0.5, z_50_FR2[1,:], '-bo', mfc='none', markersize=s)
        plt.plot(x_100/0.5, z_100_FR2[1,:], '-g*', mfc='none', markersize=s)
        #plt.plot(x_150/0.5, z_150_FR2[1,:], '-.c', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR2[1,:], '-..r', mfc='none', markersize=s)
        plt.plot(x_400/0.5, z_400_FR2[1,:], '--k', mfc='none', markersize=s)
        #plt.plot(x_500/0.5, z_500_FR2[1,:], '--k', mfc='none', markersize=s)
        plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('Reference', '50 CV', '100 CV', '200 CV', '400 CV'))
        plt.title('Results for FR-P1')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_FR2_zoom.png')
        plt.close()

        plt.figure(2)
        s=3
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s+1)
        plt.plot(x_50/0.5, z_50_FR3[1,:], '-bo', mfc='none', markersize=s)
        plt.plot(x_100/0.5, z_100_FR3[1,:], '-g*', mfc='none', markersize=s)
        #plt.plot(x_150/0.5, z_150_FR3[1,:], '-.c', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR3[1,:], '-..r', mfc='none', markersize=s)
        plt.plot(x_400/0.5, z_400_FR3[1,:], '--k', mfc='none', markersize=s)
        #plt.plot(x_500/0.5, z_500_FR2[1,:], '--k', mfc='none', markersize=s)
        #plt.xlim(0.6,1.75)
        #plt.ylim(0., 0.4)
        plt.legend(('Reference', '50 CV', '100 CV', '200 CV', '400 CV'))
        plt.title('Results for FR-P2')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_FR3.png')
        plt.close()

        plt.figure(2)
        s=3
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s+1)
        plt.plot(x_50/0.5, z_50_FR3[1,:], '-bo', mfc='none', markersize=s)
        plt.plot(x_100/0.5, z_100_FR3[1,:], '-g*', mfc='none', markersize=s)
        #plt.plot(x_150/0.5, z_150_FR3[1,:], '-.c', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR3[1,:], '-..r', mfc='none', markersize=s)
        plt.plot(x_400/0.5, z_400_FR3[1,:], '--k', mfc='none', markersize=s)
        #plt.plot(x_500/0.5, z_500_FR2[1,:], '--k', mfc='none', markersize=s)
        plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('Reference', '50 CV', '100 CV', '200 CV', '400 CV'))
        plt.title('Results for FR-P2')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_FR3_zoom.png')
        plt.close()

        plt.figure(3)
        s=3
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s+1)
        plt.plot(x_50/0.5, z_50_FR4[1,:], '-bo', mfc='none', markersize=s)
        plt.plot(x_100/0.5, z_100_FR4[1,:], '-g*', mfc='none', markersize=s)
        #plt.plot(x_150/0.5, z_150_FR4[1,:], '-.c', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR4[1,:], '-..r', mfc='none', markersize=s)
        plt.plot(x_400/0.5, z_400_FR4[1,:], '--k', mfc='none', markersize=s)
        #plt.plot(x_300/0.5, z_300_FR4[1,:], '--k', mfc='none', markersize=s)
        #plt.plot(x_400/0.5, z_400_FR4[1,:], '--k', mfc='none', markersize=s)
        #plt.xlim(0.6,1.75)
        #plt.ylim(0., 0.4)
        plt.legend(('Reference', '50 CV', '100 CV', '200 CV', '400 CV'))
        plt.title('Results for FR-P3')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_FR4.png')
        plt.close()

        plt.figure(3)
        s=3
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s+1)
        plt.plot(x_50/0.5, z_50_FR4[1,:], '-bo', mfc='none', markersize=s)
        plt.plot(x_100/0.5, z_100_FR4[1,:], '-g*', mfc='none', markersize=s)
        #plt.plot(x_150/0.5, z_150_FR4[1,:], '-.c', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR4[1,:], '-..r', mfc='none', markersize=s)
        plt.plot(x_400/0.5, z_400_FR4[1,:], '--k', mfc='none', markersize=s)
        #plt.plot(x_500/0.5, z_500_FR4[1,:], '--k', mfc='none', markersize=s)
        plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('Reference', '50 CV', '100 CV', '200 CV', '400 CV'))
        plt.title('Results for FR-P3')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_FR4_zoom.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_ans_Sg, Sg_ans, '-k', mfc='none', markersize=s)
        plt.plot(x_200/0.5, Sg_200_FR2, '-gv', mfc='none', markersize=s)
        plt.plot(x_200/0.5, Sg_200_FR3, '-yo', mfc='none', markersize=s)
        plt.plot(x_200/0.5, Sg_200_FR4, '-.r', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        #plt.xlim(1.1,1.6)
        #plt.ylim(0.18, 0.3)
        plt.legend(('Reference', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.title('Results for 200x1x1 mesh')
        plt.grid()
        plt.ylabel('$Sg$')
        #plt.xlim(0.4, 0.9)
        plt.xlabel('Distance')
        plt.savefig('results/compositional/SPE_paper/4k_Sg_Nekoui_case1_NVCM_200.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_ans_Sg, Sg_ans, '-k', mfc='none', markersize=s)
        plt.plot(x_100/0.5, Sg_100_FR2, '-gv', mfc='none', markersize=s)
        plt.plot(x_100/0.5, Sg_100_FR3, '-yo', mfc='none', markersize=s)
        plt.plot(x_100/0.5, Sg_100_FR4, '-.r', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        #plt.xlim(1.1,1.6)
        #plt.ylim(0.18, 0.3)
        plt.legend(('Reference', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.title('Results for 100x1x1 mesh')
        plt.grid()
        plt.ylabel('$Sg$')
        #plt.xlim(0.4, 0.9)
        plt.xlabel('Distance')
        plt.savefig('results/compositional/SPE_paper/4k_Sg_Nekoui_case1_NVCM_100.png')
        plt.close()


        plt.figure(1)
        s=3
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s+1)
        plt.plot(x_100/0.5, z_100_FOU[1,:], '-ro', mfc='none', markersize=s)
        plt.plot(x_100/0.5, z_100_MUSCL[1,:], '-bs', mfc='none', markersize=s)
        #plt.plot(x_100/0.5, z_100_FR2[1,:], '-gv', mfc='none', markersize=s)
        #plt.plot(x_100/0.5, z_100_FR3[1,:], '-.g', mfc='none', markersize=s)
        plt.plot(x_100/0.5, z_100_FR4[1,:], '-r*', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('Reference', 'FOU (100 CV)', 'MUSCL (100 CV)', 'FR-P3 (100 CV)'))
        #plt.title('Results for 100x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_100_1.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_100/0.5, z_100_FOU[1,:], '-ro', mfc='none', markersize=s)
        plt.plot(x_100/0.5, z_100_MUSCL[1,:], '-bs', mfc='none', markersize=s)
        #plt.plot(x_100/0.5, z_100_FR2[1,:], '-gv', mfc='none', markersize=s)
        plt.plot(x_100/0.5, z_100_FR3[1,:], '-.g', mfc='none', markersize=s)
        #plt.plot(x_100/0.5, z_100_FR4[1,:], '-.c', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s+1)
        plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('FOU (100 CV)', 'MUSCL (100 CV)', 'FR-P2 (100 CV)', 'Reference'))
        #plt.title('Results for 100x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_100_2.png')
        plt.close()


        plt.figure(1)
        s=3
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FOU[1,:], '-ro', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_MUSCL[1,:], '-bs', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_FR2[1,:], '-gv', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_FR3[1,:], '-.g', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR4[1,:], '-r*', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('Reference', 'FOU (200 CV)', 'MUSCL (200 CV)', 'FR-P3 (200 CV)'))
        #plt.title('Results for 200x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        #plt.xlim(0., 2)
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_200_1.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_200/0.5, z_200_FOU[1,:], '-ro', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_MUSCL[1,:], '-bs', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_FR2[1,:], '-gv', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR3[1,:], '-.g', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_FR4[1,:], '-.c', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s)
        plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('FOU (200 CV)', 'MUSCL (200 CV)', 'FR-P2 (200 CV)', 'Reference'))
        #plt.title('Results for 200x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        #plt.xlim(0., 2)
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_200_2.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s)

        plt.plot(x_2000/0.5, z_2000_FOU[1,:], '--b', mfc='none', markersize=s)
        plt.plot(x_350/0.5, z_350_MUSCL[1,:], '-.m', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_FR2[1,:], '-gv', mfc='none', markersize=s)
        #plt.plot(x_250/0.5, z_250_FR3[1,:], '-mo', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR4[1,:], '-r*', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('Reference', 'FOU (2000 CV)', 'MUSCL (350 CV)', 'FR-P3 (200 CV)'))
        #plt.title('Results for 250x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        #plt.xlim(0., 2)
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_BOA1.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s)

        plt.plot(x_2000/0.5, z_2000_FOU[1,:], '--b', mfc='none', markersize=s)
        plt.plot(x_300/0.5, z_300_MUSCL[1,:], '-.m', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_FR2[1,:], '-gv', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR3[1,:], '-y<', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_FR4[1,:], '-.c', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('Reference', 'FOU (2000 CV)', 'MUSCL (300 CV)', 'FR-P2   (200 CV)'))
        #plt.title('Results for 250x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{CO_2}$')
        #plt.xlim(0., 2)
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_CO2_Nekoui_case1_NVCM_BOA2.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_ans_03t, Sg_ans_03t, '-k', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR2[...,0], '-gv', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR3[...,0], '-yo', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR4[...,0], '-.r', mfc='none', markersize=s)
        plt.legend(('Reference', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.title('Results for 200x1x1 mesh (t=636s)')
        plt.grid()
        plt.ylabel('$Sg$')
        plt.xlim(0., 0.3)
        plt.xlabel('Distance')
        plt.savefig('results/compositional/SPE_paper/4k_Sg_Nekoui_case1_NVCM_FRs_200_03t.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_ans_06t, Sg_ans_06t, '-k', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR2[...,1], '-gv', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR3[...,1], '-yo', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR4[...,1], '-.r', mfc='none', markersize=s)
        plt.legend(('Reference', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.title('Results for 200x1x1 mesh (t=1272s)')
        plt.grid()
        plt.ylabel('$Sg$')
        plt.xlim(0., 0.5)
        plt.xlabel('Distance')
        plt.savefig('results/compositional/SPE_paper/4k_Sg_Nekoui_case1_NVCM_FRs_200_06t.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_ans_09t, Sg_ans_09t, '-k', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR2[...,2], '-gv', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR3[...,2], '-yo', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR4[...,2], '-.r', mfc='none', markersize=s)
        plt.legend(('Reference', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.title('Results for 200x1x1 mesh (t=1908s)')
        plt.grid()
        plt.ylabel('$Sg$')
        plt.xlim(0., 0.85)
        plt.xlabel('Distance')
        plt.savefig('results/compositional/SPE_paper/4k_Sg_Nekoui_case1_NVCM_FRs_200_09t.png')
        plt.close()


        plt.figure(1)
        s=3
        plt.plot(x_ans_03t, Sg_ans_03t, '-k', mfc='none', markersize=s)
        plt.plot(x_ans_06t, Sg_ans_06t, '--k', mfc='none', markersize=s)
        plt.plot(x_ans_09t, Sg_ans_09t, '-..k', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR2[...,0], '-gv', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR2[...,1], '-yo', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR2[...,2], '-rs', mfc='none', markersize=s)
        plt.legend(('Reference (t=636s)', 'Reference (t=1272s)', 'Reference (t=1908s)', \
                'FR-P1 (t=636s)', 'FR-P1 (t=1272s)', 'FR-P1 (t=1908s)'))
        plt.title('Results for 200x1x1 mesh FR-P1')
        plt.grid()
        plt.ylabel('$Sg$')
        plt.xlim(0., 0.85)
        plt.xlabel('Distance')
        plt.savefig('results/compositional/SPE_paper/4k_Sg_Nekoui_case1_NVCM_FR2_200_ts.png')
        plt.close()


        plt.figure(1)
        s=3
        plt.plot(x_ans_03t, Sg_ans_03t, '-k', mfc='none', markersize=s)
        plt.plot(x_ans_06t, Sg_ans_06t, '--k', mfc='none', markersize=s)
        plt.plot(x_ans_09t, Sg_ans_09t, '-..k', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR3[...,0], '-gv', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR3[...,1], '-yo', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR3[...,2], '-rs', mfc='none', markersize=s)
        plt.legend(('Reference (t=636s)', 'Reference (t=1272s)', 'Reference (t=1908s)', \
                'FR-P2 (t=636s)', 'FR-P2 (t=1272s)', 'FR-P2 (t=1908s)'))
        plt.title('Results for 200x1x1 mesh FR-P2')
        plt.grid()
        plt.ylabel('$Sg$')
        plt.xlim(0., 0.85)
        plt.xlabel('Distance')
        plt.savefig('results/compositional/SPE_paper/4k_Sg_Nekoui_case1_NVCM_FR3_200_ts.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_ans_03t, Sg_ans_03t, '-k', mfc='none', markersize=s)
        plt.plot(x_ans_06t, Sg_ans_06t, '--k', mfc='none', markersize=s)
        plt.plot(x_ans_09t, Sg_ans_09t, '-..k', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR4[...,0], '-gv', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR4[...,1], '-yo', mfc='none', markersize=s)
        plt.plot(x_200, Sg_3t_200_FR4[...,2], '-rs', mfc='none', markersize=s)
        plt.legend(('Reference (t=636s)', 'Reference (t=1272s)', 'Reference (t=1908s)', \
                'FR-P3 (t=636s)', 'FR-P3 (t=1272s)', 'FR-P3 (t=1908s)'))
        plt.title('Results for 200x1x1 mesh FR-P3')
        plt.grid()
        plt.ylabel('$Sg$')
        plt.xlim(0., 0.85)
        plt.xlabel('Distance')
        plt.savefig('results/compositional/SPE_paper/4k_Sg_Nekoui_case1_NVCM_FR4_200_ts.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_200, P_3t_200_FR4[...,0], '-gv', mfc='none', markersize=s)
        plt.plot(x_200, P_3t_200_FR4[...,1], '-yo', mfc='none', markersize=s)
        plt.plot(x_200, P_3t_200_FR4[...,2], '-.r', mfc='none', markersize=s)
        plt.legend(('Reference (t=636s)', 'Reference (t=1272s)', 'Reference (t=1908s)', \
                'FR-P3 (t=636s)', 'FR-P3 (t=1272s)', 'FR-P3 (t=1908s)'))
        plt.title('Results for 200x1x1 mesh FR-P3')
        plt.grid()
        plt.ylabel('$P$')
        plt.xlim(0., 0.85)
        plt.xlabel('Distance')
        plt.savefig('results/compositional/SPE_paper/4k_P_Nekoui_case1_NVCM_FR4_200_ts.png')
        plt.close()

        import pdb; pdb.set_trace()

        plt.figure(1)
        s=3
        plt.plot(x_200/0.5, z_200_FOU[0,:], '-ro', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_MUSCL[0,:], '-bs', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_FR2[0,:], '-gv', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR3[0,:], '-.c', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR4[0,:], '--m', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        #plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s)
        #plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('FOU (200 CV)', 'MUSCL (200 CV)', 'FR-P2 (200 CV)', 'Reference'))
        #plt.title('Results for 200x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C1}$')
        #plt.xlim(0., 2)
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_C1_Nekoui_case1_NVCM_200.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_1000/0.5, z_1000_FOU[2,:], '-ro', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_MUSCL[2,:], '-bs', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_FR2[2,:], '-gv', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR3[2,:], '-.c', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR4[2,:], '--m', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        #plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s)
        #plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('FOU (200 CV)', 'MUSCL (200 CV)', 'FR-P2 (200 CV)', 'Reference'))
        #plt.title('Results for 200x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C4}$')
        #plt.xlim(0., 2)
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_C4_Nekoui_case1_NVCM_200.png')
        plt.close()

        plt.figure(1)
        s=3
        plt.plot(x_1000/0.5, z_1000_FOU[3,:], '-ro', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_MUSCL[3,:], '-bs', mfc='none', markersize=s)
        #plt.plot(x_200/0.5, z_200_FR2[3,:], '-gv', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR3[3,:], '-.c', mfc='none', markersize=s)
        plt.plot(x_200/0.5, z_200_FR4[3,:], '--m', mfc='none', markersize=s)
        #plt.plot(x_8000, z_8000_FOU[1,:], '-k', mfc='none', markersize=s)
        #plt.plot(x_ans_zCO2, zCO2_ans, '-k', mfc='none', markersize=s)
        #plt.xlim(1,1.7)
        #plt.ylim(0., 0.4)
        plt.legend(('FOU (200 CV)', 'MUSCL (200 CV)', 'FR-P2 (200 CV)', 'Reference'))
        #plt.title('Results for 200x1x1 mesh')
        plt.grid()
        plt.ylabel('$z_{C10}$')
        #plt.xlim(0., 2)
        plt.xlabel(r'$x/\tau$')
        plt.savefig('results/compositional/SPE_paper/4k_C10_Nekoui_case1_NVCM_200.png')
        plt.close()

        import pdb; pdb.set_trace()
