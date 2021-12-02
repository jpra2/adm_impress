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
        '''x_axis_xC30 = np.linspace(0, 21.4, 100)
        xC30 = np.zeros(len(x_axis_xC30))
        x_axis_xC3 = np.array([21.5484, 21.7646, 21.9171, 21.968, 22.1459, 22.5655, 24.8665, 26.2904,
            26.5446])
        xC3 = np.array([2.867384e-4, 0.00258065, 0.00659498, 0.0219355, 0.0840143, 0.317563, 0.308817, 0.303656,
            0.])'''

        x_axis_xCH4 = np.array([0, 23.1051, 23.1051, 23.1476, 23.9119, 24.0817, 24.4958, 24.9522,  25.7909,
                        26.8737, 27.7495, 27.7866, 27.8822, 27.9989, 28.0945, 28.259,  28.4289, 28.5881,
                        28.7686, 28.9756, 29.241, 29.6125, 29.8779, 30.483,50])
        xCH4 = np.array([0, 1.743896e-4, 0.196886, 0.332735, 0.332735, 0.331166, 0.329248, 0.325237, 0.321749, 0.318784, 0.318261,
                        0.307623, 0.239437, 0.179098, 0.139686, 0.0953911, 0.0636522, 0.0430742, 0.0273792,
                        0.0162182, 0.00767314, 0.00348779, 0.00156951, 0.0,0])

        #x_axis_xC31 = np.linspace(x_axis_xC3[-1], 50, 100)
        #xC31 = np.zeros(len(x_axis_xC31))

        '''x_axis_xC3 = np.concatenate((x_axis_xC30,x_axis_xC3))
        x_axis_xC3 = np.concatenate((x_axis_xC3,x_axis_xC31))
        xC3 = np.concatenate((xC30,xC3))
        xC3 = np.concatenate((xC3,xC31))'''


        f = interp1d(x_axis_xCH4,xCH4)

        """---------------------- Convergence Study -------------------------"""

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_IMPEC_FOU_NEW_230.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_200 = data[6]
            Sg_FOU_200 = data[7]
            xkj_FOU_200 = data[13]

            xkj_FOU_200[:,1,Sg_FOU_200==0] = 0
            xkj_FOU_200[:,0,So_FOU_200==0] = 0
            n = 200
            x_200 = np.linspace(0+50/(2*n),50-50/(2*n),n)
            t_FOU_200 = data[2]

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_500_IMPEC_FOU_1516.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_500 = data[6]
            Sg_FOU_500 = data[7]
            xkj_FOU_500 = data[13]

            xkj_FOU_500[:,1,Sg_FOU_500==0] = 0
            xkj_FOU_500[:,0,So_FOU_500==0] = 0
            n = 500
            x_500 = np.linspace(0+50/(2*n),50-50/(2*n),n)
            t_FOU_500 = data[2]

            #e250_L1_FOU = (sum(abs(f(x_250)-xkj_IMPEC_250[0,0,:]))*(1/n))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_IMPEC_MUSCL_LLF_NEW_1165.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_200 = data[6]
            Sg_MUSCL_200 = data[7]
            z_MUSCL_200 = data[10]
            xkj_MUSCL_200 = data[13]

            xkj_MUSCL_200[:,1,Sg_MUSCL_200==0] = 0
            xkj_MUSCL_200[:,0,So_MUSCL_200==0] = 0
            n = 200
            x_200 = np.linspace(0+50/(2*n),50-50/(2*n),n)
            t_MUSCL_200 = data[2]
            #e500_L1_FOU = (sum(abs(f(x_500)-xkj_IMPEC_500))*(1/n))
        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_250_IMPEC_MUSCL_2394.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[1:2]:
            So_MUSCL_250 = data[6]
            Sg_MUSCL_250 = data[7]
            xkj_MUSCL_250 = data[13]
            xkj_MUSCL_250[:,1,Sg_MUSCL_250==0] = 0
            xkj_MUSCL_250[:,0,So_MUSCL_250==0] = 0
            n=250
            x_250 = np.linspace(0+1/(2*n),50-1/(2*n),n)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_500_IMPEC_MUSCL_3278.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[1:2]:
            So_MUSCL_500 = data[6]
            Sg_MUSCL_500 = data[7]
            xkj_MUSCL_500 = data[13]
            xkj_MUSCL_500[:,1,Sg_MUSCL_500==0] = 0
            xkj_MUSCL_500[:,0,So_MUSCL_500==0] = 0
            n=500
            x_500 = np.linspace(0+1/(2*n),50-1/(2*n),n)


        plt.figure(1)
        plt.plot(x_200, xkj_FOU_200[0,0], 'r')
        plt.plot(x_200, xkj_MUSCL_200[0,0], 'b')
        #plt.plot(x_250, xkj_MUSCL_250[0,0], 'g')
        plt.plot(x_500, xkj_FOU_500[0,0], 'g')
        plt.plot(x_CMG, xCH4_CMG, 'k')
        plt.plot(x_axis_xCH4, xCH4, '-y')
        plt.xlim(20,30)
        plt.ylim(-0.05, 0.35)
        plt.grid()
        plt.legend(( 'FOU-200', 'MUSCL-200', 'FOU-500', 'Referência CMG', 'Referência MM'), loc=10)
        plt.ylabel('Fração molar do metano na fase líquida')
        plt.xlabel('Distância em x')
        plt.savefig('results/compositional/TCC2/3k_methane_x_MM_IMPEC_200.png')


        import pdb; pdb.set_trace()
