import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
x_CMG = np.loadtxt('x_3k_Firoo_CMG.txt')
xCH4_CMG = np.loadtxt('xCH4_3k_Firoo.txt')
yC3H8_CMG = np.loadtxt('yC3H8_3k_Firoo.txt')
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


        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_50_FOU_1684.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw_945.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_50_FOU = data[6]
            Sg_50_FOU = data[7]
            z_50_FOU = data[10]
            xkj_50_FOU = data[13]
            P_50_FOU = data[4]

            xkj_50_FOU[:,1,Sg_50_FOU==0] = 0
            xkj_50_FOU[:,0,So_50_FOU==0] = 0
            n=50
            x_50 = np.linspace(0+50/(2*n),50-50/(2*n),n)


        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_50_LLF_3453.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_50_LLF = data[6]
            Sg_50_LLF = data[7]
            z_50_LLF = data[10]
            xkj_50_LLF = data[13]
            P_50_LLF = data[4]
            xkj_50_LLF[:,1,Sg_50_LLF==0] = 0
            xkj_50_LLF[:,0,So_50_LLF==0] = 0
            n=50
            x_50 = np.linspace(0+50/(2*n),50-50/(2*n),n)

        plt.figure(1)
        plt.plot(x_50, xkj_50_LLF[0,0,:], 'b')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.xlim(20,35)
        plt.grid()
        plt.legend(('LLF-50', 'Reference'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_LLF.png')
        import pdb; pdb.set_trace()

        '''plt.figure(4)
        #plt.plot(x_500, xkj_500[2,1,:], 'b')
        plt.plot(x_500, xkj_500_MUSCL[2,1,:], 'm')
        plt.plot(x_axis_xC3, xC3, 'k')
        plt.plot(x_CMG, yC3H8_CMG, 'y')
        #plt.plot(x_500,xkj_500_FR2_euler[2,1,:],'b')
        plt.plot(x_500,xkj_500_FR3_euler[2,1,:],'g')

        plt.plot(x_5000, xkj_5000[2,1,:], 'r')
        plt.legend(('MUSCL 500', 'Reference', 'FR2-500', 'FR3-500', 'FOU-5000'))
        plt.grid()
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_propane_y.png')'''
        import pdb; pdb.set_trace()
