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
        x_axis_xC30 = np.linspace(0, 21.4, 100)
        xC30 = np.zeros(len(x_axis_xC30))
        x_axis_xC3 = np.array([21.5484, 21.7646, 21.9171, 21.968, 22.1459, 22.5655, 24.8665, 26.2904,
            26.5446])
        xC3 = np.array([2.867384e-4, 0.00258065, 0.00659498, 0.0219355, 0.0840143, 0.317563, 0.308817, 0.303656,
            0.])

        x_axis_xC31 = np.linspace(x_axis_xC3[-1], 50, 100)
        xC31 = np.zeros(len(x_axis_xC31))

        x_axis_xC3 = np.concatenate((x_axis_xC30,x_axis_xC3))
        x_axis_xC3 = np.concatenate((x_axis_xC3,x_axis_xC31))
        xC3 = np.concatenate((xC30,xC3))
        xC3 = np.concatenate((xC3,xC31))

        datas = np.load('flying/results_Hoteit_Firoo_3k_200_IMPEC_FOU_847.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[1:2]:
            So_FOU_200 = data[6]
            Sg_FOU_200 = data[7]
            xkj_FOU_200 = data[13]
            xkj_FOU_200[:,1,Sg_FOU_200==0] = 0
            xkj_FOU_200[:,0,So_FOU_200==0] = 0
            n=200
            x_200 = np.linspace(0+1/(2*n),50-1/(2*n),n)


        datas = np.load('flying/results_Hoteit_Firoo_3k_200_IMPEC_MUSCL_1408.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[1:2]:
            So_MUSCL_200 = data[6]
            Sg_MUSCL_200 = data[7]
            xkj_MUSCL_200 = data[13]
            xkj_MUSCL_200[:,1,Sg_MUSCL_200==0] = 0
            xkj_MUSCL_200[:,0,So_MUSCL_200==0] = 0
            n=200
            x_200 = np.linspace(0+1/(2*n),50-1/(2*n),n)

        
        plt.figure(1)
        #plt.plot(x_500, xkj_500[2,1,:], 'b')
        plt.plot(x_200, xkj_FOU_200[2,1,:], 'y')
        plt.plot(x_200, xkj_MUSCL_200[2,1,:], 'b')
        plt.plot(x_axis_xC3, xC3, 'k')
        plt.plot(x_CMG, yC3H8_CMG, '-r')
        plt.legend(('FOU 200', 'MUSCL 200', 'Reference'))
        plt.grid()
        plt.xlim()
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/TCC2/3k_propane_y_HF.png')

        import pdb; pdb.set_trace()
