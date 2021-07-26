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
        x_axis_xC30 = np.linspace(0, 21.4, 100)
        xC30 = np.zeros(len(x_axis_xC30))
        x_axis_xC3 = np.array([21.5484, 21.7646, 21.9171, 21.968, 22.1459, 22.5655, 24.8665, 26.2904,
            26.5446])
        xC3 = np.array([2.867384e-4, 0.00258065, 0.00659498, 0.0219355, 0.0840143, 0.317563, 0.308817, 0.303656,
            0.])

        '''x_axis_xCH4 = np.array([23.1051, 23.1051, 23.1476, 23.9119, 24.0817, 24.4958, 24.9522,  25.7909,
                        26.8737, 27.7495, 27.7866, 27.8822, 27.9989, 28.0945, 28.259,  28.4289, 28.5881,
                        28.7686, 28.9756, 29.241, 29.6125, 29.8779, 30.483])
        xCH4 = np.array([ 1.743896e-4, 0.196886, 0.332735, 0.332735, 0.331166, 0.329248, 0.325237, 0.321749, 0.318784, 0.318261,
                        0.307623, 0.239437, 0.179098, 0.139686, 0.0953911, 0.0636522, 0.0430742, 0.0273792,
                        0.0162182, 0.00767314, 0.00348779, 0.00156951, 0.0])'''

        x_axis_xC31 = np.linspace(x_axis_xC3[-1], 50, 100)
        xC31 = np.zeros(len(x_axis_xC31))

        x_axis_xC3 = np.concatenate((x_axis_xC30,x_axis_xC3))
        x_axis_xC3 = np.concatenate((x_axis_xC3,x_axis_xC31))
        xC3 = np.concatenate((xC30,xC3))
        xC3 = np.concatenate((xC3,xC31))


        datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw_852.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_500 = data[6]
            Sg_500 = data[7]
            z_500 = data[10]
            xkj_500 = data[13]
            P_500 = data[4]

            xkj_500[:,1,Sg_500==0] = 0
            xkj_500[:,0,So_500==0] = 0
            x_500 = np.linspace(0,50,len(So_500))

        datas = np.load('flying/results_Hoteit_Firoo_3k_5000_upw_7382.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_5000 = data[6]
            Sg_5000 = data[7]
            z_5000 = data[10]
            xkj_5000 = data[13]
            P_5000 = data[4]

            xkj_5000[:,1,Sg_5000==0] = 0
            xkj_5000[:,0,So_5000==0] = 0
            x_5000 = np.linspace(0,50,len(So_5000))



        plt.figure(1)
        plt.plot(x_500, So_500, 'k', x_5000, So_5000, 'b')
        plt.grid()
        plt.ylabel('Oil Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_oil_saturation.png')

        plt.figure(7)
        plt.plot(x_500, Sg_500)
        plt.grid()
        plt.ylabel('Gas Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_gas_saturation.png')

        plt.figure(2)
        plt.plot(x_500, xkj_500[0,1,:], 'b')
        plt.grid()
        plt.ylabel('Methane molar fraction in vapor phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x.png')


        plt.figure(3)
        plt.plot(x_500, xkj_500[0,0,:], 'b')
        plt.grid()
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_y.png')

        plt.figure(4)
        plt.plot(x_500, xkj_500[2,1,:], 'b')
        plt.plot(x_axis_xC3, xC3, 'k')
        plt.plot(x_CMG, yC3H8_CMG, 'y')
        plt.plot(x_5000, xkj_5000[2,1,:], 'r')
        plt.legend(('FOU 500', 'Reference', 'CMG-5000', 'FOU-5000'))
        plt.grid()
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_propane_y.png')

        plt.figure(5)
        plt.plot(x_500, xkj_500[2,0,:])
        plt.grid()
        plt.ylabel('Propane molar fraction in liquid phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_propane_x.png')


        plt.figure(6)
        plt.plot(x_500, P_500)
        plt.grid()
        plt.ylabel('Pressure')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_pressure.png')
        import pdb; pdb.set_trace()
