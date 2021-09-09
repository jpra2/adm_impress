import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
x_CMG = np.loadtxt('x_3k_MM_CMG.txt')
xCH4_CMG = np.loadtxt('xCH4_3k_MM_CMG.txt')
for  arq in arquivos:
    if  arq.startswith(name):
        '''-------------------------MUSCL LLF RESULTS------------------------'''
        x_axis_xCH40 = np.linspace(0, 23.1051, 100)
        xCH40 = np.zeros(len(x_axis_xCH40))
        '''x_axis_xCH4 = np.array([22.3749, 22.4711, 22.6019, 22.767, 22.8144, 23.5318, 25.4384, 26.7638,
                        27.1177, 27.1588, 27.2937, 27.4269])
        xCH4 = np.array([0., 0.0582055, 0.216829, 0.34198, 0.346629, 0.346815, 0.343282, 0.336216,
                        0.122362, 0.0643422, 0.0126453, 0.])'''

        x_axis_xCH4 = np.array([23.1051, 23.1051, 23.1476, 23.9119, 24.0817, 24.4958, 24.9522,  25.7909,
                        26.8737, 27.7495, 27.7866, 27.8822, 27.9989, 28.0945, 28.259,  28.4289, 28.5881,
                        28.7686, 28.9756, 29.241, 29.6125, 29.8779, 30.483])
        xCH4 = np.array([ 1.743896e-4, 0.196886, 0.332735, 0.332735, 0.331166, 0.329248, 0.325237, 0.321749, 0.318784, 0.318261,
                        0.307623, 0.239437, 0.179098, 0.139686, 0.0953911, 0.0636522, 0.0430742, 0.0273792,
                        0.0162182, 0.00767314, 0.00348779, 0.00156951, 0.0])

        x_axis_xCH41 = np.linspace(x_axis_xCH4[-1], 50, 100)
        xCH41 = np.zeros(len(x_axis_xCH41))

        x_axis_xCH4 = np.concatenate((x_axis_xCH40,x_axis_xCH4))
        x_axis_xCH4 = np.concatenate((x_axis_xCH4,x_axis_xCH41))
        xCH4 = np.concatenate((xCH40,xCH4))
        xCH4 = np.concatenate((xCH4,xCH41))

        datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw_82170.npy', allow_pickle=True)
        datas = np.load('flying/results_Hoteit_Firoo_3k_200_FR_278977.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw_16187.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_5000_upw_321408.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So = data[6]
            Sg = data[7]
            z = data[10]
            xkj = data[13]
            P = data[4]

            xkj[:,1,Sg==0] = 0
            xkj[:,0,So==0] = 0
            x = np.linspace(0,50,len(So))

        datas = np.load('flying/results_Hoteit_Firoo_3k_5000_upw_321408.npy', allow_pickle=True)
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
        plt.plot(x, So)
        plt.grid()
        plt.ylabel('Oil Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_oil_saturation.png')

        plt.figure(7)
        plt.plot(x, Sg)
        plt.grid()
        plt.ylabel('Gas Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_gas_saturation.png')

        plt.figure(2)
        plt.plot(x, xkj[0,0,:], 'b', x_5000, xkj_5000[0,0,:], 'r')
        plt.plot(x_CMG, xCH4_CMG, 'y')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.legend(('FOU 500', 'FOU 5000', 'CMG-5000', 'Reference'))
        plt.grid()
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x.png')

        plt.figure(8)
        #plt.plot(x, xkj[0,0,:], 'b', x_5000, xkj_5000[0,0,:], 'r')
        plt.plot(x_CMG, xCH4_CMG, 'y')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.legend(('CMG-5000', 'Reference'))
        plt.grid()
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM.png')

        plt.figure(3)
        plt.plot(x, xkj[0,1,:])
        plt.grid()
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_y.png')

        plt.figure(4)
        plt.plot(x, xkj[2,0,:])
        plt.legend(('FOU', 'Reference'))
        plt.grid()
        plt.ylabel('Propane molar fraction in liquid phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_propane_x.png')

        plt.figure(5)
        plt.plot(x, xkj[2,0,:])
        plt.grid()
        plt.ylabel('Propane molar fraction in liquid phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_propane_y.png')


        plt.figure(6)
        plt.plot(x, P)
        plt.grid()
        plt.ylabel('Pressure')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_pressure.png')
        import pdb; pdb.set_trace()
