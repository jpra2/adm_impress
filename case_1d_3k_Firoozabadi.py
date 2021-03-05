import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)

for  arq in arquivos:
    if  arq.startswith(name):
        '''-------------------------MUSCL LLF RESULTS------------------------'''
        x_axis_xCH40 = np.linspace(0, 23.075, 100)
        xCH40 = np.zeros(len(x_axis_xCH40))
        '''x_axis_xCH4 = np.array([22.3749, 22.4711, 22.6019, 22.767, 22.8144, 23.5318, 25.4384, 26.7638,
                        27.1177, 27.1588, 27.2937, 27.4269])
        xCH4 = np.array([0., 0.0582055, 0.216829, 0.34198, 0.346629, 0.346815, 0.343282, 0.336216,
                        0.122362, 0.0643422, 0.0126453, 0.])'''

        x_axis_xCH4 = np.array([23.075, 23.1167, 23.6458, 24.075, 24.425, 24.9292, 25.6833,  26.4958,
                        27.6583, 27.9125, 27.9833, 28.1333, 28.3042, 28.4958, 28.7542,  28.9958,   29.1958,
                        29.4792, 29.675, 29.9458, 30.3417])
        xCH4 = np.array([ 0., 0.333538, 0.332718, 0.331214, 0.329983, 0.326838, 0.322735, 0.320137,
                        0.32041, 0.217573, 0.171214, 0.12335, 0.089573, 0.0508718, 0.0265299, 0.0144957,
                        0.00861538, 0.00382906, 0.00259829, 0.00136752, 0.0])

        x_axis_xCH41 = np.linspace(x_axis_xCH4[-1], 50, 100)
        xCH41 = np.zeros(len(x_axis_xCH41))

        x_axis_xCH4 = np.concatenate((x_axis_xCH40,x_axis_xCH4))
        x_axis_xCH4 = np.concatenate((x_axis_xCH4,x_axis_xCH41))
        xCH4 = np.concatenate((xCH40,xCH4))
        xCH4 = np.concatenate((xCH4,xCH41))
        datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw_82170.npy', allow_pickle=True)
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
            plt.plot(x, xkj[0,0,:])
            plt.plot(x_axis_xCH4, xCH4)
            plt.legend(('FOU', 'Reference'))
            plt.grid()
            plt.ylabel('Methane molar fraction in liquid phase')
            plt.xlabel('Dimensionless distance')
            plt.savefig('results/compositional/3k_methane_x.png')

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
