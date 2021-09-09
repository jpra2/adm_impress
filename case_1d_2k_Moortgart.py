import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
x_CMG = np.loadtxt('x_2k_Moort_CMG.txt')
Sg_CMG = np.loadtxt('Sg_2k_Moort_CMG.txt')
So_CMG = np.loadtxt('So_2k_Moort_CMG.txt')
yC3H8_CMG = np.loadtxt('yC3H8_2k_Moort_CMG.txt')
yC1H4_CMG = np.loadtxt('yCH4_2k_Moort_CMG.txt')
#xC1H4_CMG = np.loadtxt('xC3H8_2k_Moort_CMG.txt')
Sw_CMG = 1 - So_CMG - Sg_CMG
for  arq in arquivos:
    if  arq.startswith(name):
        '''-------------------------MUSCL LLF RESULTS------------------------'''
        x_So = np.array([0, 4.50351, 9.966788, 13.1414, 13.4367, 13.4736, 13.8058, 18.6416, 20.6349, 20.9671,
                        21.3363, 23.3296, 28.0177, 33.481, 38.0952, 42.8202, 52.9716, 63.2337,
                        68.8815, 73.4219, 78.0731, 88.3721, 98.6342])

        So_ans = np.array([0, 0, 0, 4.208754e-4, 0.18771, 0.450758, 0.614057, 0.613636, 0.614057, 0.664562,
                    0.825758, 0.824074, 0.824916, 0.824495, 0.824495, 0.824916, 0.824916, 0.824916,
                    0.824074, 0.825758, 0.825758, 0.824495, 0.826178])

        y_C3 = np.array([0, 0, 0, 0, 0.0238683, 0.0423868, 0.360288, 0.360082, 0, 0, 0, 0, 0])
        x_yC3 = np.array([0, 4.48087, 9.94536, 12.6776, 13.0783, 13.0783, 13.5155, 18.1785, 20.9107,
                        21.5665, 26.5209, 61.4208, 100])
        y_C1 = 1 - y_C3
        #import pdb; pdb.set_trace()
        x_Sw = np.array([0, 0.403374, 0.880088, 2.78695, 11.0744, 14.0814, 14.5581, 20.2054, 21.3788,
                        22.5156, 26.6227, 31.3531, 35.9736, 41.4375, 45.9846, 51.7052, 56.399,
                        62.0462, 66.7033, 71.6538, 76.8243, 81.4815, 88.8522, 95.8929, 98.3865])

        Sw_ans = np.array([0.0913008, 0.149675, 0.160732, 0.16065, 0.16047, 0.160732, 0.164065, 0.164065,
                            0.17252, 0.17374, 0.173659, 0.173659, 0.173659, 0.173659, 0.173659, 0.173659,
                            0.173659, 0.173659, 0.173659, 0.173659, 0.173659, 0.173659, 0.173659, 0.173659,
                            0.173659])

        x_zC1 = np.array([0, 0.842697, 1.68539, 5.55556, 10.0499, 13.764, 14.1698, 14.3571, 14.5755,
                        15.9176, 18.2896, 19.0699, 22.7216, 23.1586, 23.3146, 23.5643, 23.6267, 28.4644,
                        100])

        zC1_ans = np.array([0.422401, 0.243011, 0.227419, 0.228674, 0.229211, 0.229211, 0.216129, 0.190502,
                            0.171505, 0.164516, 0.16828, 0.16828, 0.167921, 0.15233, 0.0675627, 0.00824373,
                            0, 0, 0])

        datas = np.load('flying/results_Moortgat_3k_256_upw_70444.npy', allow_pickle=True)
        #datas = np.load('flying/results_Moortgat_3k_256_upw_69915.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_5000_upw_321408.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So256 = data[6]
            Sg256 = data[7]
            Sw256 = data[5]
            Nk = data[12]
            z256_w = Nk/np.sum(Nk,axis=0)
            z256 = data[10]
            xkj256 = data[13]
            P = data[4]

            xkj256[:,1,Sg256==0] = 0
            xkj256[:,0,So256==0] = 0
            x256 = np.linspace(0,100,len(So256))


        datas = np.load('flying/results_Moortgat_3k_512_upw_75736.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            So512 = data[6]
            Sg512 = data[7]
            Sw512 = data[5]
            z512 = data[10]
            Nk512 = data[12]
            z512_w = Nk512/np.sum(Nk512,axis=0)
            xkj512 = data[13]
            P = data[4]
            xkj512[:,1,Sg512==0] = 0
            xkj512[:,0,So512==0] = 0
            x512 = np.linspace(0,100,len(So512))

        '''datas = np.load('flying/results_Moortgat_3k_1024_upw_78928.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            So1024 = data[6]
            Sg1024 = data[7]
            Sw1024 = data[5]
            z1024 = data[10]
            Nk1024 = data[12]
            z1024_w = Nk1024/np.sum(Nk1024,axis=0)
            xkj1024 = data[13]
            P = data[4]

            xkj1024[:,1,Sg1024==0] = 0
            xkj1024[:,0,So1024==0] = 0
            x1024 = np.linspace(0,100,len(So1024))'''

        plt.figure(1)
        plt.plot(x256, So256, 'r', x512, So512, 'b')#, x1024, So1024, 'g')
        plt.plot(x_So, So_ans, 'k')
        plt.plot(x_CMG, So_CMG, 'y')
        plt.grid()
        plt.legend(('FOU - 256 elements', 'FOU - 512 elements',
            'Moortgart - 256 elements', 'CMG - 5000 elements'))
        plt.ylabel('Oil Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_Moortgart_oil_saturation.png')

        plt.figure(2)
        plt.plot(x256, Sw256, 'r', x512, Sw512, 'b')
        plt.plot(x_Sw, Sw_ans, 'k')
        plt.plot(x_CMG, Sw_CMG, 'y')
        plt.grid()
        plt.legend(('FOU - 256 elements', 'FOU - 512 elements', 'Moortgart - 256 elements', 'CMG - 5000 elements'))
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_Moortgart_water_saturation.png')

        plt.figure(3)
        plt.plot(x256, xkj256[1,1,:], 'r', x512, xkj512[1,1,:], 'b')
        plt.plot(x_yC3, y_C3, 'k')
        plt.plot(x_CMG, yC3H8_CMG, 'y')
        plt.grid()
        plt.legend(('FOU - 256 elements', 'FOU - 512 elements', 'Moortgart - 256 elements', 'CMG - 5000 elements'))
        plt.ylabel('Propane mole fraction in gas phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_Moortgart_yC3H8.png')

        plt.figure(5)
        plt.plot(x256, xkj256[0,1,:], 'r', x512, xkj512[0,1,:], 'b')
        #plt.plot(x_yC3, y_C3, 'k')
        #plt.plot(x_CMG, xC1H4_CMG, 'y')
        plt.grid()
        plt.legend(('FOU - 256 elements', 'FOU - 512 elements'))
        plt.ylabel('Propane mole fraction in gas phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_Moortgart_xCH4.png')

        plt.figure(6)

        plt.plot(x256, z256_w[0,:], 'y', x_zC1, zC1_ans, 'k', x512, z512_w[0,:])
        plt.grid()
        plt.legend(('FOU - 256 elements', 'Moortgat et. al. - 256 elements', 'FOU - 512 elements'))
        plt.ylabel('Overall molar fraction of C1')
        plt.xlabel('Distance [m]')
        plt.savefig('results/compositional/3k_Moortgat_zC1.png')

        plt.figure(4)
        plt.plot(x_So, So_ans, 'k')
        plt.plot(x256, So256, 'r')
        plt.plot(x_CMG, So_CMG, 'y')
        plt.grid()
        plt.legend(('Moortgat et. al. - 256 elements', 'FOU - 256 elements', 'CMG - 5000 elements'))
        plt.ylabel('Oil Saturation')
        plt.xlabel('Distance [m]')
        plt.savefig('results/compositional/Moortgart_Sun_Firooz_oil_saturation.png')

        plt.figure(5)
        plt.plot(x256, xkj256[0,1,:], 'r', x512, xkj512[0,1,:], 'b')
        plt.plot(x_yC3, y_C1, 'k')
        plt.plot(x_CMG, yC1H4_CMG, 'y')
        plt.grid()
        plt.legend(('FOU - 256 elements', 'FOU - 512 elements', 'Moortgat et. al. -256 elements'))
        plt.ylabel('Propane mole fraction in gas phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_Moortgart_yCH4.png')
        import pdb; pdb.set_trace()
