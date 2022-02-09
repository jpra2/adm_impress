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

        datas = np.load('flying/results_Hoteit_Firoo_3k_200_FR3_5175.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_200 = data[6]
            Sg_FR2_200 = data[7]
            oil_P_FR2 = data[8]
            z_FR2_200 = data[10]
            xkj_FR2_200 = data[13]
            P_FR2_200 = data[4]
            Nk_FR2_SPs = data[12]
            Nk_FR2 = data[12].sum(axis=-1)/2
            xkj_FR2_200[:,1,Sg_FR2_200==0] = 0
            xkj_FR2_200[:,0,So_FR2_200==0] = 0
            x_200 = np.linspace(0,50,len(So_FR2_200))

        datas = np.load('flying/results_Hoteit_Firoo_3k_200_FR3_7023.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR3_200 = data[6]
            Sg_FR3_200 = data[7]
            oil_P_FR3 = data[8]
            z_FR3_200 = data[10]
            xkj_FR3_200 = data[13]
            P_FR3_200 = data[4]
            Nk_FR3_SPs = data[12]
            w = np.array([1/3, 4/3, 1/3])
            Nk_FR3 = (data[12]*w[np.newaxis,np.newaxis,:]).sum(axis=-1)/2

            xkj_FR3_200[:,1,Sg_FR3_200<1e-14] = 0
            xkj_FR3_200[:,0,So_FR3_200<1e-14] = 0
            x_200 = np.linspace(0,50,len(So_FR3_200))

        datas = np.load('flying/results_Hoteit_Firoo_3k_200_FOU_4800.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_200 = data[6]
            Sg_200 = data[7]
            z_200 = data[10]
            xkj_200 = data[13]
            P_200 = data[4]

            xkj_200[:,1,Sg_200==0] = 0
            xkj_200[:,0,So_200==0] = 0
            x_200 = np.linspace(0,50,len(So_200))

        datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw_945.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_500 = data[6]
            Sg_500 = data[7]
            z_500 = data[10]
            xkj_500 = data[13]
            P_500 = data[4]

            xkj_500[:,1,Sg_500==0] = 0
            xkj_500[:,0,So_500==0] = 0
            x_500 = np.linspace(0,50,len(So_500))


        datas = np.load('flying/results_Hoteit_Firoo_3k_500_MUSCL_28322.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_500_MUSCL = data[6]
            Sg_500_MUSCL = data[7]
            z_500_MUSCL = data[10]
            xkj_500_MUSCL = data[13]
            P_500_MUSCL = data[4]

            xkj_500_MUSCL[:,1,Sg_500_MUSCL==0] = 0
            xkj_500_MUSCL[:,0,So_500_MUSCL==0] = 0
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

        datas = np.load('flying/results_Hoteit_Firoo_3k_128_IMPSAT_2383.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_128_IMPSAT = data[6]
            Sg_128_IMPSAT = data[7]
            z_128_IMPSAT = data[10]
            xkj_128_IMPSAT = data[13]
            P_128_IMPSAT = data[4]

            xkj_128_IMPSAT[:,1,Sg_128_IMPSAT<1e-13] = 0
            xkj_128_IMPSAT[:,0,So_128_IMPSAT<1e-13] = 0
            x_128 = np.linspace(0,50,len(So_128_IMPSAT))

        datas = np.load('flying/results_Hoteit_Firoo_3k_256_IMPSAT_3084.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_256_IMPSAT = data[6]
            Sg_256_IMPSAT = data[7]
            z_256_IMPSAT = data[10]
            xkj_256_IMPSAT = data[13]
            P_256_IMPSAT = data[4]

            xkj_256_IMPSAT[:,1,Sg_256_IMPSAT<1e-12] = 0
            xkj_256_IMPSAT[:,0,So_256_IMPSAT<1e-12] = 0
            x_256 = np.linspace(0,50,len(So_256_IMPSAT))

        datas = np.load('flying/results_Hoteit_Firoo_3k_512_IMPSAT_2507.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_512_IMPSAT = data[6]
            Sg_512_IMPSAT = data[7]
            z_512_IMPSAT = data[10]
            xkj_512_IMPSAT = data[13]
            P_512_IMPSAT = data[4]

            xkj_512_IMPSAT[:,1,Sg_512_IMPSAT<1e-12] = 0
            xkj_512_IMPSAT[:,0,So_512_IMPSAT<1e-12] = 0
            x_512 = np.linspace(0,50,len(So_512_IMPSAT))

        datas = np.load('flying/results_Hoteit_Firoo_3k_128_IMPEC_6868.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_128_IMPEC = data[6]
            Sg_128_IMPEC = data[7]
            z_128_IMPEC = data[10]
            xkj_128_IMPEC = data[13]
            P_128_IMPEC = data[4]

            xkj_128_IMPEC[:,1,Sg_128_IMPEC<1e-13] = 0
            xkj_128_IMPEC[:,0,So_128_IMPEC<1e-13] = 0
            x_128 = np.linspace(0,50,len(So_128_IMPEC))

        '''datas = np.load('flying/results_Hoteit_Firoo_3k_256_IMPEC_6112.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_256_IMPEC = data[6]
            Sg_256_IMPEC = data[7]
            z_256_IMPEC = data[10]
            xkj_256_IMPEC = data[13]
            P_256_IMPEC = data[4]

            xkj_256_IMPEC[:,1,Sg_256_IMPEC==0] = 0
            xkj_256_IMPEC[:,0,So_256_IMPEC==0] = 0
            x_256 = np.linspace(0,50,len(So_256_IMPEC))'''

        plt.figure(1)
        #plt.plot(x_500, xkj_500[2,1,:], 'b')
        plt.plot(x_500, xkj_500_MUSCL[2,1,:], 'm')
        plt.plot(x_axis_xC3, xC3, 'k')
        plt.plot(x_5000, xkj_5000[2,1,:], 'r')
        plt.plot(x_200, xkj_FR2_200[2,1,:], 'y')
        plt.plot(x_200, xkj_FR3_200[2,1,:], 'b')
        plt.plot(x_200, xkj_200[2,1,:], 'g')
        plt.legend(('MUSCL 500', 'Reference', 'FOU 5000', 'FR P1 200', 'FR P2 200', 'FOU-200'))
        plt.grid()
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_propane_y_HF.png')

        plt.figure(2)
        #plt.plot(x_500, xkj_500[2,1,:], 'b')
        plt.plot(x_500[200:400], xkj_500_MUSCL[2,1,200:400], 'm')
        plt.plot(x_axis_xC3, xC3, 'k')
        plt.plot(x_5000[2000:4000], xkj_5000[2,1,2000:4000], 'r')
        #plt.plot(x_200, xkj_FR2_200[2,1,:], 'y')
        #plt.plot(x_200, xkj_FR3_200[2,1,:], 'b')
        #plt.plot(x_200, xkj_200[2,1,:], 'g')
        plt.legend(('MUSCL 500', 'Reference', 'FOU 5000'))
        plt.grid()
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_propane_y_HF_MUSCL.png')

        plt.figure(3)
        #plt.plot(x_500, xkj_500[2,1,:], 'b')
        plt.plot(x_128, xkj_128_IMPSAT[2,1,:], 'm')
        #plt.plot(x_128, xkj_128_IMPEC[2,1,:], 'r')
        plt.plot(x_256, xkj_256_IMPSAT[2,1,:], 'b')
        plt.plot(x_512, xkj_512_IMPSAT[2,1,:], 'g')
        #plt.plot(x_axis_xC3, xC3, 'k')
        plt.plot(x_5000, xkj_5000[2,1,:], 'k')
        plt.xlim(20,30)
        plt.legend(('IMPSAT 128', 'IMPSAT 256', 'IMPSAT 512','FOU 5000'))
        plt.grid()
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_propane_y_HF_IMPSAT.png')
        import pdb; pdb.set_trace()
