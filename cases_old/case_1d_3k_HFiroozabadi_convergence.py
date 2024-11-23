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

        f = interp1d(x_axis_xC3,xC3)

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

        '-------------------------------- FOU ---------------------------------'
        datas = np.load('flying/results_Hoteit_Firoo_3k_128_FOU_5737.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_128 = data[6]
            Sg_FOU_128 = data[7]
            z_FOU_128 = data[10]
            xkj_FOU_128 = data[13]
            P_FOU_128 = data[4]

            xkj_FOU_128[:,1,Sg_FOU_128==0] = 0
            xkj_FOU_128[:,0,So_FOU_128==0] = 0
            x_128 = np.linspace(0,50,len(So_FOU_128))

        datas = np.load('flying/results_Hoteit_Firoo_3k_256_FOU_4720.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_256 = data[6]
            Sg_FOU_256 = data[7]
            z_FOU_256 = data[10]
            xkj_FOU_256 = data[13]
            P_FOU_256 = data[4]

            xkj_FOU_256[:,1,Sg_FOU_256==0] = 0
            xkj_FOU_256[:,0,So_FOU_256==0] = 0
            x_256 = np.linspace(0,50,len(So_FOU_256))

        datas = np.load('flying/results_Hoteit_Firoo_3k_512_FOU_3871.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_512 = data[6]
            Sg_FOU_512 = data[7]
            z_FOU_512 = data[10]
            xkj_FOU_512 = data[13]
            P_FOU_512 = data[4]

            xkj_FOU_512[:,1,Sg_FOU_512==0] = 0
            xkj_FOU_512[:,0,So_FOU_512==0] = 0
            x_512 = np.linspace(0,50,len(So_FOU_512))

        '-------------------------------- MUSCL -------------------------------'

        datas = np.load('flying/results_Hoteit_Firoo_3k_128_MUSCL_3181.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_128 = data[6]
            Sg_MUSCL_128 = data[7]
            z_MUSCL_128 = data[10]
            xkj_MUSCL_128 = data[13]
            P_MUSCL_128 = data[4]

            xkj_MUSCL_128[:,1,Sg_MUSCL_128==0] = 0
            xkj_MUSCL_128[:,0,So_MUSCL_128==0] = 0
            x_128 = np.linspace(0,50,len(So_MUSCL_128))

        datas = np.load('flying/results_Hoteit_Firoo_3k_256_MUSCL_2812.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_256 = data[6]
            Sg_MUSCL_256 = data[7]
            z_MUSCL_256 = data[10]
            xkj_MUSCL_256 = data[13]
            P_MUSCL_256 = data[4]

            xkj_MUSCL_256[:,1,Sg_MUSCL_256==0] = 0
            xkj_MUSCL_256[:,0,So_MUSCL_256==0] = 0
            x_256 = np.linspace(0,50,len(So_MUSCL_256))

        datas = np.load('flying/results_Hoteit_Firoo_3k_512_MUSCL_5491.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_512 = data[6]
            Sg_MUSCL_512 = data[7]
            z_MUSCL_512 = data[10]
            xkj_MUSCL_512 = data[13]
            P_MUSCL_512 = data[4]

            xkj_MUSCL_512[:,1,Sg_MUSCL_512==0] = 0
            xkj_MUSCL_512[:,0,So_MUSCL_512==0] = 0
            x_512 = np.linspace(0,50,len(So_MUSCL_512))

        '------------------------------- FR P1 --------------------------------'

        datas = np.load('flying/results_Hoteit_Firoo_3k_16_FR2_5500.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_16 = data[6]
            Sg_FR2_16 = data[7]
            z_FR2_16 = data[10]
            xkj_FR2_16 = data[13]
            P_FR2_16 = data[4]

            xkj_FR2_16[:,1,Sg_FR2_16==0] = 0
            xkj_FR2_16[:,0,So_FR2_16==0] = 0
            x_16 = np.linspace(0,50,len(So_FR2_16))

        datas = np.load('flying/results_Hoteit_Firoo_3k_32_FR2_3159.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_32 = data[6]
            Sg_FR2_32 = data[7]
            z_FR2_32 = data[10]
            xkj_FR2_32 = data[13]
            P_FR2_32 = data[4]

            xkj_FR2_32[:,1,Sg_FR2_32==0] = 0
            xkj_FR2_32[:,0,So_FR2_32==0] = 0
            x_32 = np.linspace(0,50,len(So_FR2_32))

        datas = np.load('flying/results_Hoteit_Firoo_3k_64_FR2_2482.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_64 = data[6]
            Sg_FR2_64 = data[7]
            z_FR2_64 = data[10]
            xkj_FR2_64 = data[13]
            P_FR2_64 = data[4]

            xkj_FR2_64[:,1,Sg_FR2_64==0] = 0
            xkj_FR2_64[:,0,So_FR2_64==0] = 0
            x_64 = np.linspace(0,50,len(So_FR2_64))

        datas = np.load('flying/results_Hoteit_Firoo_3k_128_FR2_2900.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_128 = data[6]
            Sg_FR2_128 = data[7]
            z_FR2_128 = data[10]
            xkj_FR2_128 = data[13]
            P_FR2_128 = data[4]

            xkj_FR2_128[:,1,Sg_FR2_128==0] = 0
            xkj_FR2_128[:,0,So_FR2_128==0] = 0
            x_128 = np.linspace(0,50,len(So_FR2_128))

        datas = np.load('flying/results_Hoteit_Firoo_3k_256_FR2_2413.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_256 = data[6]
            Sg_FR2_256 = data[7]
            z_FR2_256 = data[10]
            xkj_FR2_256 = data[13]
            P_FR2_256 = data[4]

            xkj_FR2_256[:,1,Sg_FR2_256==0] = 0
            xkj_FR2_256[:,0,So_FR2_256==0] = 0
            x_256 = np.linspace(0,50,len(So_FR2_256))

        datas = np.load('flying/results_Hoteit_Firoo_3k_512_FR2_4149.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_512 = data[6]
            Sg_FR2_512 = data[7]
            z_FR2_512 = data[10]
            xkj_FR2_512 = data[13]
            P_FR2_512 = data[4]

            xkj_FR2_512[:,1,Sg_FR2_512==0] = 0
            xkj_FR2_512[:,0,So_FR2_512==0] = 0
        '------------------------------- FR P2 --------------------------------'

        datas = np.load('flying/results_Hoteit_Firoo_3k_128_FR3_5026.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR3_128 = data[6]
            Sg_FR3_128 = data[7]
            z_FR3_128 = data[10]
            xkj_FR3_128 = data[13]
            P_FR3_128 = data[4]

            xkj_FR3_128[:,1,Sg_FR3_128==0] = 0
            xkj_FR3_128[:,0,So_FR3_128==0] = 0

        datas = np.load('flying/results_Hoteit_Firoo_3k_512_FR3_6445.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR3_512 = data[6]
            Sg_FR3_512 = data[7]
            z_FR3_512 = data[10]
            xkj_FR3_512 = data[13]
            P_FR3_512 = data[4]

            xkj_FR3_512[:,1,Sg_FR3_512==0] = 0
            xkj_FR3_512[:,0,So_FR3_512==0] = 0
            x_512 = np.linspace(0,50,512)


        plt.figure(1)
        plt.plot(x_5000,xkj_5000[2,1,:], 'k')
        plt.plot(x_128,xkj_FR2_128[2,1,:], 'b')
        plt.plot(x_128,xkj_MUSCL_128[2,1,:], 'g')
        plt.plot(x_128,xkj_FOU_128[2,1,:], 'r')
        #plt.plot(x_CMG, yC3H8_CMG, 'r')
        plt.grid()
        plt.title('Results for mesh: 128x1x1')
        plt.legend(('Reference Solution', 'FR P1', 'MUSCL', 'FOU'))
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_propane_y_HF_128.png')

        plt.figure(2)
        plt.plot(x_5000,xkj_5000[2,1,:], 'k')
        plt.plot(x_256,xkj_FR2_256[2,1,:], 'b')
        plt.plot(x_256,xkj_MUSCL_256[2,1,:], 'g')
        plt.plot(x_256,xkj_FOU_256[2,1,:], 'r')
        #plt.plot(x_CMG, yC3H8_CMG, 'r')
        plt.grid()
        plt.legend(('Reference Solution', 'FR P1', 'MUSCL', 'FOU'))
        plt.title('Results for mesh: 256x1x1')
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_propane_y_HF_256.png')

        plt.figure(3)
        plt.plot(x_5000,xkj_5000[2,1,:], 'k')
        plt.plot(x_512,xkj_FR2_512[2,1,:], 'b')
        plt.plot(x_512,xkj_MUSCL_512[2,1,:], 'g')
        plt.plot(x_512,xkj_FOU_512[2,1,:], 'r')
        #plt.plot(x_CMG, yC3H8_CMG, 'r')
        plt.grid()
        plt.title('Results for mesh: 512x1x1')
        plt.legend(('Reference Solution', 'FR P1', 'MUSCL', 'FOU'))
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_propane_y_HF_512.png')

        plt.figure(4)
        plt.plot(x_5000,xkj_5000[2,1,:], 'k')
        plt.plot(x_128,xkj_FR2_128[2,1,:], 'b')
        plt.plot(x_128,xkj_MUSCL_128[2,1,:], 'g')
        plt.plot(x_128,xkj_FOU_128[2,1,:], 'r')
        #plt.plot(x_CMG, yC3H8_CMG, 'r')
        plt.xlim(20,28)
        plt.ylim(.25, .35)
        plt.grid()
        plt.title('Results for mesh: 128x1x1')
        plt.legend(('Reference Solution', 'FR P1', 'MUSCL', 'FOU'))
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_propane_y_HF_128_zoom.png')

        plt.figure(5)
        plt.plot(x_5000,xkj_5000[2,1,:], 'k')
        plt.plot(x_256,xkj_FR2_256[2,1,:], 'b')
        plt.plot(x_256,xkj_MUSCL_256[2,1,:], 'g')
        plt.plot(x_256,xkj_FOU_256[2,1,:], 'r')
        #plt.plot(x_CMG, yC3H8_CMG, 'r')
        plt.xlim(20,28)
        plt.ylim(.25, .35)
        plt.grid()
        plt.legend(('Reference Solution', 'FR P1', 'MUSCL', 'FOU'))
        plt.title('Results for mesh: 256x1x1')
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_propane_y_HF_256_zoom.png')

        plt.figure(6)
        plt.plot(x_5000,xkj_5000[2,1,:], 'k')
        plt.plot(x_512,xkj_FR2_512[2,1,:], 'b')
        plt.plot(x_512,xkj_MUSCL_512[2,1,:], 'g')
        plt.plot(x_512,xkj_FOU_512[2,1,:], 'r')
        #plt.plot(x_CMG, yC3H8_CMG, 'r')
        plt.xlim(20,28)
        plt.ylim(.25,.35)
        plt.grid()
        plt.title('Results for mesh: 512x1x1')
        plt.legend(('Reference Solution', 'FR P1', 'MUSCL', 'FOU'))
        plt.ylabel('Propane molar fraction in gas phase')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/3k_propane_y_HF_512_zoom.png')
        import pdb; pdb.set_trace()
