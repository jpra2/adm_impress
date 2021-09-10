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

        x_axis_xCH4 = np.array([23.1051, 23.1051, 23.1476, 23.9119, 24.0817, 24.4958, 24.9522,  25.7909,
                        26.8737, 27.7495, 27.7866, 27.8822, 27.9989, 28.0945, 28.259,  28.4289, 28.5881,
                        28.7686, 28.9756, 29.241, 29.6125, 29.8779, 30.483])
        xCH4 = np.array([ 1.743896e-4, 0.196886, 0.332735, 0.332735, 0.331166, 0.329248, 0.325237, 0.321749, 0.318784, 0.318261,
                        0.307623, 0.239437, 0.179098, 0.139686, 0.0953911, 0.0636522, 0.0430742, 0.0273792,
                        0.0162182, 0.00767314, 0.00348779, 0.00156951, 0.0])

        #x_axis_xC31 = np.linspace(x_axis_xC3[-1], 50, 100)
        #xC31 = np.zeros(len(x_axis_xC31))

        '''x_axis_xC3 = np.concatenate((x_axis_xC30,x_axis_xC3))
        x_axis_xC3 = np.concatenate((x_axis_xC3,x_axis_xC31))
        xC3 = np.concatenate((xC30,xC3))
        xC3 = np.concatenate((xC3,xC31))'''


        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_50_FOU_7767.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw_945.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_50 = data[6]
            Sg_50 = data[7]
            z_50 = data[10]
            xkj_50 = data[13]
            P_50 = data[4]

            xkj_50[:,1,Sg_50==0] = 0
            xkj_50[:,0,So_50==0] = 0
            x_50 = np.linspace(0,50,len(So_50))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FOU_2557.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw_945.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_200 = data[6]
            Sg_200 = data[7]
            z_200 = data[10]
            xkj_200 = data[13]
            P_200 = data[4]
            Nk_200 = data[12]

            xkj_200[:,1,Sg_200==0] = 0
            xkj_200[:,0,So_200==0] = 0
            x_200 = np.linspace(0,50,len(So_200))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_MUSCL_8208.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw_945.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_200 = data[6]
            Sg_MUSCL_200 = data[7]
            z_MUSCL_200 = data[10]
            xkj_MUSCL_200 = data[13]
            P_MUSCL_200 = data[4]

            xkj_MUSCL_200[:,1,Sg_MUSCL_200==0] = 0
            xkj_MUSCL_200[:,0,So_MUSCL_200==0] = 0
            x_MUSCL_200 = np.linspace(0,50,len(So_MUSCL_200))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_UPW_LLF_7176.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_LLF_200 = data[6]
            Sg_LLF_200 = data[7]
            z_LLF_200 = data[10]
            xkj_LLF_200 = data[13]
            P_LLF_200 = data[4]

            xkj_LLF_200[:,1,Sg_LLF_200==0] = 0
            xkj_LLF_200[:,0,So_LLF_200==0] = 0
            x_LLF_200 = np.linspace(0,50,len(So_LLF_200))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_UPW_LLF_ROE_7195.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_LLF_ROE_200 = data[6]
            Sg_LLF_ROE_200 = data[7]
            z_LLF_ROE_200 = data[10]
            xkj_LLF_ROE_200 = data[13]
            P_LLF_ROE_200 = data[4]

            xkj_LLF_ROE_200[:,1,Sg_LLF_ROE_200==0] = 0
            xkj_LLF_ROE_200[:,0,So_LLF_ROE_200==0] = 0
            x_LLF_ROE_200 = np.linspace(0,50,len(So_LLF_ROE_200))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_UPW_LLF_ROE_7710.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_LLF_ROE2_200 = data[6]
            Sg_LLF_ROE2_200 = data[7]
            z_LLF_ROE2_200 = data[10]
            xkj_LLF_ROE2_200 = data[13]
            P_LLF_ROE2_200 = data[4]

            xkj_LLF_ROE2_200[:,1,Sg_LLF_ROE2_200==0] = 0
            xkj_LLF_ROE2_200[:,0,So_LLF_ROE2_200==0] = 0
            x_LLF_ROE2_200 = np.linspace(0,50,len(So_LLF_ROE2_200))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FR2_7754.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_200 = data[6]
            Sg_FR2_200 = data[7]
            z_FR2_200 = data[10]
            xkj_FR2_200 = data[13]
            P_FR2_200 = data[4]
            Nk_FR2_200 = data[12].sum(axis=-1)/2

            xkj_FR2_200[:,1,Sg_FR2_200==0] = 0
            xkj_FR2_200[:,0,So_FR2_200==0] = 0
            x_FR2_200 = np.linspace(0,50,len(So_FR2_200))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FR3_9241.npy', allow_pickle=True) #20194
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FR3_9615.npy', allow_pickle=True) #20194
        for data in datas[datas.shape[0]-1:]:
            So_FR3_200 = data[6]
            Sg_FR3_200 = data[7]
            z_FR3_200 = data[10]
            xkj_FR3_200 = data[13]
            P_FR3_200 = data[4]

            xkj_FR3_200[:,1,Sg_FR3_200==0] = 0
            xkj_FR3_200[:,0,So_FR3_200==0] = 0

            Nk_FR3_200 = (data[12]*np.array([[[1/3,4/3,1/3]]])).sum(axis=-1)/2
            x_FR3_200 = np.linspace(0,50,len(So_FR3_200))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_500_FOU_4172.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_500 = data[6]
            Sg_500 = data[7]
            z_500 = data[10]
            xkj_500 = data[13]
            P_500 = data[4]

            xkj_500[:,1,Sg_500==0] = 0
            xkj_500[:,0,So_500==0] = 0
            x_500 = np.linspace(0,50,len(So_500))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_IMPSAT_2721.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_IMPSAT_200 = data[6]
            Sg_IMPSAT_200 = data[7]
            z_IMPSAT_200 = data[10]
            xkj_IMPSAT_200 = data[13]
            P_IMPSAT_200 = data[4]

            xkj_IMPSAT_200[:,1,Sg_IMPSAT_200==0] = 0
            xkj_IMPSAT_200[:,0,So_IMPSAT_200==0] = 0
            x_IMPSAT_200 = np.linspace(0,50,len(So_IMPSAT_200))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_500_IMPSAT_4278.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_IMPSAT_500 = data[6]
            Sg_IMPSAT_500 = data[7]
            z_IMPSAT_500 = data[10]
            xkj_IMPSAT_500 = data[13]
            P_IMPSAT_500 = data[4]

            xkj_IMPSAT_500[:,1,Sg_IMPSAT_500==0] = 0
            xkj_IMPSAT_500[:,0,So_IMPSAT_500==0] = 0
            x_IMPSAT_500 = np.linspace(0,50,len(So_IMPSAT_500))

        #datas = np.load('flying/results_Hoteit_Firoo_3k_5000_upw_7382.npy', allow_pickle=True)
        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_5000_IMPEC_19082.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_5000 = data[6]
            Sg_5000 = data[7]
            z_5000 = data[10]
            xkj_5000 = data[13]
            P_5000 = data[4]
            xkj_5000[:,1,Sg_5000==0] = 0
            xkj_5000[:,0,So_5000==0] = 0
            x_5000 = np.linspace(0,50,len(So_5000))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_5000_FOU_ROE_17147.npy', allow_pickle=True)
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
        #plt.plot(x_5000[2000:3500], xkj_5000[0,0,2000:3500], 'g')
        #plt.plot(x_500[200:350], xkj_500[0,0,200:350], 'b')
        plt.plot(x_FR2_200[80:140], xkj_FR2_200[0,0,80:140], 'y')
        plt.plot(x_FR3_200[80:140], xkj_FR3_200[0,0,80:140], 'g')
        plt.plot(x_200[80:140], xkj_200[0,0,80:140], 'r')
        plt.plot(x_MUSCL_200[80:140], xkj_MUSCL_200[0,0,80:140], 'b')
        #plt.plot(x_LLF_200[80:140], xkj_LLF_200[0,0,80:140], 'b')
        #plt.plot(x_LLF_ROE_200[80:140], xkj_LLF_ROE_200[0,0,80:140], 'g')
        #plt.plot(x_LLF_ROE2_200[80:140], xkj_LLF_ROE2_200[0,0,80:140], 'y')
        #plt.plot(x_50[20:35], xkj_50[0,0,20:35], 'm')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.grid()
        plt.legend(( 'FR P1-200', 'FR P2-203', 'FOU-200', 'MUSCL-200', 'Reference'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM.png')

        plt.figure(2)
        plt.plot(x_500[200:350], xkj_500[0,0,200:350], 'b')
        #plt.plot(x_5000[2000:3500], xkj_5000[0,0,2000:3500], 'g')
        #plt.plot(x_200[80:140], xkj_200[0,0,80:140], 'r')
        plt.plot(x_IMPSAT_500[200:350], xkj_IMPSAT_500[0,0,200:350], 'y')
        #plt.plot(x_IMPSAT_200[80:140], xkj_IMPSAT_200[0,0,80:140], 'm')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.grid()
        plt.legend(( 'IMPEC-500', 'IMPSAT-500', 'Reference'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_IMPSAT.png')

        plt.figure(3)
        plt.plot(x_5000[2000:3500], xkj_5000[0,0,2000:3500], 'c')
        #plt.plot(x_500[200:350], xkj_500[0,0,200:350], 'b')
        plt.plot(x_FR2_200[80:140], xkj_FR2_200[0,0,80:140], 'y')
        plt.plot(x_FR3_200[80:140], xkj_FR3_200[0,0,80:140], 'g')
        #plt.plot(x_MUSCL_200[80:140], xkj_MUSCL_200[0,0,80:140], 'b')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.grid()
        plt.legend(('FOU 5000', 'FR P1-200', 'FR P2-200', 'Reference'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_FR.png')

        plt.figure(4)
        plt.plot(x_FR2_200, P_FR2_200, 'y')
        plt.plot(x_FR3_200, P_FR3_200, 'g')
        plt.plot(x_200, P_200, 'b')
        plt.grid()
        plt.legend(( 'FR P1-200', 'FR P2-200', 'Reference'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_P.png')

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
