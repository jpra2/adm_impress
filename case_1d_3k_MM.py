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


        x_axis_xCH4_LLF_50 = np.array([0, 22.4628, 23.4607, 24.517, 25.5255, 26.5658, \
            27.569, 28.5722, 29.6178, 30.6316, 31.6295, 32.6805, 33.689, 34.6868, 50])
        xCH4_LLF_50 = np.array([0, 0, 0.326972, 0.3261, 0.324009, 0.322789, 0.321046, \
            0.320871, 0.292462, 0.228497, 0.17342, 0.127582, 0.0916776, 0.0650109, 0])

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_50_FOU_wa_2644.npy', allow_pickle=True)
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

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_512_FOU_wa_3206.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_512_FOU = data[6]
            Sg_512_FOU = data[7]
            z_512_FOU = data[10]
            xkj_512_FOU = data[13]
            P_512_FOU = data[4]

            xkj_512_FOU[:,1,Sg_512_FOU==0] = 0
            xkj_512_FOU[:,0,So_512_FOU==0] = 0
            n=512
            x_512 = np.linspace(0+50/(2*n),50-50/(2*n),n)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_5000_FOU_19254.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_5000_FOU = data[6]
            Sg_5000_FOU = data[7]
            z_5000_FOU = data[10]
            xkj_5000_FOU = data[13]
            P_5000_FOU = data[4]

            xkj_5000_FOU[:,1,Sg_5000_FOU==0] = 0
            xkj_5000_FOU[:,0,So_5000_FOU==0] = 0
            n=5000
            x_5000 = np.linspace(0+50/(2*n),50-50/(2*n),n)


        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_32_FOU_2184.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_32_FOU = data[6]
            Sg_32_FOU = data[7]
            z_32_FOU = data[10]
            xkj_32_FOU = data[13]
            P_32_FOU = data[4]

            xkj_32_FOU[:,1,Sg_32_FOU==0] = 0
            xkj_32_FOU[:,0,So_32_FOU==0] = 0
            n=32
            x_32 = np.linspace(0+50/(2*n),50-50/(2*n),n)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_50_LLF_wa_20097.npy', allow_pickle=True)
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

        '''datas = np.load('flying/results_case2_Moshiri_Manzari_3k_50_LLF2__3466.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_50_LLF = data[6]
            Sg_50_LLF = data[7]
            z_50_LLF = data[10]
            xkj_50_LLF = data[13]
            P_50_LLF = data[4]
            xkj_50_LLF[:,1,Sg_50_LLF==0] = 0
            xkj_50_LLF[:,0,So_50_LLF==0] = 0
            n=50
            x_50 = np.linspace(0+50/(2*n),50-50/(2*n),n)'''

        '''plt.figure(1)
        plt.plot(x_50, xkj_50_LLF[0,0,:], 'b')
        plt.plot(x_50, xkj_512_FOU[0,0,:], 'r')
        plt.plot(x_32, xkj_32_FOU[0,0,:], 'g')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk')
        plt.xlim(20,35)
        plt.grid()
        plt.legend(('LLF 50', 'FOU 512', 'Reference', 'Reference LLF 50'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_LLF.png')'''

        plt.figure(1)
        plt.plot(x_512, xkj_512_FOU[0,0,:], 'b')
        plt.plot(x_5000, xkj_5000_FOU[0,0,:], 'r')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.plot(x_CMG, xCH4_CMG, 'y')
        plt.xlim(20,35)
        plt.grid()
        plt.legend(('FOU 512', 'FOU 5000', 'Reference', 'CMG 5000'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_FOU.png')
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
