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

        #x_axis_xC31 = np.linspace(x_axis_xC3[-1], 50, 100)
        #xC31 = np.zeros(len(x_axis_xC31))

        '''x_axis_xC3 = np.concatenate((x_axis_xC30,x_axis_xC3))
        x_axis_xC3 = np.concatenate((x_axis_xC3,x_axis_xC31))
        xC3 = np.concatenate((xC30,xC3))
        xC3 = np.concatenate((xC3,xC31))'''


        f = interp1d(x_axis_xCH4,xCH4)

        """---------------------- Convergence Study -------------------------"""

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_8_FOU_8235.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_8 = data[6]
            Sg_FOU_8 = data[7]
            z_FOU_8 = data[10]
            xkj_FOU_8 = data[13]

            xkj_FOU_8[:,1,Sg_FOU_8==0] = 0
            xkj_FOU_8[:,0,So_FOU_8==0] = 0
            n = 8
            x_8 = np.linspace(0+50/(2*n),50-50/(2*n),len(So_FOU_8))

            e8_L1_FOU = (sum(abs(f(x_8)-Sg_FOU_8))*(1/n))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_16_FOU_25660.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_16 = data[6]
            Sg_FOU_16 = data[7]
            z_FOU_16 = data[10]
            xkj_FOU_16 = data[13]

            xkj_FOU_16[:,1,Sg_FOU_16==0] = 0
            xkj_FOU_16[:,0,So_FOU_16==0] = 0
            n = 16
            x_16 = np.linspace(0+50/(2*n),50-50/(2*n),len(So_FOU_16))

            e16_L1_FOU = (sum(abs(f(x_16)-Sg_FOU_16))*(1/n))
            R16_L1_FOU = math.log(e8_L1_FOU/e16_L1_FOU,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_32_FOU_5625.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_32 = data[6]
            Sg_FOU_32 = data[7]
            z_FOU_32 = data[10]
            xkj_FOU_32 = data[13]

            xkj_FOU_32[:,1,Sg_FOU_32==0] = 0
            xkj_FOU_32[:,0,So_FOU_32==0] = 0
            n = 32
            x_32 = np.linspace(0+50/(2*n),50-50/(2*n),len(So_FOU_32))

            e32_L1_FOU = (sum(abs(f(x_32)-Sg_FOU_32))*(1/n))
            R32_L1_FOU = math.log(e16_L1_FOU/e32_L1_FOU,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_64_FOU_5874.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_64 = data[6]
            Sg_FOU_64 = data[7]
            z_FOU_64 = data[10]
            xkj_FOU_64 = data[13]

            xkj_FOU_64[:,1,Sg_FOU_64==0] = 0
            xkj_FOU_64[:,0,So_FOU_64==0] = 0
            n = 64
            x_64 = np.linspace(0,50,len(So_FOU_64))

            e64_L1_FOU = (sum(abs(f(x_64)-Sg_FOU_64))*(1/n))
            R64_L1_FOU = math.log(e32_L1_FOU/e64_L1_FOU,2)

        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_128_FOU_4284.npy', allow_pickle=True)
        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_128_FOU_8590.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_128_FOU_6360.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_128 = data[6]
            Sg_FOU_128 = data[7]
            z_FOU_128 = data[10]
            xkj_FOU_128 = data[13]

            xkj_FOU_128[:,1,Sg_FOU_128==0] = 0
            xkj_FOU_128[:,0,So_FOU_128==0] = 0
            n = 128
            x_128 = np.linspace(0+50/(2*n),50-50/(2*n),len(So_FOU_128))

            e128_L1_FOU = (sum(abs(f(x_128)-Sg_FOU_128))*(1/n))
            R128_L1_FOU = math.log(e64_L1_FOU/e128_L1_FOU,2)

        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_256_FOU_3036.npy', allow_pickle=True)
        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_256_FOU_cc_7171.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_256_FOU_upw_2924.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_256_FOU_ha_4096.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_256_LLF_2792.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            So_FOU_256 = data[6]
            Sg_FOU_256 = data[7]
            z_FOU_256 = data[10]
            xkj_FOU_256 = data[13]

            xkj_FOU_256[:,1,Sg_FOU_256==0] = 0
            xkj_FOU_256[:,0,So_FOU_256==0] = 0
            x_256 = np.linspace(0+50/(2*n),50-50/(2*n),len(So_FOU_256))
            n = 256
            e256_L1_FOU = (sum(abs(f(x_256)-Sg_FOU_256))*(1/n))
            R256_L1_FOU = math.log(e128_L1_FOU/e256_L1_FOU,2)
            #import pdb; pdb.set_trace()

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_512_FOU_2642.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_512_FOU_5482.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_512_FOU_cc_5578.npy', allow_pickle=True)
        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_512_FOU_ha_3199.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_512_FOU_7615.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_512 = data[6]
            Sg_FOU_512 = data[7]
            z_FOU_512 = data[10]
            xkj_FOU_512 = data[13]

            xkj_FOU_512[:,1,Sg_FOU_512==0] = 0
            xkj_FOU_512[:,0,So_FOU_512==0] = 0
            n = 512
            x_512 = np.linspace(0+50/(2*n),50-50/(2*n),len(So_FOU_512))

            e512_L1_FOU = (sum(abs(f(x_512)-Sg_FOU_512))*(1/n))
            R512_L1_FOU = math.log(e256_L1_FOU/e512_L1_FOU,2)

        '''MUSCL'''

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_8_MUSCL_10209.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_8 = data[6]
            Sg_MUSCL_8 = data[7]
            z_MUSCL_8 = data[10]
            xkj_MUSCL_8 = data[13]

            xkj_MUSCL_8[:,1,Sg_MUSCL_8==0] = 0
            xkj_MUSCL_8[:,0,So_MUSCL_8==0] = 0
            n = 8
            e8_L1_MUSCL = (sum(abs(f(x_8)-Sg_MUSCL_8))*(1/n))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_16_MUSCL_4644.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_16 = data[6]
            Sg_MUSCL_16 = data[7]
            z_MUSCL_16 = data[10]
            xkj_MUSCL_16 = data[13]

            xkj_MUSCL_16[:,1,Sg_MUSCL_16==0] = 0
            xkj_MUSCL_16[:,0,So_MUSCL_16==0] = 0
            n = 16
            e16_L1_MUSCL = (sum(abs(f(x_16)-Sg_MUSCL_16))*(1/n))
            R16_L1_MUSCL = math.log(e8_L1_MUSCL/e16_L1_MUSCL,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_32_MUSCL_4983.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_32 = data[6]
            Sg_MUSCL_32 = data[7]
            z_MUSCL_32 = data[10]
            xkj_MUSCL_32 = data[13]

            xkj_MUSCL_32[:,1,Sg_MUSCL_32==0] = 0
            xkj_MUSCL_32[:,0,So_MUSCL_32==0] = 0
            n = 32
            e32_L1_MUSCL = (sum(abs(f(x_32)-Sg_MUSCL_32))*(1/n))
            R32_L1_MUSCL = math.log(e16_L1_MUSCL/e32_L1_MUSCL,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_64_MUSCL_3602.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_64 = data[6]
            Sg_MUSCL_64 = data[7]
            z_MUSCL_64 = data[10]
            xkj_MUSCL_64 = data[13]

            xkj_MUSCL_64[:,1,Sg_MUSCL_64==0] = 0
            xkj_MUSCL_64[:,0,So_MUSCL_64==0] = 0
            n = 64
            e64_L1_MUSCL = (sum(abs(f(x_64)-Sg_MUSCL_64))*(1/n))
            R64_L1_MUSCL = math.log(e32_L1_MUSCL/e64_L1_MUSCL,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_128_MUSCL_4449.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_128_MUSCL_8493.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_128 = data[6]
            Sg_MUSCL_128 = data[7]
            z_MUSCL_128 = data[10]
            xkj_MUSCL_128 = data[13]

            xkj_MUSCL_128[:,1,Sg_MUSCL_128==0] = 0
            xkj_MUSCL_128[:,0,So_MUSCL_128==0] = 0
            n = 128
            e128_L1_MUSCL = (sum(abs(f(x_128)-Sg_MUSCL_128))*(1/n))
            R128_L1_MUSCL = math.log(e64_L1_MUSCL/e128_L1_MUSCL,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_256_MUSCL_3880.npy', allow_pickle=True)
        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_256_MUSCL_cc_3674.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_256 = data[6]
            Sg_MUSCL_256 = data[7]
            z_MUSCL_256 = data[10]
            xkj_MUSCL_256 = data[13]

            xkj_MUSCL_256[:,1,Sg_MUSCL_256==0] = 0
            xkj_MUSCL_256[:,0,So_MUSCL_256==0] = 0
            n = 256
            e256_L1_MUSCL = (sum(abs(f(x_256)-Sg_MUSCL_256))*(1/n))
            R256_L1_MUSCL = math.log(e128_L1_MUSCL/e256_L1_MUSCL,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_512_MUSCL_7231.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_512 = data[6]
            Sg_MUSCL_512 = data[7]
            z_MUSCL_512 = data[10]
            xkj_MUSCL_512 = data[13]

            xkj_MUSCL_512[:,1,Sg_MUSCL_512==0] = 0
            xkj_MUSCL_512[:,0,So_MUSCL_512==0] = 0
            n = 512
            e512_L1_MUSCL = (sum(abs(f(x_512)-Sg_MUSCL_512))*(1/n))
            R512_L1_MUSCL = math.log(e256_L1_MUSCL/e512_L1_MUSCL,2)


        '''FR P1'''
        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_8_FR2_10332.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_8 = data[6]
            Sg_FR2_8 = data[7]
            z_FR2_8 = data[10]
            xkj_FR2_8 = data[13]

            xkj_FR2_8[:,1,Sg_FR2_8==0] = 0
            xkj_FR2_8[:,0,So_FR2_8==0] = 0
            n = 8
            e8_L1_FR2 = (sum(abs(f(x_8)-Sg_FR2_8))*(1/n))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_16_FR2_4869.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_16 = data[6]
            Sg_FR2_16 = data[7]
            z_FR2_16 = data[10]
            xkj_FR2_16 = data[13]

            xkj_FR2_16[:,1,Sg_FR2_16==0] = 0
            xkj_FR2_16[:,0,So_FR2_16==0] = 0
            n = 16
            e16_L1_FR2 = (sum(abs(f(x_16)-Sg_FR2_16))*(1/n))
            R16_L1_FR2 = math.log(e8_L1_FR2/e16_L1_FR2,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_32_FR2_5302.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_32 = data[6]
            Sg_FR2_32 = data[7]
            z_FR2_32 = data[10]
            xkj_FR2_32 = data[13]

            xkj_FR2_32[:,1,Sg_FR2_32==0] = 0
            xkj_FR2_32[:,0,So_FR2_32==0] = 0
            n = 32
            e32_L1_FR2 = (sum(abs(f(x_32)-Sg_FR2_32))*(1/n))
            R32_L1_FR2 = math.log(e16_L1_FR2/e32_L1_FR2,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_64_FR2_3791.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_64 = data[6]
            Sg_FR2_64 = data[7]
            z_FR2_64 = data[10]
            xkj_FR2_64 = data[13]

            xkj_FR2_64[:,1,Sg_FR2_64==0] = 0
            xkj_FR2_64[:,0,So_FR2_64==0] = 0
            n = 64
            e64_L1_FR2 = (sum(abs(f(x_64)-Sg_FR2_64))*(1/n))
            R64_L1_FR2 = math.log(e32_L1_FR2/e64_L1_FR2,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_128_FR2_4394.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_128 = data[6]
            Sg_FR2_128 = data[7]
            z_FR2_128 = data[10]
            xkj_FR2_128 = data[13]

            xkj_FR2_128[:,1,Sg_FR2_128==0] = 0
            xkj_FR2_128[:,0,So_FR2_128==0] = 0
            n = 128
            e128_L1_FR2 = (sum(abs(f(x_128)-Sg_FR2_128))*(1/n))
            R128_L1_FR2 = math.log(e64_L1_FR2/e128_L1_FR2,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_256_FR2_3526.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_256 = data[6]
            Sg_FR2_256 = data[7]
            z_FR2_256 = data[10]
            xkj_FR2_256 = data[13]

            xkj_FR2_256[:,1,Sg_FR2_256==0] = 0
            xkj_FR2_256[:,0,So_FR2_256==0] = 0
            n = 256
            e256_L1_FR2 = (sum(abs(f(x_256)-Sg_FR2_256))*(1/n))
            R256_L1_FR2 = math.log(e128_L1_FR2/e256_L1_FR2,2)

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_512_FR2_6012.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_512 = data[6]
            Sg_FR2_512 = data[7]
            z_FR2_512 = data[10]
            xkj_FR2_512 = data[13]

            xkj_FR2_512[:,1,Sg_FR2_512==0] = 0
            xkj_FR2_512[:,0,So_FR2_512==0] = 0
            n = 512
            e512_L1_FR2 = (sum(abs(f(x_512)-Sg_FR2_512))*(1/n))
            R512_L1_FR2 = math.log(e256_L1_FR2/e512_L1_FR2,2)


        plt.figure(1)
        plt.plot(x_8, xkj_FOU_8[0,0], 'c')
        plt.plot(x_16, xkj_FOU_16[0,0], 'g')
        plt.plot(x_32, xkj_FOU_32[0,0], 'r')
        plt.plot(x_64, xkj_FOU_64[0,0], 'y')
        plt.plot(x_128, xkj_FOU_128[0,0], 'c')
        plt.plot(x_256, xkj_FOU_256[0,0], 'm')
        plt.plot(x_512, xkj_FOU_512[0,0], 'b')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.xlim(20,30)
        plt.ylim(-0.05, 0.35)
        plt.grid()
        plt.legend(( 'FOU-8','FOU-16','FOU-32','FOU-64','FOU-128','FOU-256','FOU-512', 'Reference'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_FOU_conv.png')

        plt.figure(2)
        plt.plot(x_8, xkj_MUSCL_8[0,0], 'c')
        plt.plot(x_16, xkj_MUSCL_16[0,0], 'g')
        plt.plot(x_32, xkj_MUSCL_32[0,0], 'r')
        plt.plot(x_64, xkj_MUSCL_64[0,0], 'y')
        plt.plot(x_128, xkj_MUSCL_128[0,0], 'c')
        plt.plot(x_256, xkj_MUSCL_256[0,0], 'm')
        plt.plot(x_512, xkj_MUSCL_512[0,0], 'b')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.grid()
        plt.legend(( 'MUSCL-8','MUSCL-16','MUSCL-32','MUSCL-64','MUSCL-128','MUSCL-256','MUSCL-512', 'Reference'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_MUSCL_conv.png')

        plt.figure(3)
        plt.plot(x_8, xkj_FR2_8[0,0], 'c')
        plt.plot(x_16, xkj_FR2_16[0,0], 'g')
        plt.plot(x_32, xkj_FR2_32[0,0], 'r')
        plt.plot(x_64, xkj_FR2_64[0,0], 'y')
        plt.plot(x_128, xkj_FR2_128[0,0], 'c')
        plt.plot(x_256, xkj_FR2_256[0,0], 'm')
        plt.plot(x_512, xkj_FR2_512[0,0], 'b')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.grid()
        plt.legend(( 'FR2-8','FR2-16','FR2-32','FR2-64','FR2-128','FR2-256','FR2-512', 'Reference'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_FR2_conv.png')

        plt.figure(4)
        plt.plot(x_512, xkj_FOU_512[0,0,], 'm')
        plt.plot(x_512, xkj_FR2_512[0,0, ], 'b')
        plt.plot(x_512, xkj_MUSCL_512[0,0,], 'g')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.plot(x_CMG, xCH4_CMG, 'y')
        plt.xlim(20,30)
        plt.ylim(-0.05,0.35)
        plt.grid()
        plt.legend(( 'FOU-512','FR P1-512','MUSCL-512','Reference',
            'CMG 5000'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_512_comp.png')

        plt.figure(5)
        plt.plot(x_256, xkj_FOU_256[0,0, :], 'm')
        plt.plot(x_256, xkj_FR2_256[0,0, :], 'b')
        plt.plot(x_256, xkj_MUSCL_256[0,0,:], 'g')
        plt.plot(x_axis_xCH4[1:-1], xCH4[1:-1], 'k')
        plt.xlim(20,30)
        plt.ylim(-0.05,0.35)
        plt.grid()
        plt.legend(( 'FOU-256','FR P1-256','MUSCL-256','Reference'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_256_comp.png')

        plt.figure(6)
        plt.plot(x_128, xkj_FOU_128[0,0, :], 'm')
        plt.plot(x_128, xkj_FR2_128[0,0, :], 'b')
        plt.plot(x_128, xkj_MUSCL_128[0,0,:], 'g')
        plt.plot(x_axis_xCH4, xCH4, 'k')
        plt.grid()
        plt.xlim(20,30)
        plt.ylim(-0.05,0.35)
        plt.legend(( 'FOU-128','FR P1-128','MUSCL-128','Reference'))
        plt.ylabel('Methane molar fraction in liquid phase')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_128_comp.png')

        import pdb; pdb.set_trace()


        '''plt.figure(2)
        N = np.log10(np.array([8,16,32,64,128,256,512]))
        e_L1_FOU = np.log10(np.array([e8_L1_FOU, e16_L1_FOU, e32_L1_FOU, e64_L1_FOU, e128_L1_FOU, \
            e256_L1_FOU, e512_L1_FOU]))
        plt.plot(N, e_L1_FOU)
        plt.grid()
        plt.legend(( 'IMPEC-500', 'IMPSAT-500', 'Reference'))
        plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{10}({E}_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.title('Case2 - Moshiri and Manzari\'s paper')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/3k_methane_x_MM_logE.png')'''


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
