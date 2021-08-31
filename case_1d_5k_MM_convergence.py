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

        x_MOC = np.array([0, 0.0615711, 0.121019, 0.181529, 0.214437, 0.214968, \
            0.396497, 0.544586, 0.725584, 0.726115, 0.835987, 0.835987, 1.0552, \
            1.05573, 1.09395, 1.09395, 1.5])
        Sg_MOC = np.array([1, 1, 1, 1, 1, 0.947941, 0.947941, 0.947941, 0.947552, \
        0.877622, 0.877622, 0.772727, 0.734266, 0.563326, 0.563326, 0, 0])

        f = interp1d(x_MOC,Sg_MOC)

        """---------------------- Convergence Study -------------------------"""

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FOU_85.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_8 = data[6]
            Sg_FOU_8 = data[7]
            z_FOU_8 = data[10]
            xkj_FOU_8 = data[13]

            xkj_FOU_8[:,1,Sg_FOU_8==0] = 0
            xkj_FOU_8[:,0,So_FOU_8==0] = 0
            x_8 = np.linspace(0,1.5,len(So_FOU_8))
            n = 8
            e8_L1_FOU = (sum(abs(f(x_8)-Sg_FOU_8))*(1/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FOU_159.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_16 = data[6]
            Sg_FOU_16 = data[7]
            z_FOU_16 = data[10]
            xkj_FOU_16 = data[13]

            xkj_FOU_16[:,1,Sg_FOU_16==0] = 0
            xkj_FOU_16[:,0,So_FOU_16==0] = 0
            x_16 = np.linspace(0,1.5,len(So_FOU_16))
            n = 16
            e16_L1_FOU = (sum(abs(f(x_16)-Sg_FOU_16))*(1/n))
            R16_L1_FOU = math.log(e8_L1_FOU/e16_L1_FOU,2)

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FOU_289.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_32 = data[6]
            Sg_FOU_32 = data[7]
            z_FOU_32 = data[10]
            xkj_FOU_32 = data[13]

            xkj_FOU_32[:,1,Sg_FOU_32==0] = 0
            xkj_FOU_32[:,0,So_FOU_32==0] = 0
            x_32 = np.linspace(0,1.5,len(So_FOU_32))
            n = 32
            e32_L1_FOU = (sum(abs(f(x_32)-Sg_FOU_32))*(1/n))
            R32_L1_FOU = math.log(e16_L1_FOU/e32_L1_FOU,2)

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FOU_518.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_64 = data[6]
            Sg_FOU_64 = data[7]
            z_FOU_64 = data[10]
            xkj_FOU_64 = data[13]

            xkj_FOU_64[:,1,Sg_FOU_64==0] = 0
            xkj_FOU_64[:,0,So_FOU_64==0] = 0
            x_64 = np.linspace(0,1.5,len(So_FOU_64))
            n = 64
            e64_L1_FOU = (sum(abs(f(x_64)-Sg_FOU_64))*(1/n))
            R64_L1_FOU = math.log(e32_L1_FOU/e64_L1_FOU,2)

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FOU_258.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_128 = data[6]
            Sg_FOU_128 = data[7]
            z_FOU_128 = data[10]
            xkj_FOU_128 = data[13]

            xkj_FOU_128[:,1,Sg_FOU_128==0] = 0
            xkj_FOU_128[:,0,So_FOU_128==0] = 0
            x_128 = np.linspace(0,1.5,len(So_FOU_128))
            n = 128
            e128_L1_FOU = (sum(abs(f(x_128)-Sg_FOU_128))*(1/n))
            R128_L1_FOU = math.log(e64_L1_FOU/e128_L1_FOU,2)

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FOU_498.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_256 = data[6]
            Sg_FOU_256 = data[7]
            z_FOU_256 = data[10]
            xkj_FOU_256 = data[13]

            xkj_FOU_256[:,1,Sg_FOU_256==0] = 0
            xkj_FOU_256[:,0,So_FOU_256==0] = 0
            x_256 = np.linspace(0,1.5,len(So_FOU_256))
            n = 256
            e256_L1_FOU = (sum(abs(f(x_256)-Sg_FOU_256))*(1/n))
            R256_L1_FOU = math.log(e128_L1_FOU/e256_L1_FOU,2)

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_FOU_976.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_512 = data[6]
            Sg_FOU_512 = data[7]
            z_FOU_512 = data[10]
            xkj_FOU_512 = data[13]

            xkj_FOU_512[:,1,Sg_FOU_512==0] = 0
            xkj_FOU_512[:,0,So_FOU_512==0] = 0
            x_512 = np.linspace(0,1.5,len(So_FOU_512))
            n = 512
            e512_L1_FOU = (sum(abs(f(x_512)-Sg_FOU_512))*(1/n))
            R512_L1_FOU = math.log(e256_L1_FOU/e512_L1_FOU,2)


        import pdb; pdb.set_trace()

        plt.figure(1)
        plt.plot(x_MOC, Sg_MOC, 'g')
        plt.plot(x_100, Sg_FOU_100, 'r')
        plt.plot(x_200, Sg_FOU_200, 'y')
        plt.plot(x_500, Sg_FOU_500, 'b')
        plt.legend(('MOC', 'FOU-100', 'FOU-200', 'FOU-500'))
        plt.grid()
        plt.ylabel('Sg')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/5k_Sg.png')

        plt.figure(2)
        #plt.plot(x_MOC, _MOC, 'g')
        plt.plot(x_100, z_FOU_100[1,:], 'r')
        plt.plot(x_200, z_FOU_200[1,:], 'y')
        plt.plot(x_500, z_FOU_500[1,:], 'b')
        plt.legend(('MOC', 'FOU-100', 'FOU-200', 'FOU-500'))
        plt.grid()
        plt.ylabel('Sg')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/5k_zC1.png')
