import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

        x_zC1_MOC = np.array([1, 1.0427, 1.04483, 1.04583, 1.04724, 1.04793, \
            1.04878, 1.05016, 1.05194, 1.05326, 1.0547, 1.05611, 1.05808, 1.06021, \
            1.06231, 1.06432, 1.06642, 1.06767, 1.0698, 1.07206, 1.07415, 1.07572, \
            1.07823, 1.08064, 1.08596, 1.0927, 1.09305, 1.09323, 1.09342, 1.09348, \
            1.09364, 1.09389, 1.09417, 1.09455, 1.09499, 1.09574, 1.09668, 1.09828, \
            1.09975, 1.10636, 1.12607, 1.13434, 1.14223, 1.14981])

        zC1_MOC = np.array([0, 2.978236e-4, 0.00476518, 0.00863688, 0.0178694, \
            0.0254639, 0.0379725, 0.070882, 0.133276, 0.20252, 0.282635, 0.353666, \
            0.43795, 0.506598, 0.551271, 0.57614, 0.593265, 0.600561, 0.608007, \
            0.612772, 0.615006, 0.616495, 0.617239, 0.617984, 0.617835, 0.617835, \
            0.612176, 0.60622, 0.593414, 0.561249, 0.502726, 0.384341, 0.297973, \
            0.285762, 0.275934, 0.264765, 0.257617, 0.251959, 0.250767, 0.249725, \
            0.249576, 0.249725, 0.249576, 0.250023])

        x_zC1_LLF = np.array([1.0021, 1.00968, 1.01732, 1.02494, 1.03258, 1.04016, \
            1.04759, 1.0552, 1.06288, 1.07021, 1.07794, 1.08515, 1.0933, 1.1005, \
            1.10824, 1.1157, 1.1234, 1.13105, 1.1386, 1.1459, 1.14972])
        zC1_LLF = np.array([0.0820504, 0.0972394, 0.114513, 0.133276, 0.152188, \
            0.172589, 0.193436, 0.215178, 0.237216, 0.258809, 0.280699, 0.301546, \
            0.322841, 0.333414, 0.310779, 0.249725, 0.216667, 0.222325, 0.22858, \
            0.233196, 0.234536])

        x_zC1_DW = np.array([1.00244, 1.00987, 1.01745, 1.025, 1.03261, 1.03997, \
            1.04774, 1.05523, 1.06275, 1.07046, 1.07804, 1.08543, 1.09295, 1.10063, \
            1.10827, 1.11275, 1.11563, 1.1234, 1.13073, 1.13847, 1.1459])
        zC1_DW = np.array([0.0479496, 0.0635853, 0.0817526, 0.103047, 0.129107, \
            0.158442, 0.188969, 0.223517, 0.259107, 0.295441, 0.33118, 0.366472, \
            0.397892, 0.39819, 0.304078, 0.268041, 0.241535, 0.241237, 0.241833, \
            0.241535, 0.243918])

        x_zC1_MDW = np.array([1.00244, 1.00987, 1.01745, 1.025, 1.03261, 1.03997, \
            1.04774, 1.05523, 1.06275, 1.07046, 1.07804, 1.08543, 1.09314, 1.10063, \
            1.10661, 1.10946, 1.11291, 1.11626, 1.12356, 1.13083, 1.13838])
        zC1_MDW = np.array([0.0479496, 0.0635853, 0.0817526, 0.103047, 0.129107, \
            0.158442, 0.188969, 0.223517, 0.259107, 0.295441, 0.33118, 0.366472, \
            0.39134, 0.379129, 0.335052, 0.303036, 0.271317, 0.237068, 0.237365, \
            0.238706, 0.241535])

        x_zC1_ROE = np.array([1.00244, 1.00987, 1.01745, 1.025, 1.03261, 1.03997, \
            1.04774, 1.05523, 1.06275, 1.07046, 1.07804, 1.08543, 1.09314, 1.10066, \
            1.10815, 1.11579, 1.12343, 1.13086, 1.13844, 1.14599])
        zC1_ROE = np.array([0.0479496, 0.0635853, 0.0817526, 0.103047, 0.129107, \
            0.158442, 0.188969, 0.223517, 0.259107, 0.295441, 0.33118, 0.366472, \
            0.397892, 0.419038, 0.347858, 0.235876, 0.237663, 0.238855, 0.240641, \
            0.242279])

        f = interp1d(x_MOC,Sg_MOC)

        """---------------------- Convergence Study -------------------------"""

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_4000_FOU_8726.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_4000 = data[6]
            Sg_FOU_4000 = data[7]
            zCO2_FOU_4000 = data[10][0]
            x_4000 = np.linspace(0,1.5,4000)
            f = interp1d(x_4000,zCO2_FOU_4000)

        """----------------------------- FOU --------------------------------"""

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FOU_37.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_8 = data[6]
            Sg_FOU_8 = data[7]
            zC1_FOU_8 = data[10][1]
            zCO2_FOU_8 = data[10][0]
            xkj_FOU_8 = data[13]

            xkj_FOU_8[:,1,Sg_FOU_8==0] = 0
            xkj_FOU_8[:,0,So_FOU_8==0] = 0
            x_8 = np.linspace(0,1.5,len(So_FOU_8))
            n = 8
            e8_L1_FOU = (sum(abs(f(x_8)- data[10][0]))*(1/n))
            t8_FOU = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FOU_36.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_16 = data[6]
            Sg_FOU_16 = data[7]
            zC1_FOU_16 = data[10][1]
            zCO2_FOU_16 = data[10][0]
            xkj_FOU_16 = data[13]

            xkj_FOU_16[:,1,Sg_FOU_16==0] = 0
            xkj_FOU_16[:,0,So_FOU_16==0] = 0
            x_16 = np.linspace(0,1.5,len(So_FOU_16))
            n = 16
            e16_L1_FOU = (sum(abs(f(x_16)-data[10][0]))*(1/n))
            R16_L1_FOU = math.log(e8_L1_FOU/e16_L1_FOU,2)
            t16_FOU = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FOU_118.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_32 = data[6]
            Sg_FOU_32 = data[7]
            zC1_FOU_32 = data[10][1]
            zCO2_FOU_32 = data[10][0]
            xkj_FOU_32 = data[13]

            xkj_FOU_32[:,1,Sg_FOU_32==0] = 0
            xkj_FOU_32[:,0,So_FOU_32==0] = 0
            x_32 = np.linspace(0,1.5,len(So_FOU_32))
            n = 32
            e32_L1_FOU = (sum(abs(f(x_32)-data[10][0]))*(1/n))
            R32_L1_FOU = math.log(e16_L1_FOU/e32_L1_FOU,2)
            e32_L2_FOU = np.sqrt(np.sum((f(x_32)-data[10][0])**2) * 1 / n)

            t32_FOU = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FOU_142.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_64 = data[6]
            Sg_FOU_64 = data[7]
            zC1_FOU_64 = data[10][1]
            zCO2_FOU_64 = data[10][0]
            xkj_FOU_64 = data[13]

            xkj_FOU_64[:,1,Sg_FOU_64==0] = 0
            xkj_FOU_64[:,0,So_FOU_64==0] = 0
            x_64 = np.linspace(0,1.5,len(So_FOU_64))
            n = 64
            e64_L1_FOU = (sum(abs(f(x_64)-data[10][0]))*(1/n))
            R64_L1_FOU = math.log(e32_L1_FOU/e64_L1_FOU,2)
            e64_L2_FOU = np.sqrt(np.sum((f(x_64)-data[10][0])**2) * 1 / n)

            t64_FOU = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FOU_284.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_128 = data[6]
            Sg_FOU_128 = data[7]
            zC1_FOU_128 = data[10][1]
            zCO2_FOU_128 = data[10][0]
            xkj_FOU_128 = data[13]

            xkj_FOU_128[:,1,Sg_FOU_128==0] = 0
            xkj_FOU_128[:,0,So_FOU_128==0] = 0
            x_128 = np.linspace(0,1.5,len(So_FOU_128))
            n = 128
            e128_L1_FOU = (sum(abs(f(x_128)-data[10][0]))*(1/n))
            R128_L1_FOU = math.log(e64_L1_FOU/e128_L1_FOU,2)
            e128_L2_FOU = np.sqrt(np.sum((f(x_128)-data[10][0])**2) * 1 / n)

            t128_FOU = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_200_FOU_1125.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_200 = data[6]
            Sg_FOU_200 = data[7]
            zC1_FOU_200 = data[10][1]
            xkj_FOU_200 = data[13]

            xkj_FOU_200[:,1,Sg_FOU_200==0] = 0
            xkj_FOU_200[:,0,So_FOU_200==0] = 0
            x_200 = np.linspace(0,1.5,len(So_FOU_200))


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FOU_579.npy', allow_pickle=True)
        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FOU_t_404.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_256 = data[6]
            Sg_FOU_256 = data[7]
            zC1_FOU_256 = data[10][1]
            zCO2_FOU_256 = data[10][0]
            xkj_FOU_256 = data[13]

            xkj_FOU_256[:,1,Sg_FOU_256==0] = 0
            xkj_FOU_256[:,0,So_FOU_256==0] = 0
            x_256 = np.linspace(0,1.5,len(So_FOU_256))
            n = 256
            e256_L1_FOU = (sum(abs(f(x_256)-data[10][0]))*(1/n))
            R256_L1_FOU = math.log(e128_L1_FOU/e256_L1_FOU,2)
            e256_L2_FOU = np.sqrt(np.sum((f(x_256)-data[10][0])**2) * 1 / n)

            t256_FOU = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_FOU_976.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_512 = data[6]
            Sg_FOU_512 = data[7]
            zC1_FOU_512 = data[10][1]
            zCO2_FOU_512 = data[10][0]
            xkj_FOU_512 = data[13]

            xkj_FOU_512[:,1,Sg_FOU_512==0] = 0
            xkj_FOU_512[:,0,So_FOU_512==0] = 0
            x_512 = np.linspace(0,1.5,len(So_FOU_512))
            n = 512
            e512_L1_FOU = (sum(abs(f(x_512)-data[10][0]))*(1/n))
            R512_L1_FOU = math.log(e256_L1_FOU/e512_L1_FOU,2)
            e512_L2_FOU = np.sqrt(np.sum((f(x_512)-data[10][0])**2) * 1 / n)

            t512_FOU = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_1024_FOU_2288.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_1024 = data[6]
            Sg_FOU_1024 = data[7]
            zC1_FOU_1024 = data[10][1]
            zCO2_FOU_1024 = data[10][0]

            x_1024 = np.linspace(0,1.5,len(So_FOU_1024))
            n = 1024
            e1024_L1_FOU = (sum(abs(f(x_1024)-data[10][0]))*(1/n))
            R1024_L1_FOU = math.log(e512_L1_FOU/e1024_L1_FOU,2)
            e1024_L2_FOU = np.sqrt(np.sum((f(x_1024)-data[10][0])**2) * 1 / n)

            t1024_FOU = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_2048_FOU_5139.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_2048 = data[6]
            Sg_FOU_2048 = data[7]
            zC1_FOU_2048 = data[10][1]
            zCO2_FOU_2048 = data[10][0]

            x_2048 = np.linspace(0,1.5,len(So_FOU_2048))
            n = 2048
            e2048_L1_FOU = (sum(abs(f(x_2048)-data[10][0]))*(1/n))
            e2048_L2_FOU = np.sqrt(np.sum((f(x_2048)-data[10][0])**2) * 1 / n)

            R2048_L1_FOU = math.log(e1024_L1_FOU/e2048_L1_FOU,2)
            t2048_FOU = data[2]

        """----------------------------- MUSCL --------------------------------"""

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_MUSCL_LLF_89.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_8 = data[7]
            zC1_MUSCL_8 = data[10][1]
            zCO2_MUSCL_8 = data[10][0]
            n = 8
            e8_L1_MUSCL = (sum(abs(f(x_8)-data[10][0]))*(1/n))
            t8_MUSCL = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_MUSCL_108.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_16 = data[7]
            zC1_MUSCL_16 = data[10][1]
            n = 16
            e16_L1_MUSCL = (sum(abs(f(x_16)-data[10][0]))*(1/n))
            R16_L1_MUSCL = math.log(e8_L1_MUSCL/e16_L1_MUSCL,2)
            t16_MUSCL = data[2]

        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_MUSCL_149.npy', allow_pickle=True)
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_MUSCL_MDW_197.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_32 = data[7]
            zC1_MUSCL_32 = data[10][1]
            zCO2_MUSCL_32 = data[10][0]
            n = 32
            e32_L1_MUSCL = (sum(abs(f(x_32)-data[10][0]))*(1/n))
            R32_L1_MUSCL = math.log(e16_L1_MUSCL/e32_L1_MUSCL,2)
            e32_L2_MUSCL = np.sqrt(np.sum((f(x_32)-data[10][0])**2) * 1 / n)

            t32_MUSCL = data[2]

        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_MUSCL_352.npy', allow_pickle=True)
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_MUSCL_MDW_438.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_64 = data[7]
            zC1_MUSCL_64 = data[10][1]
            zCO2_MUSCL_64 = data[10][0]
            n = 64
            e64_L1_MUSCL = (sum(abs(f(x_64)-data[10][0]))*(1/n))
            R64_L1_MUSCL = math.log(e32_L1_MUSCL/e64_L1_MUSCL,2)
            e64_L2_MUSCL = np.sqrt(np.sum((f(x_64)-data[10][0])**2) * 1 / n)

            t64_MUSCL = data[2]

        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_MUSCL_LLF_994.npy', allow_pickle=True)
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_MUSCL_MDW_922.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_128 = data[7]
            zC1_MUSCL_128 = data[10][1]
            zCO2_MUSCL_128 = data[10][0]
            n = 128
            e128_L1_MUSCL = (sum(abs(f(x_128)-data[10][0]))*(1/n))
            R128_L1_MUSCL = math.log(e64_L1_MUSCL/e128_L1_MUSCL,2)
            e128_L2_MUSCL = np.sqrt(np.sum((f(x_128)-data[10][0])**2) * 1 / n)

            t128_MUSCL = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_200_MUSCL_LLF_t_1609.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_200 = data[6]
            Sg_MUSCL_200 = data[7]
            z_MUSCL_200 = data[10]
            xkj_MUSCL_200 = data[13]

            xkj_MUSCL_200[:,1,Sg_MUSCL_200==0] = 0
            xkj_MUSCL_200[:,0,So_MUSCL_200==0] = 0
            t_MUSCL_200 = data[2]

        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_MUSCL_LLF_t_2044.npy', allow_pickle=True)
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_MUSCL_MDW_1832.npy', allow_pickle=True)
        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_MUSCL_LLF_t_1434.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_256 = data[7]
            zC1_MUSCL_256 = data[10][1]
            zCO2_MUSCL_256 = data[10][0]
            n = 256
            e256_L1_MUSCL = (sum(abs(f(x_256)-data[10][0]))*(1/n))
            e256_L2_MUSCL = np.sqrt(np.sum((f(x_256)-data[10][0])**2) * 1 / n)
            R256_L1_MUSCL = math.log(e128_L1_MUSCL/e256_L1_MUSCL,2)
            t256_MUSCL = data[2]

        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_MUSCL_LLF_9198.npy', allow_pickle=True)
        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_MUSCL_LLF_4260.npy', allow_pickle=True)
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_MUSCL_MDW_3662.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_512 = data[7]
            zC1_MUSCL_512 = data[10][1]
            zCO2_MUSCL_512 = data[10][0]
            n = 512
            e512_L2_MUSCL = np.sqrt(np.sum((f(x_512)-data[10][0])**2) * 1 / n)
            e512_L1_MUSCL = (sum(abs(f(x_512)-data[10][0]))*(1/n))
            R512_L1_MUSCL = math.log(e256_L1_MUSCL/e512_L1_MUSCL,2)
            t512_MUSCL = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_1024_MUSCL_7654.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_1024 = data[7]
            zC1_MUSCL_1024 = data[10][1]
            zCO2_MUSCL_1024 = data[10][0]
            n = 1024
            x_1024 = np.linspace(0, 1.5, 1024)
            e1024_L1_MUSCL = (sum(abs(f(x_1024)-data[10][0]))*(1/n))
            R1024_L1_MUSCL = math.log(e512_L1_MUSCL/e1024_L1_MUSCL,2)
            t1024_MUSCL = data[2]

        '''datas = np.load('flying/results_case1_Moshiri_Manzari_5k_2048_LLF_7654.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_2048 = data[7]
            zC1_MUSCL_2048 = data[10][1]
            zCO2_MUSCL_2048 = data[10][0]
            n = 2048
            e2048_L1_MUSCL = (sum(abs(f(x_2048)-data[10][0]))*(1/n))
            #R2048_L1_MUSCL = math.log(e512_L1_MUSCL/e2048_L1_MUSCL,2)
            t2048_MUSCL = data[2]'''

        size = 5
        plt.figure(1)
        plt.plot(x_MOC, Sg_MOC, 'k')
        #plt.plot(x_8, Sg_MUSCL_8, '-r', mfc='none')
        #plt.plot(x_16, Sg_MUSCL_16, '-y', mfc='none')
        plt.plot(x_32, Sg_MUSCL_32, '-go', mfc='none', markersize=size)
        plt.plot(x_64, Sg_MUSCL_64, '-bp', mfc='none', markersize=size)
        plt.plot(x_128, Sg_MUSCL_128, '-mD', mfc='none', markersize=size)
        #plt.plot(x_200, Sg_MUSCL_200, 'y', mfc='none', markersize=size)
        plt.plot(x_256, Sg_MUSCL_256, '-cs', mfc='none', markersize=size)
        #plt.plot(x_512, Sg_MUSCL_512, '-rv', mfc='none', markersize=size)
        #plt.plot(x_1024, Sg_MUSCL_1024, '-y')
        plt.legend(('MOC', 'MUSCL-32', 'MUSCL-64', \
            'MUSCL-128', 'MUSCL-256'))
        plt.grid()
        plt.ylabel('Sg')
        plt.xlabel('Distância')
        plt.savefig('results/compositional/TCC2/5k_Sg_MUSCL.png')

        plt.figure(2)
        plt.plot(x_MOC, Sg_MOC, 'k')
        #plt.plot(x_8, Sg_FOU_8, '-r', mfc='none')
        #plt.plot(x_16, Sg_FOU_16, '-y', mfc='none')
        plt.plot(x_32, Sg_FOU_32, '-go', mfc='none', markersize=size)
        plt.plot(x_64, Sg_FOU_64, '-bp', mfc='none', markersize=size)
        plt.plot(x_128, Sg_FOU_128, '-mD', mfc='none', markersize=size)
        plt.plot(x_256, Sg_FOU_256, '-cs', mfc='none', markersize=size)
        #plt.plot(x_512, Sg_FOU_512, '-rv', mfc='none', markersize=size)
        plt.legend(('MOC', 'FOU-32', 'FOU-64', \
            'FOU-128', 'FOU-256'))
        plt.grid()
        plt.ylabel('Sg')
        plt.xlabel('Distância')
        plt.savefig('results/compositional/TCC2/5k_Sg_FOU.png')

        plt.figure(3)
        plt.plot(x_MOC, Sg_MOC, 'k')
        #plt.plot(x_8, Sg_MUSCL_8, '-r', mfc='none')
        #plt.plot(x_16, Sg_MUSCL_16, '-y', mfc='none')
        plt.plot(x_32, Sg_MUSCL_32, '-go', mfc='none', markersize=size)
        plt.plot(x_64, Sg_MUSCL_64, '-bp', mfc='none', markersize=size)
        plt.plot(x_128, Sg_MUSCL_128, '-mD', mfc='none', markersize=size)
        #plt.plot(x_200, Sg_MUSCL_200, 'ys', mfc='none', markersize=size)
        plt.plot(x_256, Sg_MUSCL_256, '-c<', mfc='none', markersize=size)
        #plt.plot(x_512, Sg_MUSCL_512, '-r*', mfc='none', markersize=6)
        plt.legend(('MOC', 'MUSCL-32', 'MUSCL-64', \
            'MUSCL-128', 'MUSCL-256'))
        plt.grid()
        plt.xlim((0.5, 1.2))
        plt.ylim((0.4, 1))
        plt.ylabel('Sg')
        plt.xlabel('Distância')
        plt.savefig('results/compositional/TCC2/5k_Sg_MUSCL_zoom.png')

        plt.figure(4)
        plt.plot(x_MOC, Sg_MOC, 'k')
        #plt.plot(x_8, Sg_FOU_8, '-r', mfc='none')
        #plt.plot(x_16, Sg_FOU_16, '-y', mfc='none')
        plt.plot(x_32, Sg_FOU_32, '-go', mfc='none', markersize=size)
        plt.plot(x_64, Sg_FOU_64, '-bp', mfc='none', markersize=size)
        plt.plot(x_128, Sg_FOU_128, '-mD', mfc='none', markersize=size)
        plt.plot(x_256, Sg_FOU_256, '-cs', mfc='none', markersize=size)
        #plt.plot(x_512, Sg_FOU_512, '-rv', mfc='none', markersize=size)
        plt.legend(('MOC', 'FOU-32', 'FOU-64', \
            'FOU-128', 'FOU-256'))
        plt.grid()
        plt.xlim((0.5, 1.2))
        plt.ylim((0.4, 1))
        plt.ylabel('Sg')
        plt.xlabel('Distância')
        plt.savefig('results/compositional/TCC2/5k_Sg_FOU_zoom.png')

        plt.figure(5)
        plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        #plt.plot(x_8, Sg_MUSCL_8, '-r', mfc='none')
        #plt.plot(x_16, Sg_MUSCL_16, '-y', mfc='none')
        plt.plot(x_32, zC1_MUSCL_32, '-go', mfc='none', markersize=size)
        plt.plot(x_64, zC1_MUSCL_64, '-bp', mfc='none', markersize=size)
        plt.plot(x_128, zC1_MUSCL_128, '-mD', mfc='none', markersize=size)
        plt.plot(x_256, zC1_MUSCL_256, '-cs', mfc='none', markersize=size)
        plt.plot(x_512, zC1_MUSCL_512, '-rv', mfc='none', markersize=size)
        plt.legend(('MOC', 'MUSCL-32', 'MUSCL-64', \
            'MUSCL-128', 'MUSCL-256', 'MUSCL-512'))
        plt.grid()
        plt.xlim((1, 1.14))
        plt.ylim((0, 0.65))
        plt.ylabel('$z_{C_1}$')
        plt.xlabel('Distância')
        plt.savefig('results/compositional/TCC2/5k_zC1_MUSCL.png')

        plt.figure(6)
        plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        #plt.plot(x_8, Sg_FOU_8, '-r', mfc='none')
        #plt.plot(x_16, Sg_FOU_16, '-y', mfc='none')
        plt.plot(x_32, zC1_FOU_32, '-go', mfc='none', markersize=size)
        plt.plot(x_64, zC1_FOU_64, '-bp', mfc='none', markersize=size)
        plt.plot(x_128, zC1_FOU_128, '-mD', mfc='none', markersize=size)
        plt.plot(x_256, zC1_FOU_256, '-cs', mfc='none', markersize=size)
        plt.plot(x_512, zC1_FOU_512, '-rv', mfc='none', markersize=size)
        plt.legend(('MOC', 'FOU-32', 'FOU-64', \
            'FOU-128', 'FOU-256', 'FOU-512'))
        plt.grid()
        plt.xlim((1, 1.14))
        plt.ylim((0, 0.65))
        plt.ylabel('$x_{C_1}')
        plt.xlabel('Distância')
        plt.savefig('results/compositional/TCC2/5k_zC1_FOU.png')

        plt.figure(7)
        x = np.log10(np.array([32,64,128,256,512]))
        eL1_FOU = np.log10(np.array([e32_L1_FOU, e64_L1_FOU, \
            e128_L1_FOU, e256_L1_FOU, e512_L1_FOU]))
        eL1_MUSCL = np.log10(np.array([e32_L1_MUSCL, e64_L1_MUSCL, \
            e128_L1_MUSCL, e256_L1_MUSCL, e512_L1_MUSCL]))
        y = -x
        plt.plot(x, eL1_FOU, '-ro', mfc='none', markersize=size)
        plt.plot(x, eL1_MUSCL, '-gs', mfc='none', markersize=size)
        plt.plot(x,y,'k')
        plt.legend(('FOU', 'MUSCL', 'Primeira Ordem'))
        plt.grid()
        plt.ylabel('$log(E_{L1})$')
        plt.xlabel('log($n_b$)')
        plt.savefig('results/compositional/TCC2/5k_eL1_zCO2.png')

        plt.figure(8)
        x = (np.array([32,64,128,256,512]))
        t_FOU = (np.array([t32_FOU, t64_FOU, \
            t128_FOU, t256_FOU, t512_FOU]))
        t_MUSCL = (np.array([t32_MUSCL, t64_MUSCL, \
            t128_MUSCL, t256_MUSCL, t512_MUSCL]))
        plt.plot(x, t_FOU, '-ro', mfc='none', markersize=size)
        plt.plot(x, t_MUSCL, '-gs', mfc='none', markersize=size)
        plt.legend(('FOU', 'MUSCL', 'Primeira Ordem'))
        plt.grid()
        plt.ylabel('Tempo computacional [s]')
        plt.xlabel('Número de volumes de controle')
        plt.savefig('results/compositional/TCC2/5k_time.png')

        plt.figure(9)
        x = (np.array([32,64,128,256,512]))
        t_FOU = (np.array([t32_FOU, t64_FOU, \
            t128_FOU, t256_FOU, t512_FOU, t1024_FOU]))
        t_MUSCL = (np.array([t32_MUSCL, t64_MUSCL, \
            t128_MUSCL, t256_MUSCL, t512_MUSCL]))
        eL1_FOU = (np.array([e32_L1_FOU, e64_L1_FOU, \
            e128_L1_FOU, e256_L1_FOU, e512_L1_FOU, e1024_L1_FOU]))
        eL1_MUSCL = (np.array([e32_L1_MUSCL, e64_L1_MUSCL, \
            e128_L1_MUSCL, e256_L1_MUSCL, e512_L1_MUSCL]))
        size_FOU = 7
        size_MUSCL = 7
        #plt.plot(t_FOU,eL1_FOU, '-ro', mfc='none', markersize=size)
        #plt.plot(t_MUSCL,eL1_MUSCL, '-gs', mfc='none', markersize=size)
        plt.plot(t32_FOU, e32_L1_FOU, '-ro', mfc='none', markersize=size_FOU)
        plt.plot(t32_MUSCL, e32_L1_MUSCL, '-go', mfc='g', markersize=size_MUSCL)
        plt.plot(t64_FOU, e64_L1_FOU, '-rs', mfc='none', markersize=size_FOU)
        plt.plot(t64_MUSCL, e64_L1_MUSCL, '-gs', mfc='g', markersize=size_MUSCL)
        plt.plot(t128_FOU, e128_L1_FOU, '-rv', mfc='none', markersize=size_FOU)
        plt.plot(t128_MUSCL, e128_L1_MUSCL, '-gv', mfc='g', markersize=size_MUSCL)
        plt.plot(t256_FOU, e256_L1_FOU, '-rp', mfc='none', markersize=size_FOU)
        plt.plot(t256_MUSCL, e256_L1_MUSCL, '-gp', mfc='g', markersize=size_MUSCL)
        plt.plot(t512_FOU, e512_L1_FOU, '-rD', mfc='none', markersize=size_FOU)
        plt.plot(t512_MUSCL, e512_L1_MUSCL, '-gD', mfc='g', markersize=size_MUSCL)
        plt.plot(t1024_FOU, e1024_L1_FOU, '-r*', mfc='none', markersize=size_FOU)
        plt.plot(t2048_FOU, e2048_L1_FOU, '-r<', mfc='none', markersize=size_FOU)
        legend_elements = [Line2D([0], [0], color='g', label='MUSCL'),
                   Line2D([0], [0], color='w', markerfacecolor='g', marker='o', markeredgecolor='g', label='32 CV', markersize=size_MUSCL),#, mfc='none'),
                   Line2D([0], [0], marker='s', color='w', markerfacecolor='g', label='64 CV', markeredgecolor='g', markersize=size_MUSCL),#, mfc='none'),
                   Line2D([0], [0], marker='v', color='w', markerfacecolor='g',label='128 CV',markeredgecolor='g', markersize=size_MUSCL),#, mfc='none'),
                   Line2D([0], [0], marker='p', color='w', markerfacecolor='g', label='256 CV',markeredgecolor='g', markersize=size_MUSCL),#, mfc='none'),
                   Line2D([0], [0], marker='D', color='w', markerfacecolor='g', label='512 CV',markeredgecolor='g', markersize=size_MUSCL),#, mfc='none'),
                   Line2D([0], [0], marker='D', color='w', markerfacecolor='none', label=' ',markeredgecolor='w', markersize=size_MUSCL),#, mfc='none'),
                   Line2D([0], [0], marker='D', color='w', markerfacecolor='none', label=' ',markeredgecolor='w', markersize=size_MUSCL),#, mfc='none'),
                   Line2D([0], [0], color='r', label='FOU'),
                   Line2D([0], [0], color='w', marker='o', markeredgecolor='r', label='32 CV', markersize=size_FOU, mfc='none'),
                   Line2D([0], [0], marker='s', color='w', label='64 CV', markeredgecolor='r', markersize=size_FOU, mfc='none'),
                   Line2D([0], [0], marker='v', color='w', label='128 CV',markeredgecolor='r', markersize=size_FOU, mfc='none'),
                   Line2D([0], [0], marker='p', color='w', label='256 CV',markeredgecolor='r', markersize=size_FOU, mfc='none'),
                   Line2D([0], [0], marker='D', color='w', label='512 CV',markeredgecolor='r', markersize=size_FOU, mfc='none'),
                   Line2D([0], [0], marker='*', color='w', label='1024 CV',markeredgecolor='r', markersize=size_FOU, mfc='none'),
                   Line2D([0], [0], marker='<', color='w', label='2048 CV',markeredgecolor='r', markersize=size_FOU, mfc='none'),]
        # Create the figure
        plt.legend(handles=legend_elements, ncol=2,)
        #plt.legend(('FOU', 'MUSCL', 'First Order'))
        plt.grid()
        plt.xlabel('Tempo computacional [s]')
        plt.ylabel('$E_{L1}$')
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig('results/compositional/TCC2/5k_time_eL1.png')


        plt.figure(10)
        plt.plot(x_4000, zCO2_FOU_4000, 'k')
        #plt.plot(x_8, Sg_MUSCL_8, '-r', mfc='none')
        #plt.plot(x_16, Sg_MUSCL_16, '-y', mfc='none')
        plt.plot(x_32, zCO2_MUSCL_32, '-go', mfc='none', markersize=size)
        plt.plot(x_64, zCO2_MUSCL_64, '-bp', mfc='none', markersize=size)
        plt.plot(x_128, zCO2_MUSCL_128, '-mD', mfc='none', markersize=size)
        plt.plot(x_256, zCO2_MUSCL_256, '-cs', mfc='none', markersize=size)
        #plt.plot(x_512, zCO2_MUSCL_512, '-rv', mfc='none', markersize=size)
        plt.legend(('FOU-4000', 'MUSCL-32', 'MUSCL-64', \
            'MUSCL-128', 'MUSCL-256'), loc=3)
        plt.grid()
        plt.xlim((0, 1.1))
        plt.ylim((0.75, 1))
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel('Distância')
        plt.savefig('results/compositional/TCC2/5k_zCO2_MUSCL.png')

        plt.figure(11)
        plt.plot(x_4000, zCO2_FOU_4000, 'k')
        #plt.plot(x_8, Sg_FOU_8, '-r', mfc='none')
        #plt.plot(x_16, Sg_FOU_16, '-y', mfc='none')
        plt.plot(x_32, zCO2_FOU_32, '-go', mfc='none', markersize=size)
        plt.plot(x_64, zCO2_FOU_64, '-bp', mfc='none', markersize=size)
        plt.plot(x_128, zCO2_FOU_128, '-mD', mfc='none', markersize=size)
        plt.plot(x_256, zCO2_FOU_256, '-cs', mfc='none', markersize=size)
        #plt.plot(x_512, zCO2_FOU_512, '-rv', mfc='none', markersize=size)
        plt.legend(('FOU-4000', 'FOU-32', 'FOU-64', \
            'FOU-128', 'FOU-256'), loc=3)
        plt.grid()
        plt.xlim((0, 1.1))
        plt.ylim((0.75, 1))
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel('Distância')
        plt.savefig('results/compositional/TCC2/5k_zCO2_FOU.png')

        import pdb; pdb.set_trace()
