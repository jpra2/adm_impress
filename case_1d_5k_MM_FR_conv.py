import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
L=1.5

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

        '''datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8000_FOU_8726.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_8000 = data[6]
            Sg_FOU_8000 = data[7]
            zCO2_FOU_8000 = data[10][0]
            x_8000 = np.linspace(0,1.5,8000)
            f = interp1d(x_8000,zCO2_FOU_8000)'''

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8000_FOU_25227.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_8000 = data[6]
            Sg_FOU_8000 = data[7]
            zCO2_FOU_8000 = data[10][0]
            zC1_FOU_8000 = data[10][1]
            x_8000 = np.linspace(0,1.5,8000)
            f = interp1d(x_8000,zCO2_FOU_8000)
            fC1 = interp1d(x_8000,zC1_FOU_8000)

        '------------------------------- FR -----------------------------------'
        'FR2'
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FR2_CFL09m_RK3_MLPmod_36.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_8 = data[7]
            z_FR2_8 = data[10]
            t_FR2_8 = data[2]
            n=8
            x_8 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e8_L1_FR2 = (sum(abs(f(x_8)-data[10][0]))*(L/n))
            e8_L1_FR2_C1 = (sum(abs(fC1(x_8)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FR2_CFL09m_RK3_MLPmod_72.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_16 = data[7]
            z_FR2_16 = data[10]
            t_FR2_16 = data[2]
            n=16
            x_16 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e16_L1_FR2 = (sum(abs(f(x_16)-data[10][0]))*(L/n))
            e16_L1_FR2_C1 = (sum(abs(fC1(x_16)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FR2_CFL09m_RK3_MLPmod_145.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_32 = data[7]
            z_FR2_32 = data[10]
            t_FR2_32 = data[2]
            n=32
            x_32 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e32_L1_FR2 = (sum(abs(f(x_32)-data[10][0]))*(L/n))
            e32_L1_FR2_C1 = (sum(abs(fC1(x_32)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FR2_CFL09_RK3_MLPmod_295.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_64 = data[7]
            z_FR2_64 = data[10]
            t_FR2_64 = data[2]
            n=64
            x_64 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e64_L1_FR2 = (sum(abs(f(x_64)-data[10][0]))*(L/n))
            e64_L1_FR2_C1 = (sum(abs(fC1(x_64)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR2_CFL09m_RK3_MLPmod_589.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_128 = data[7]
            z_FR2_128 = data[10]
            t_FR2_128 = data[2]
            n=128
            x_128 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e128_L1_FR2 = (sum(abs(f(x_128)-data[10][0]))*(L/n))
            e128_L1_FR2_C1 = (sum(abs(fC1(x_128)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_160_FR2_CFL09m_RK3_MLPmod_734.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_160 = data[7]
            z_FR2_160 = data[10]
            t_FR2_160 = data[2]
            n=160
            x_160 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e160_L1_FR2 = (sum(abs(f(x_160)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FR2_CFL09m_RK3_MLPmod_1176.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_256 = data[7]
            z_FR2_256 = data[10]
            t_FR2_256 = data[2]
            n=256
            x_256 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e256_L1_FR2 = (sum(abs(f(x_256)-data[10][0]))*(L/n))
            e256_L1_FR2_C1 = (sum(abs(fC1(x_256)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_FR2_CFL09m_RK3_MLPmod_2920.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_512 = data[7]
            z_FR2_512 = data[10]
            t_FR2_512 = data[2]
            n=512
            x_512 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e512_L1_FR2 = (sum(abs(f(x_512)-data[10][0]))*(L/n))
            e512_L1_FR2_C1 = (sum(abs(fC1(x_512)-data[10][1]))*(L/n))


        'FR3'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FR3_CFL09m_RK3_MLPmod_49.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_8 = data[7]
            z_FR3_8 = data[10]
            t_FR3_8 = data[2]
            n=8
            x_8 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e8_L1_FR3 = (sum(abs(f(x_8)-data[10][0]))*(L/n))
            e8_L1_FR3_C1 = (sum(abs(fC1(x_8)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FR3_CFL09m_RK3_MLPmod_100.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_16 = data[7]
            z_FR3_16 = data[10]
            t_FR3_16 = data[2]
            n=16
            x_16 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e16_L1_FR3 = (sum(abs(f(x_16)-data[10][0]))*(L/n))
            e16_L1_FR3_C1 = (sum(abs(fC1(x_16)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FR3_CFL09m_RK3_MLPmod_207.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_32 = data[7]
            z_FR3_32 = data[10]
            t_FR3_32 = data[2]
            n=32
            x_32 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e32_L1_FR3 = (sum(abs(f(x_32)-data[10][0]))*(L/n))
            e32_L1_FR3_C1 = (sum(abs(fC1(x_32)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FR3_CFL09m_RK3_MLPmod_512.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_64 = data[7]
            z_FR3_64 = data[10]
            t_FR3_64 = data[2]
            n=64
            x_64 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e64_L1_FR3 = (sum(abs(f(x_64)-data[10][0]))*(L/n))
            e64_L1_FR3_C1 = (sum(abs(fC1(x_64)-data[10][1]))*(L/n))

        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR3_CFL09m_RK3_MLPmod_897.npy', allow_pickle=True)
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR3_CFL09m_RK3_MLPmod_914.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_128 = data[7]
            z_FR3_128 = data[10]
            t_FR3_128 = data[2]
            n=128
            x_128 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e128_L1_FR3 = (sum(abs(f(x_128)-data[10][0]))*(L/n))
            e128_L1_FR3_C1 = (sum(abs(fC1(x_128)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_160_FR3_CFL09m_RK3_MLPmod_1151.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_160 = data[7]
            z_FR3_160 = data[10]
            t_FR3_160 = data[2]
            n=160
            x_160 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e160_L1_FR3 = (sum(abs(f(x_160)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FR3_CFL09m_RK3_MLPmod_1895.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_256 = data[7]
            z_FR3_256 = data[10]
            t_FR3_256 = data[2]
            n=256
            x_256 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e256_L1_FR3 = (sum(abs(f(x_256)-data[10][0]))*(L/n))
            e256_L1_FR3_C1 = (sum(abs(fC1(x_256)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_FR3_CFL09m_RK3_MLPmod_4346.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_512 = data[7]
            z_FR3_512 = data[10]
            t_FR3_512 = data[2]
            n=512
            x_512 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e512_L1_FR3 = (sum(abs(f(x_512)-data[10][0]))*(L/n))
            e512_L1_FR3_C1 = (sum(abs(fC1(x_512)-data[10][1]))*(L/n))

        'FR4'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FR4_CFL09m_RK3_MLPmod_60.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_8 = data[7]
            z_FR4_8 = data[10]
            t_FR4_8 = data[2]
            n=8
            x_8 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e8_L1_FR4 = (sum(abs(f(x_8)-data[10][0]))*(L/n))
            e8_L1_FR4_C1 = (sum(abs(fC1(x_8)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FR4_CFL09m_RK3_MLPmod_122.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_16 = data[7]
            z_FR4_16 = data[10]
            t_FR4_16 = data[2]
            n=16
            x_16 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e16_L1_FR4 = (sum(abs(f(x_16)-data[10][0]))*(L/n))
            e16_L1_FR4_C1 = (sum(abs(fC1(x_16)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FR4_CFL09m_RK3_MLPmod_278.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_32 = data[7]
            z_FR4_32 = data[10]
            t_FR4_32 = data[2]
            n=32
            x_32 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e32_L1_FR4 = (sum(abs(f(x_32)-data[10][0]))*(L/n))
            e32_L1_FR4_C1 = (sum(abs(fC1(x_32)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FR4_CFL09m_RK3_MLPmod_584.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_64 = data[7]
            z_FR4_64 = data[10]
            t_FR4_64 = data[2]
            n=64
            x_64 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e64_L1_FR4 = (sum(abs(f(x_64)-data[10][0]))*(L/n))
            e64_L1_FR4_C1 = (sum(abs(fC1(x_64)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR4_CFL09m_RK3_MLPmod_1249.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_128 = data[7]
            z_FR4_128 = data[10]
            t_FR4_128 = data[2]
            n=128
            x_128 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e128_L1_FR4 = (sum(abs(f(x_128)-data[10][0]))*(L/n))
            e128_L1_FR4_C1 = (sum(abs(fC1(x_128)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_160_FR4_CFL09m_RK3_MLPmod_1579.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_160 = data[7]
            z_FR4_160 = data[10]
            t_FR4_160 = data[2]
            n=160
            x_160 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e160_L1_FR4 = (sum(abs(f(x_160)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FR4_CFL09m_RK3_MLPmod_2606.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_256 = data[7]
            z_FR4_256 = data[10]
            t_FR4_256 = data[2]
            n=256
            x_256 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e256_L1_FR4 = (sum(abs(f(x_256)-data[10][0]))*(L/n))
            e256_L1_FR4_C1 = (sum(abs(fC1(x_256)-data[10][1]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_FR4_CFL09m_RK3_MLPmod_5207.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_512 = data[7]
            z_FR4_512 = data[10]
            t_FR4_512 = data[2]
            n=512
            x_512 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e512_L1_FR4 = (sum(abs(f(x_512)-data[10][0]))*(L/n))
            e512_L1_FR4_C1 = (sum(abs(fC1(x_512)-data[10][1]))*(L/n))

        'MUSCL + LLF'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_MUSCL_CFL09_RK3_VA_31.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_8 = data[7]
            z_MUSCL_8 = data[10]
            n = 8
            e8_L1_MUSCL = (sum(abs(f(x_8)-data[10][0]))*(L/n))
            e8_L1_MUSCL_C1 = (sum(abs(fC1(x_8)-data[10][1]))*(L/n))
            t_MUSCL_8 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_MUSCL_CFL09_RK3_VA_66.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_16 = data[7]
            z_MUSCL_16 = data[10]
            n = 16
            e16_L1_MUSCL = (sum(abs(f(x_16)-data[10][0]))*(L/n))
            e16_L1_MUSCL_C1 = (sum(abs(fC1(x_16)-data[10][1]))*(L/n))
            t_MUSCL_16 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_MUSCL_CFL09_RK3_VA_137.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_32 = data[7]
            z_MUSCL_32 = data[10]
            n = 32
            e32_L1_MUSCL = (sum(abs(f(x_32)-data[10][0]))*(L/n))
            e32_L1_MUSCL_C1 = (sum(abs(fC1(x_32)-data[10][1]))*(L/n))
            t_MUSCL_32 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_MUSCL_CFL09_RK3_VA_284.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_64 = data[7]
            z_MUSCL_64 = data[10]
            n = 64
            e64_L1_MUSCL = (sum(abs(f(x_64)-data[10][0]))*(L/n))
            e64_L1_MUSCL_C1 = (sum(abs(fC1(x_64)-data[10][1]))*(L/n))
            t_MUSCL_64 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_MUSCL_CFL09_RK3_VA_575.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_128 = data[7]
            z_MUSCL_128 = data[10]
            n = 128
            e128_L1_MUSCL = (sum(abs(f(x_128)-data[10][0]))*(L/n))
            e128_L1_MUSCL_C1 = (sum(abs(fC1(x_128)-data[10][1]))*(L/n))
            t_MUSCL_128 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_MUSCL_CFL09_RK3_VA_1154.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_256 = data[7]
            z_MUSCL_256 = data[10]
            n = 256
            e256_L1_MUSCL = (sum(abs(f(x_256)-data[10][0]))*(L/n))
            e256_L1_MUSCL_C1 = (sum(abs(fC1(x_256)-data[10][1]))*(L/n))
            t_MUSCL_256 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_MUSCL_CFL09_RK3_VA_2897.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_512 = data[7]
            z_MUSCL_512 = data[10]
            n = 512
            e512_L1_MUSCL = (sum(abs(f(x_512)-data[10][0]))*(L/n))
            e512_L1_MUSCL_C1 = (sum(abs(fC1(x_512)-data[10][1]))*(L/n))
            t_MUSCL_512 = data[2]

        'FOU'
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FOU_CFL09_11.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_8 = data[6]
            Sg_FOU_8 = data[7]
            z_FOU_8 = data[10]
            n = 8
            e8_L1_FOU = (sum(abs(f(x_8)- data[10][0]))*(L/n))
            e8_L1_FOU_C1 = (sum(abs(fC1(x_8)-data[10][1]))*(L/n))
            t_FOU_8 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FOU_CFL09_21.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_16 = data[6]
            Sg_FOU_16 = data[7]
            z_FOU_16 = data[10]
            n = 16
            e16_L1_FOU = (sum(abs(f(x_16)-data[10][0]))*(L/n))
            e16_L1_FOU_C1 = (sum(abs(fC1(x_16)-data[10][1]))*(L/n))
            t_FOU_16 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FOU_CFL09_42.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_32 = data[6]
            Sg_FOU_32 = data[7]
            z_FOU_32 = data[10]
            n = 32
            e32_L1_FOU = (sum(abs(f(x_32)-data[10][0]))*(L/n))
            e32_L1_FOU_C1 = (sum(abs(fC1(x_32)-data[10][1]))*(L/n))
            t_FOU_32 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FOU_CFL09_85.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_64 = data[6]
            Sg_FOU_64 = data[7]
            z_FOU_64 = data[10]
            n = 64
            e64_L1_FOU = (sum(abs(f(x_64)-data[10][0]))*(L/n))
            e64_L1_FOU_C1 = (sum(abs(fC1(x_64)-data[10][1]))*(L/n))
            t_FOU_64 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FOU_CFL09_173.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_128 = data[6]
            Sg_FOU_128 = data[7]
            z_FOU_128 = data[10]
            n = 128
            e128_L1_FOU = (sum(abs(f(x_128)-data[10][0]))*(L/n))
            e128_L1_FOU_C1 = (sum(abs(fC1(x_128)-data[10][1]))*(L/n))
            t_FOU_128 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FOU_CFL09_357.npy', allow_pickle=True)
        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FOU_t_404.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_256 = data[6]
            Sg_FOU_256 = data[7]
            z_FOU_256 = data[10]
            n = 256
            e256_L1_FOU = (sum(abs(f(x_256)-data[10][0]))*(L/n))
            e256_L1_FOU_C1 = (sum(abs(fC1(x_256)-data[10][1]))*(L/n))
            t_FOU_256 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_FOU_CFL09_726.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_512 = data[6]
            Sg_FOU_512 = data[7]
            z_FOU_512 = data[10]
            x_512 = np.linspace(0,1.5,len(So_FOU_512))
            n = 512
            e512_L1_FOU = (sum(abs(f(x_512)-data[10][0]))*(L/n))
            e512_L1_FOU_C1 = (sum(abs(fC1(x_512)-data[10][1]))*(L/n))
            t_FOU_512 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_1024_FOU_CFL09_1556.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_1024 = data[6]
            Sg_FOU_1024 = data[7]
            z_FOU_1024 = data[10]

            x_1024 = np.linspace(0,1.5,len(So_FOU_1024))
            n = 1024
            e1024_L1_FOU = (sum(abs(f(x_1024)-data[10][0]))*(L/n))
            t_FOU_1024 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_2048_FOU_5139.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_2048 = data[6]
            Sg_FOU_2048 = data[7]
            z_FOU_2048 = data[10]

            x_2048 = np.linspace(0,1.5,len(So_FOU_2048))
            n = 2048
            e2048_L1_FOU = (sum(abs(f(x_2048)-data[10][0]))*(L/n))
            t_FOU_2048 = data[2]


        plt.figure(1)
        #plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        plt.plot(x_8000, zC1_FOU_8000, 'k')
        plt.plot(x_128, z_FOU_128[1,:], '-bo')
        plt.plot(x_128, z_MUSCL_128[1,:], '-mv')
        plt.plot(x_128, z_FR2_128[1,:], '-gs')
        plt.plot(x_128, z_FR3_128[1,:], '-y<')
        plt.plot(x_128, z_FR4_128[1,:], '-r*')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.xlim((1, 1.14))
        plt.ylim((0, 0.65))
        plt.title('NVCM example with 128x1x1 mesh')
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zC1_FR_hMLPmod_128.png')

        plt.figure(2)
        plt.plot(x_8000, zCO2_FOU_8000, 'k')
        plt.plot(x_128, z_FOU_128[0,:], '-b')
        plt.plot(x_128, z_MUSCL_128[0,:], '-m')
        plt.plot(x_128, z_FR2_128[0,:], '-g')
        plt.plot(x_128, z_FR3_128[0,:], '-y')
        plt.plot(x_128, z_FR4_128[0,:], '-r')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.xlim((0.55, 1.1))
        plt.ylim((0.75, 1))
        plt.title('NVCM example with 128x1x1 mesh')
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zCO2_FR_hMLPmod_128.png')

        plt.figure(3)
        plt.plot(x_8000, Sg_FOU_8000, 'k')
        plt.plot(x_128, Sg_FOU_128, '-b')
        plt.plot(x_128, Sg_MUSCL_128, '-m')
        #plt.plot(x_128, Sg_FR2_128, '-y')
        plt.plot(x_128, Sg_FR2_128, '-g')
        plt.plot(x_128, Sg_FR3_128, '-y')
        plt.plot(x_128, Sg_FR4_128, '-r')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.xlim((0.4, 1.2))
        #plt.ylim((0, 0.65))
        plt.title('NVCM example with 128x1x1 mesh')
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_Sg_FR_hMLPmod_128.png')


        plt.figure(4)
        e_FR2 = np.log10(np.array([e8_L1_FR2, e16_L1_FR2, e32_L1_FR2, e64_L1_FR2, e128_L1_FR2, e256_L1_FR2, e512_L1_FR2]))
        e_FR3 = np.log10(np.array([e8_L1_FR3, e16_L1_FR3, e32_L1_FR3, e64_L1_FR3, e128_L1_FR3, e256_L1_FR3, e512_L1_FR3]))
        e_FR4 = np.log10(np.array([e8_L1_FR4, e16_L1_FR4, e32_L1_FR4, e64_L1_FR4, e128_L1_FR4, e256_L1_FR4, e512_L1_FR4]))
        e_MUSCL = np.log10(np.array([e8_L1_MUSCL, e16_L1_MUSCL, e32_L1_MUSCL, e64_L1_MUSCL, e128_L1_MUSCL, e256_L1_MUSCL, e512_L1_MUSCL]))
        e_FOU = np.log10(np.array([e8_L1_FOU, e16_L1_FOU, e32_L1_FOU, e64_L1_FOU, e128_L1_FOU, e256_L1_FOU, e512_L1_FOU]))

        x = np.log10(np.array([8,16,32,64,128,256, 512]))
        y = -x
        plt.plot(x,y,'-k')

        plt.plot(x, e_FOU, '-bo', mfc='none')
        plt.plot(x, e_MUSCL, '-mv', mfc='none')
        plt.plot(x, e_FR2, '-gs', mfc='none')
        plt.plot(x, e_FR3, '-y<', mfc='none')
        plt.plot(x, e_FR4, '-r*', mfc='none')
        plt.title('Error using zCO2')
        #plt.plot(x_256, z_FR2_256[1,:], '-gv', mfc='none')
        plt.legend(('1st order','FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.ylabel('log(|$E_{L1}$|)')
        plt.xlabel('log($n_b$)')
        plt.savefig('results/compositional/FR/5k_erro_zCO2_FR_hMLPmod.png')

        plt.figure(5)
        e_FR2 = np.log10(np.array([e8_L1_FR2_C1, e16_L1_FR2_C1, e32_L1_FR2_C1, e64_L1_FR2_C1, e128_L1_FR2_C1, e256_L1_FR2_C1, e512_L1_FR2_C1]))
        e_FR3 = np.log10(np.array([e8_L1_FR3_C1, e16_L1_FR3_C1, e32_L1_FR3_C1, e64_L1_FR3_C1, e128_L1_FR3_C1, e256_L1_FR3_C1, e512_L1_FR3_C1]))
        e_FR4 = np.log10(np.array([e8_L1_FR4_C1, e16_L1_FR4_C1, e32_L1_FR4_C1, e64_L1_FR4_C1, e128_L1_FR4_C1, e256_L1_FR4_C1, e512_L1_FR4_C1]))
        e_MUSCL = np.log10(np.array([e8_L1_MUSCL_C1, e16_L1_MUSCL_C1, e32_L1_MUSCL_C1, e64_L1_MUSCL_C1, e128_L1_MUSCL_C1, e256_L1_MUSCL_C1, e512_L1_MUSCL_C1]))
        e_FOU = np.log10(np.array([e8_L1_FOU_C1, e16_L1_FOU_C1, e32_L1_FOU_C1, e64_L1_FOU_C1, e128_L1_FOU_C1, e256_L1_FOU_C1, e512_L1_FOU_C1]))

        x = np.log10(np.array([8,16,32,64,128,256, 512]))
        y = -x
        plt.plot(x[3:],y[3:],'-k')

        plt.plot(x[3:], e_FOU[3:], '-bo', mfc='none')
        plt.plot(x[3:], e_MUSCL[3:], '-mv', mfc='none')
        plt.plot(x[3:], e_FR2[3:], '-gs', mfc='none')
        plt.plot(x[3:], e_FR3[3:], '-y<', mfc='none')
        plt.plot(x[3:], e_FR4[3:], '-r*', mfc='none')
        #plt.plot(x_256, z_FR2_256[1,:], '-gv', mfc='none')
        plt.legend(('1st order','FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.title('Error using zC1')
        plt.ylabel('log(|$E_{L1}$|)')
        plt.xlabel('log($n_b$)')
        plt.savefig('results/compositional/FR/5k_erro_zC1_FR_hMLPmod.png')

        plt.figure(6)
        plt.plot(x_8000, zC1_FOU_8000, 'k')
        plt.plot(x_64, z_FR2_64[1,:], '-ys')
        plt.plot(x_64, z_FR3_64[1,:], '-bv')
        plt.plot(x_64, z_FR4_64[1,:], '-ro')
        plt.plot(x_64, z_MUSCL_64[1,:], '-g*')
        plt.legend(('Reference', 'FR-P1', 'FR-P2', 'FR-P3','MUSCL'))
        plt.grid()
        plt.xlim((1, 1.14))
        plt.ylim((0, 0.65))
        plt.title('NVCM example with 64x1x1 mesh')
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zC1_FR_hMLPmod_64.png')

        plt.figure(7)
        plt.plot(x_8000, zCO2_FOU_8000, 'k')
        plt.plot(x_64, z_FOU_64[0,:], '-b')
        plt.plot(x_64, z_MUSCL_64[0,:], '-m')
        plt.plot(x_64, z_FR2_64[0,:], '-g')
        plt.plot(x_64, z_FR3_64[0,:], '-y')
        plt.plot(x_64, z_FR4_64[0,:], '-r')
        plt.legend(('Reference','FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.xlim((0.55, 1.1))
        plt.ylim((0.75, 1))
        plt.title('NVCM example with 64x1x1 mesh')
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zCO2_FR_hMLPmod_64.png')

        plt.figure(8)
        plt.plot(x_8000, zCO2_FOU_8000, 'k')
        plt.plot(x_256, z_FOU_256[0,:], '-b')
        plt.plot(x_256, z_MUSCL_256[0,:], '-m')
        plt.plot(x_256, z_FR2_256[0,:], '-g')
        plt.plot(x_256, z_FR3_256[0,:], '-y')
        plt.plot(x_256, z_FR4_256[0,:], '-r')
        plt.legend(('Reference','FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.xlim((0.55, 1.1))
        plt.ylim((0.75, 1))
        plt.title('NVCM example with 256x1x1 mesh')
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zCO2_FR_hMLPmod_256.png')

        plt.figure(9)
        plt.plot(x_8000, Sg_FOU_8000, 'k')
        plt.plot(x_256, Sg_FOU_256, '-b')
        plt.plot(x_256, Sg_MUSCL_256, '-m')
        plt.plot(x_256, Sg_FR2_256, '-g')
        plt.plot(x_256, Sg_FR3_256, '-y')
        plt.plot(x_256, Sg_FR4_256, '-r')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.xlim((0.4, 1.2))
        #plt.ylim((0, 0.65))
        plt.title('NVCM example with 256x1x1 mesh')
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_Sg_FR_hMLPmod_256.png')

        plt.figure(10)
        #plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        plt.plot(x_8000, zC1_FOU_8000, 'k')
        plt.plot(x_256, z_FOU_256[1,:], '-bo')
        plt.plot(x_256, z_MUSCL_256[1,:], '-mv')
        plt.plot(x_256, z_FR2_256[1,:], '-gs')
        plt.plot(x_256, z_FR3_256[1,:], '-y<')
        plt.plot(x_256, z_FR4_256[1,:], '-r*')
        plt.legend(('Reference','FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.xlim((1, 1.14))
        plt.ylim((0, 0.65))
        plt.title('NVCM example with 256x1x1 mesh')
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zC1_FR_hMLPmod_256.png')

        plt.figure(11)
        plt.plot(x_8000, zCO2_FOU_8000, 'k')
        plt.plot(x_512, z_FOU_512[0,:], '-b')
        plt.plot(x_512, z_MUSCL_512[0,:], '-m')
        plt.plot(x_512, z_FR2_512[0,:], '-g')
        plt.plot(x_512, z_FR3_512[0,:], '-y')
        plt.plot(x_512, z_FR4_512[0,:], '-r')
        plt.legend(('Reference','FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.xlim((0.55, 1.1))
        plt.ylim((0.75, 1))
        plt.title('NVCM example with 512x1x1 mesh')
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zCO2_FR_hMLPmod_512.png')

        plt.figure(12)
        plt.plot(x_8000, Sg_FOU_8000, 'k')
        plt.plot(x_512, Sg_FOU_512, '-b')
        plt.plot(x_512, Sg_MUSCL_512, '-m')
        plt.plot(x_512, Sg_FR2_512, '-g')
        plt.plot(x_512, Sg_FR3_512, '-y')
        plt.plot(x_512, Sg_FR4_512, '-r')
        plt.legend(('Reference','FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.xlim((0.4, 1.2))
        #plt.ylim((0, 0.65))
        plt.title('NVCM example with 512x1x1 mesh')
        plt.ylabel('$z_{CO_2}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_Sg_FR_hMLPmod_512.png')

        plt.figure(13)
        plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        #plt.plot(x_8000, zC1_FOU_8000, 'k')
        plt.plot(x_512, z_FOU_512[1,:], '-bo')
        plt.plot(x_512, z_MUSCL_512[1,:], '-mv')
        plt.plot(x_512, z_FR2_512[1,:], '-gs')
        plt.plot(x_512, z_FR3_512[1,:], '-y<')
        plt.plot(x_512, z_FR4_512[1,:], '-r*')
        plt.legend(('Reference','FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.xlim((1, 1.14))
        plt.ylim((0, 0.65))
        plt.title('NVCM example with 512x1x1 mesh')
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zC1_FR_hMLPmod_512.png')
        import pdb; pdb.set_trace()
