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

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_4000_FOU_8726.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_4000 = data[6]
            Sg_FOU_4000 = data[7]
            zCO2_FOU_4000 = data[10][0]
            x_4000 = np.linspace(0,1.5,4000)
            f = interp1d(x_4000,zCO2_FOU_4000)

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

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FR2_CFL09m_RK3_MLPmod_72.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_16 = data[7]
            z_FR2_16 = data[10]
            t_FR2_16 = data[2]
            n=16
            x_16 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e16_L1_FR2 = (sum(abs(f(x_16)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FR2_CFL09m_RK3_MLPmod_145.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_32 = data[7]
            z_FR2_32 = data[10]
            t_FR2_32 = data[2]
            n=32
            x_32 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e32_L1_FR2 = (sum(abs(f(x_32)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FR2_CFL09_RK3_MLPmod_295.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_64 = data[7]
            z_FR2_64 = data[10]
            t_FR2_64 = data[2]
            n=64
            x_64 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e64_L1_FR2 = (sum(abs(f(x_64)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR2_CFL09_RK3_MLPmod_584.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_128 = data[7]
            z_FR2_128 = data[10]
            t_FR2_128 = data[2]
            n=128
            x_128 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e128_L1_FR2 = (sum(abs(f(x_128)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FR2_CFL09m_RK3_MLPmod_1172.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_256 = data[7]
            z_FR2_256 = data[10]
            t_FR2_256 = data[2]
            n=256
            x_256 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e256_L1_FR2 = (sum(abs(f(x_256)-data[10][0]))*(L/n))



        'FR3'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FR3_CFL09m_RK3_MLPmod_49.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_8 = data[7]
            z_FR3_8 = data[10]
            t_FR3_8 = data[2]
            n=8
            x_8 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e8_L1_FR3 = (sum(abs(f(x_8)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FR3_CFL09m_RK3_MLPmod_100.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_16 = data[7]
            z_FR3_16 = data[10]
            t_FR3_16 = data[2]
            n=16
            x_16 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e16_L1_FR3 = (sum(abs(f(x_16)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FR3_CFL09m_RK3_MLPmod_207.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_32 = data[7]
            z_FR3_32 = data[10]
            t_FR3_32 = data[2]
            n=32
            x_32 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e32_L1_FR3 = (sum(abs(f(x_32)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FR3_CFL09m_RK3_MLPmod_512.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_64 = data[7]
            z_FR3_64 = data[10]
            t_FR3_64 = data[2]
            n=64
            x_64 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e64_L1_FR3 = (sum(abs(f(x_64)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR3_CFL09m_RK3_MLPmod_1223.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_128 = data[7]
            z_FR3_128 = data[10]
            t_FR3_128 = data[2]
            n=128
            x_128 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e128_L1_FR3 = (sum(abs(f(x_128)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FR3_CFL09m_RK3_MLPmod_2612.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_256 = data[7]
            z_FR3_256 = data[10]
            t_FR3_256 = data[2]
            n=256
            x_256 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e256_L1_FR3 = (sum(abs(f(x_256)-data[10][0]))*(L/n))

        'FR3+TVB'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FR3_CFL09m_RK3_TVD_56.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_TVB_8 = data[7]
            z_FR3_TVB_8 = data[10]
            t_FR3_TVB_8 = data[2]
            n=8
            e8_L1_FR3_TVB = (sum(abs(f(x_8)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FR3_CFL09m_RK3_TVD_114.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_TVB_16 = data[7]
            z_FR3_TVB_16 = data[10]
            t_FR3_TVB_16 = data[2]
            n=16
            e16_L1_FR3_TVB = (sum(abs(f(x_16)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FR3_CFL09m_RK3_TVD_234.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_TVB_32 = data[7]
            z_FR3_TVB_32 = data[10]
            t_FR3_TVB_32 = data[2]
            n=32
            x_32 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e32_L1_FR3_TVB = (sum(abs(f(x_32)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FR3_CFL09m_RK3_TVD_481.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_TVB_64 = data[7]
            z_FR3_TVB_64 = data[10]
            t_FR3_TVB_64 = data[2]
            n=64
            x_64 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e64_L1_FR3_TVB = (sum(abs(f(x_64)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR3_CFL09m_RK3_TVD_973.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_TVB_128 = data[7]
            z_FR3_TVB_128 = data[10]
            t_FR3_TVB_128 = data[2]
            n=128
            x_128 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e128_L1_FR3_TVB = (sum(abs(f(x_128)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FR3_CFL09m_RK3_TVD_1968.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_TVB_256 = data[7]
            z_FR3_TVB_256 = data[10]
            t_FR3_TVB_256 = data[2]
            n=256
            x_256 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e256_L1_FR3_TVB = (sum(abs(f(x_256)-data[10][0]))*(L/n))


        'FR4+TVB'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FR4_CFL09m_RK3_TVD_78.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_TVB_8 = data[7]
            z_FR4_TVB_8 = data[10]
            t_FR4_TVB_8 = data[2]
            n=8
            e8_L1_FR4_TVB = (sum(abs(f(x_8)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FR4_CFL09m_RK3_TVD_159.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_TVB_16 = data[7]
            z_FR4_TVB_16 = data[10]
            t_FR4_TVB_16 = data[2]
            n=16
            e16_L1_FR4_TVB = (sum(abs(f(x_16)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FR4_CFL09m_RK3_TVD_328.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_TVB_32 = data[7]
            z_FR4_TVB_32 = data[10]
            t_FR4_TVB_32 = data[2]
            n=32
            x_32 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e32_L1_FR4_TVB = (sum(abs(f(x_32)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FR4_CFL09m_RK3_TVD_673.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_TVB_64 = data[7]
            z_FR4_TVB_64 = data[10]
            t_FR4_TVB_64 = data[2]
            n=64
            x_64 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e64_L1_FR4_TVB = (sum(abs(f(x_64)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR4_CFL09m_RK3_TVD_1362.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_TVB_128 = data[7]
            z_FR4_TVB_128 = data[10]
            t_FR4_TVB_128 = data[2]
            n=128
            x_128 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e128_L1_FR4_TVB = (sum(abs(f(x_128)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FR4_CFL09m_RK3_TVD_2771.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_TVB_256 = data[7]
            z_FR4_TVB_256 = data[10]
            t_FR4_TVB_256 = data[2]
            n=256
            x_256 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e256_L1_FR4_TVB = (sum(abs(f(x_256)-data[10][0]))*(L/n))


        'FR4'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FR4_CFL09m_RK3_MLPmod_60.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_8 = data[7]
            z_FR4_8 = data[10]
            t_FR4_8 = data[2]
            n=8
            x_8 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e8_L1_FR4 = (sum(abs(f(x_8)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FR4_CFL09m_RK3_MLPmod_122.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_16 = data[7]
            z_FR4_16 = data[10]
            t_FR4_16 = data[2]
            n=16
            x_16 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e16_L1_FR4 = (sum(abs(f(x_16)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FR4_CFL09m_RK3_MLPmod_278.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_32 = data[7]
            z_FR4_32 = data[10]
            t_FR4_32 = data[2]
            n=32
            x_32 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e32_L1_FR4 = (sum(abs(f(x_32)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FR4_CFL09m_RK3_MLPmod_640.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_64 = data[7]
            z_FR4_64 = data[10]
            t_FR4_64 = data[2]
            n=64
            x_64 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e64_L1_FR4 = (sum(abs(f(x_64)-data[10][0]))*(L/n))


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR4_CFL09m_RK3_MLPmod_1596.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_128 = data[7]
            z_FR4_128 = data[10]
            t_FR4_128 = data[2]
            n=128
            x_128 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e128_L1_FR4 = (sum(abs(f(x_128)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FR4_CFL09m_RK3_MLPmod_3232.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_256 = data[7]
            z_FR4_256 = data[10]
            t_FR4_256 = data[2]
            n=256
            x_256 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e256_L1_FR4 = (sum(abs(f(x_256)-data[10][0]))*(L/n))


        'MUSCL + LLF'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_MUSCL_CFL09_RK3_VA_31.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_8 = data[7]
            z_MUSCL_8 = data[10]
            n = 8
            e8_L1_MUSCL = (sum(abs(f(x_8)-data[10][0]))*(L/n))
            t_MUSCL_8 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_MUSCL_CFL09_RK3_VA_66.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_16 = data[7]
            z_MUSCL_16 = data[10]
            n = 16
            e16_L1_MUSCL = (sum(abs(f(x_16)-data[10][0]))*(L/n))
            t_MUSCL_16 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_MUSCL_CFL09_RK3_VA_137.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_32 = data[7]
            z_MUSCL_32 = data[10]
            n = 32
            e32_L1_MUSCL = (sum(abs(f(x_32)-data[10][0]))*(L/n))
            t_MUSCL_32 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_MUSCL_CFL09_RK3_VA_284.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_64 = data[7]
            z_MUSCL_64 = data[10]
            n = 64
            e64_L1_MUSCL = (sum(abs(f(x_64)-data[10][0]))*(L/n))
            t_MUSCL_64 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_MUSCL_CFL09_RK3_VA_575.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_128 = data[7]
            z_MUSCL_128 = data[10]
            n = 128
            e128_L1_MUSCL = (sum(abs(f(x_128)-data[10][0]))*(L/n))
            t_MUSCL_128 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_MUSCL_CFL09_RK3_VA_1154.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_256 = data[7]
            z_MUSCL_256 = data[10]
            n = 256
            e256_L1_MUSCL = (sum(abs(f(x_256)-data[10][0]))*(L/n))
            t_MUSCL_256 = data[2]

        'FR2_TVB'
        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FR2_CFL2m_35.npy', allow_pickle=True)
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FR2_CFL09m_RK3_TVD_34.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_TVB_8 = data[7]
            z_FR2_TVB_8 = data[10]
            t_FR2_TVB_8 = data[2]
            n=8
            x_8 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e8_L1_FR2_TVB = (sum(abs(f(x_8)-data[10][0]))*(L/n))


        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FR2_CFL2_33.npy', allow_pickle=True)
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FR2_CFL09m_RK3_TVD_69.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_TVB_16 = data[7]
            z_FR2_TVB_16 = data[10]
            t_FR2_TVB_16 = data[2]
            n=16
            x_16 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e16_L1_FR2_TVB = (sum(abs(f(x_16)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FR2_CFL09m_RK3_TVD_141.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_TVB_32 = data[7]
            z_FR2_TVB_32 = data[10]
            t_FR2_TVB_32 = data[2]
            n=32
            x_32 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e32_L1_FR2_TVB = (sum(abs(f(x_32)-data[10][0]))*(L/n))

        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FR2_CFL2_132.npy', allow_pickle=True)
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FR2_CFL09m_RK3_TVD_289.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_TVB_64 = data[7]
            z_FR2_TVB_64 = data[10]
            t_FR2_TVB_64 = data[2]
            n=64
            x_64 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e64_L1_FR2_TVB = (sum(abs(f(x_64)-data[10][0]))*(L/n))

        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR2_CFL2_262.npy', allow_pickle=True)
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FR2_CFL09m_RK3_TVD_584.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_TVB_128 = data[7]
            z_FR2_TVB_128 = data[10]
            t_FR2_TVB_128 = data[2]
            n=128
            x_128 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e128_L1_FR2_TVB = (sum(abs(f(x_128)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FR2_CFL09m_RK3_TVD_1170.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_TVB_256 = data[7]
            z_FR2_TVB_256 = data[10]
            t_FR2_TVB_256 = data[2]
            n=256
            x_256 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e256_L1_FR2_TVB = (sum(abs(f(x_256)-data[10][0]))*(L/n))

        'FOU'
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_8_FOU_CFL09_11.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_8 = data[6]
            Sg_FOU_8 = data[7]
            z_FOU_8 = data[10]
            n = 8
            e8_L1_FOU = (sum(abs(f(x_8)- data[10][0]))*(L/n))
            t_FOU_8 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_16_FOU_CFL09_21.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_16 = data[6]
            Sg_FOU_16 = data[7]
            z_FOU_16 = data[10]
            n = 16
            e16_L1_FOU = (sum(abs(f(x_16)-data[10][0]))*(L/n))
            t_FOU_16 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_32_FOU_CFL09_42.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_32 = data[6]
            Sg_FOU_32 = data[7]
            z_FOU_32 = data[10]
            n = 32
            e32_L1_FOU = (sum(abs(f(x_32)-data[10][0]))*(L/n))
            t_FOU_32 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_64_FOU_CFL09_85.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_64 = data[6]
            Sg_FOU_64 = data[7]
            z_FOU_64 = data[10]
            n = 64
            e64_L1_FOU = (sum(abs(f(x_64)-data[10][0]))*(L/n))
            t_FOU_64 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_128_FOU_CFL09_173.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_128 = data[6]
            Sg_FOU_128 = data[7]
            z_FOU_128 = data[10]
            n = 128
            e128_L1_FOU = (sum(abs(f(x_128)-data[10][0]))*(L/n))
            t_FOU_128 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FOU_CFL09_357.npy', allow_pickle=True)
        #datas = np.load('flying/results_case1_Moshiri_Manzari_5k_256_FOU_t_404.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_256 = data[6]
            Sg_FOU_256 = data[7]
            z_FOU_256 = data[10]
            n = 256
            e256_L1_FOU = (sum(abs(f(x_256)-data[10][0]))*(L/n))
            t_FOU_256 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_512_FOU_CFL09_726.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_512 = data[6]
            Sg_FOU_512 = data[7]
            z_FOU_512 = data[10]
            x_512 = np.linspace(0,1.5,len(So_FOU_512))
            n = 512
            e512_L1_FOU = (sum(abs(f(x_512)-data[10][0]))*(L/n))
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
        plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        plt.plot(x_128, z_FOU_128[1,:], '-ro')
        plt.plot(x_128, z_MUSCL_128[1,:], '-g*')
        #plt.plot(x_128, z_FR2_128[1,:], '-ys')
        plt.plot(x_128, z_FR2_TVB_128[1,:], '-bv')
        plt.plot(x_128, z_FR3_TVB_128[1,:], '-ys')
        plt.plot(x_128, z_FR4_TVB_128[1,:], '-m<')
        plt.legend(('MOC', 'FOU', 'MUSCL', 'FR-P1+TVB', 'FR-P2+TVB', 'FR-P3+TVB'))
        plt.grid()
        plt.xlim((1, 1.14))
        plt.ylim((0, 0.65))
        plt.title('NVCM example with 128x1x1 mesh')
        plt.ylabel('zC1')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zC1_FR_TVB_128.png')

        plt.figure(2)
        plt.plot(x_4000, zCO2_FOU_4000, 'k')
        plt.plot(x_128, z_FOU_128[0,:], '-r')
        plt.plot(x_128, z_MUSCL_128[0,:], '-g')
        #plt.plot(x_128, z_FR2_128[0,:], '-y')
        plt.plot(x_128, z_FR2_TVB_128[0,:], '-b')
        plt.plot(x_128, z_FR3_TVB_128[0,:], '-y')
        plt.plot(x_128, z_FR4_TVB_128[0,:], '-m')
        plt.legend(('MOC', 'FOU', 'MUSCL', 'FR-P1+TVB', 'FR-P2+TVB', 'FR-P3+TVB'))
        plt.grid()
        plt.xlim((0.55, 1.1))
        plt.ylim((0.75, 1))
        plt.title('NVCM example with 128x1x1 mesh')
        plt.ylabel('zCO2')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zCO2_FR_TVB_128.png')

        plt.figure(3)
        plt.plot(x_4000, Sg_FOU_4000, 'k')
        plt.plot(x_128, Sg_FOU_128, '-r')
        plt.plot(x_128, Sg_MUSCL_128, '-g')
        #plt.plot(x_128, Sg_FR2_128, '-y')
        plt.plot(x_128, Sg_FR2_TVB_128, '-b')
        plt.plot(x_128, Sg_FR3_TVB_128, '-y')
        plt.plot(x_128, Sg_FR4_TVB_128, '-m')
        plt.legend(('MOC', 'FOU', 'MUSCL', 'FR-P1+TVB', 'FR-P1+TVB', 'FR-P2+TVB', 'FR-P3+TVB'))
        plt.grid()
        plt.xlim((0.4, 1.2))
        #plt.ylim((0, 0.65))
        plt.title('NVCM example with 128x1x1 mesh')
        plt.ylabel('zCO2')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_Sg_FR_TVB_128.png')


        plt.figure(4)
        e_FR2 = np.log10(np.array([e8_L1_FR2, e16_L1_FR2, e32_L1_FR2, e64_L1_FR2, e128_L1_FR2, e256_L1_FR2]))
        e_FR3 = np.log10(np.array([e8_L1_FR3, e16_L1_FR3, e32_L1_FR3, e64_L1_FR3, e128_L1_FR3, e256_L1_FR3]))
        e_FR4 = np.log10(np.array([e8_L1_FR4, e16_L1_FR4, e32_L1_FR4, e64_L1_FR4, e128_L1_FR4, e256_L1_FR4]))
        e_MUSCL = np.log10(np.array([e8_L1_MUSCL, e16_L1_MUSCL, e32_L1_MUSCL, e64_L1_MUSCL, e128_L1_MUSCL, e256_L1_MUSCL]))
        e_FOU = np.log10(np.array([e8_L1_FOU, e16_L1_FOU, e32_L1_FOU, e64_L1_FOU, e128_L1_FOU, e256_L1_FOU]))

        x = np.log10(np.array([8,16,32,64,128,256]))
        y = -x

        plt.plot(x, e_FOU, '-m<', mfc='none')
        plt.plot(x, e_MUSCL, '-g*', mfc='none')
        plt.plot(x, e_FR2, '-ys', mfc='none')
        plt.plot(x, e_FR3, '-bv', mfc='none')
        plt.plot(x, e_FR4, '-ro', mfc='none')
        plt.plot(x,y,'-k')
        #plt.plot(x_256, z_FR2_256[1,:], '-gv', mfc='none')
        plt.legend(('FOU', 'MUSCL', 'FR-$\mathcal{P}_1$', 'FR-$\mathcal{P}_2$', 'FR-$\mathcal{P}_3$','1st order' ))
        plt.grid()
        plt.ylabel('log(|$E_{L1}$|)')
        plt.xlabel('log($n_b$)')
        plt.savefig('results/compositional/FR/5k_erro_zC1_FR_MLPmod_TVB2.png')

        plt.figure(5)
        x = np.log10(np.array([8,16,32,64,128,256]))
        x_FOU = np.log10(np.array([8,16,32,64,128,256, 512, 1024, 2048]))
        e_FR2 = np.array([e8_L1_FR2, e16_L1_FR2, e32_L1_FR2, e64_L1_FR2, e128_L1_FR2, e256_L1_FR2])
        e_FR2_TVB = np.array([e8_L1_FR2_TVB, e16_L1_FR2_TVB, e32_L1_FR2_TVB, e64_L1_FR2_TVB, e128_L1_FR2_TVB, e256_L1_FR2_TVB])
        e_FR3 = np.array([e8_L1_FR3, e16_L1_FR3, e32_L1_FR3, e64_L1_FR3, e128_L1_FR3, e256_L1_FR3])
        e_FR3_TVB = np.array([e8_L1_FR3_TVB, e16_L1_FR3_TVB, e32_L1_FR3_TVB, e64_L1_FR3_TVB, e128_L1_FR3_TVB, e256_L1_FR3_TVB])
        e_FR4 = np.array([e8_L1_FR4, e16_L1_FR4, e32_L1_FR4, e64_L1_FR4, e128_L1_FR4, e256_L1_FR4])
        e_FR4_TVB = np.array([e8_L1_FR4_TVB, e16_L1_FR4_TVB, e32_L1_FR4_TVB, e64_L1_FR4_TVB, e128_L1_FR4_TVB, e256_L1_FR4_TVB])

        e_MUSCL = np.array([e8_L1_MUSCL, e16_L1_MUSCL, e32_L1_MUSCL, e64_L1_MUSCL, e128_L1_MUSCL, e256_L1_MUSCL])
        e_FOU = np.array([e8_L1_FOU, e16_L1_FOU, e32_L1_FOU, e64_L1_FOU, e128_L1_FOU, e256_L1_FOU, e512_L1_FOU, e1024_L1_FOU, e2048_L1_FOU])
        plt.plot(x, np.log10(e_FOU[:6]), '-m>', mfc='none')
        plt.plot(x, np.log10(e_MUSCL), '-g*', mfc='none')
        #plt.plot(x, np.log10(e_FR2), '-ys', mfc='none')
        plt.plot(x, np.log10(e_FR2_TVB), '-yo', mfc='none')
        #plt.plot(x, np.log10(e_FR3), '-bs', mfc='none')
        plt.plot(x, np.log10(e_FR3_TVB), '-bv', mfc='none')
        #plt.plot(x, np.log10(e_FR4), '-rs', mfc='none')
        plt.plot(x, np.log10(e_FR4_TVB), '-rs', mfc='none')
        y = -x
        #plt.plot(x,y,'k')
        plt.legend(('FOU', 'MUSCL', 'FR-P1+hMLP', 'FR-P1+TVD', \
        'FR-P2+hMLP', 'FR-P2+TVD','FR-P3+hMLP', 'FR-P3+TVD'))
        plt.grid()
        plt.title('NVCM example')
        plt.ylabel('$log(|E_{L1}|)$')
        plt.xlabel('log($n_b$)')
        plt.savefig('results/compositional/FR/5k_EL1_comp_TVB.png')

        plt.figure(6)
        plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        plt.plot(x_32, z_FR2_32[1,:], '-ys')
        plt.plot(x_32, z_FR3_32[1,:], '-bv')
        plt.plot(x_32, z_FR4_32[1,:], '-ro')
        plt.plot(x_32, z_MUSCL_32[1,:], '-g*')
        plt.legend(('MOC', 'FR-P1', 'FR-P2', 'FR-P3','MUSCL'))
        plt.grid()
        #plt.xlim((1, 1.14))
        #plt.ylim((0, 0.65))
        plt.title('NVCM example with 32x1x1 mesh')
        plt.ylabel('zC1')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zC1_FR_TVB_32.png')
        import pdb; pdb.set_trace()
