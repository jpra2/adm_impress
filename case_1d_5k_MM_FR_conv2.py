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
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_20_FR2_CFL09m_RK3_MLPmod_90.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_20 = data[7]
            z_FR2_20 = data[10]
            t_FR2_20 = data[2]
            n=20
            x_20 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e20_L1_FR2 = (sum(abs(f(x_20)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_40_FR2_CFL09m_RK3_MLPmod_183.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_40 = data[7]
            z_FR2_40 = data[10]
            t_FR2_40 = data[2]
            n=40
            x_40 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e40_L1_FR2 = (sum(abs(f(x_40)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_80_FR2_CFL09m_RK3_MLPmod_369.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_80 = data[7]
            z_FR2_80 = data[10]
            t_FR2_80 = data[2]
            n=80
            x_80 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e80_L1_FR2 = (sum(abs(f(x_80)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_160_FR2_CFL09_RK3_MLPmod_734.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_160 = data[7]
            z_FR2_160 = data[10]
            t_FR2_160 = data[2]
            n=160
            x_160 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e160_L1_FR2 = (sum(abs(f(x_160)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_320_FR2_CFL09_RK3_MLPmod_1511.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_320 = data[7]
            z_FR2_320 = data[10]
            t_FR2_320 = data[2]
            n=320
            x_320 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e320_L1_FR2 = (sum(abs(f(x_320)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_640_FR2_CFL09m_RK3_MLPmod_4027.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR2_640 = data[7]
            z_FR2_640 = data[10]
            t_FR2_640 = data[2]
            n=640
            x_640 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e640_L1_FR2 = (sum(abs(f(x_640)-data[10][0]))*(L/n))


        'FR3'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_20_FR3_CFL09m_RK3_MLPmod_49.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_20 = data[7]
            z_FR3_20 = data[10]
            t_FR3_20 = data[2]
            n=20
            x_20 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e20_L1_FR3 = (sum(abs(f(x_20)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_40_FR3_CFL09m_RK3_MLPmod_100.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_40 = data[7]
            z_FR3_40 = data[10]
            t_FR3_40 = data[2]
            n=40
            x_40 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e40_L1_FR3 = (sum(abs(f(x_40)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_80_FR3_CFL09m_RK3_MLPmod_207.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_80 = data[7]
            z_FR3_80 = data[10]
            t_FR3_80 = data[2]
            n=80
            x_80 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e80_L1_FR3 = (sum(abs(f(x_80)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_160_FR3_CFL09m_RK3_MLPmod_512.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_160 = data[7]
            z_FR3_160 = data[10]
            t_FR3_160 = data[2]
            n=160
            x_160 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e160_L1_FR3 = (sum(abs(f(x_160)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_320_FR3_CFL09m_RK3_MLPmod_1223.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_320 = data[7]
            z_FR3_320 = data[10]
            t_FR3_320 = data[2]
            n=320
            x_320 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e320_L1_FR3 = (sum(abs(f(x_320)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_640_FR3_CFL09m_RK3_MLPmod_2612.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR3_640 = data[7]
            z_FR3_640 = data[10]
            t_FR3_640 = data[2]
            n=640
            x_640 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e640_L1_FR3 = (sum(abs(f(x_640)-data[10][0]))*(L/n))


        'FR4'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_20_FR4_CFL09m_RK3_MLPmod_60.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_20 = data[7]
            z_FR4_20 = data[10]
            t_FR4_20 = data[2]
            n=20
            x_20 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e20_L1_FR4 = (sum(abs(f(x_20)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_40_FR4_CFL09m_RK3_MLPmod_122.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_40 = data[7]
            z_FR4_40 = data[10]
            t_FR4_40 = data[2]
            n=40
            x_40 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e40_L1_FR4 = (sum(abs(f(x_40)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_80_FR4_CFL09m_RK3_MLPmod_278.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_80 = data[7]
            z_FR4_80 = data[10]
            t_FR4_80 = data[2]
            n=80
            x_80 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e80_L1_FR4 = (sum(abs(f(x_80)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_160_FR4_CFL09m_RK3_MLPmod_640.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_160 = data[7]
            z_FR4_160 = data[10]
            t_FR4_160 = data[2]
            n=160
            x_160 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e160_L1_FR4 = (sum(abs(f(x_160)-data[10][0]))*(L/n))


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_320_FR4_CFL09m_RK3_MLPmod_1596.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_320 = data[7]
            z_FR4_320 = data[10]
            t_FR4_320 = data[2]
            n=320
            x_320 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e320_L1_FR4 = (sum(abs(f(x_320)-data[10][0]))*(L/n))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_640_FR4_CFL09m_RK3_MLPmod_3232.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_FR4_640 = data[7]
            z_FR4_640 = data[10]
            t_FR4_640 = data[2]
            n=640
            x_640 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)
            e640_L1_FR4 = (sum(abs(f(x_640)-data[10][0]))*(L/n))


        'MUSCL + LLF'

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_20_MUSCL_CFL09_RK3_VA_31.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_20 = data[7]
            z_MUSCL_20 = data[10]
            n = 20
            e20_L1_MUSCL = (sum(abs(f(x_20)-data[10][0]))*(L/n))
            t_MUSCL_20 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_40_MUSCL_CFL09_RK3_VA_66.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_40 = data[7]
            z_MUSCL_40 = data[10]
            n = 40
            e40_L1_MUSCL = (sum(abs(f(x_40)-data[10][0]))*(L/n))
            t_MUSCL_40 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_80_MUSCL_CFL09_RK3_VA_137.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_80 = data[7]
            z_MUSCL_80 = data[10]
            n = 80
            e80_L1_MUSCL = (sum(abs(f(x_80)-data[10][0]))*(L/n))
            t_MUSCL_80 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_160_MUSCL_CFL09_RK3_VA_284.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_160 = data[7]
            z_MUSCL_160 = data[10]
            n = 160
            e160_L1_MUSCL = (sum(abs(f(x_160)-data[10][0]))*(L/n))
            t_MUSCL_160 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_320_MUSCL_CFL09_RK3_VA_575.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_320 = data[7]
            z_MUSCL_320 = data[10]
            n = 320
            e320_L1_MUSCL = (sum(abs(f(x_320)-data[10][0]))*(L/n))
            t_MUSCL_320 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_640_MUSCL_CFL09_RK3_VA_1154.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sg_MUSCL_640 = data[7]
            z_MUSCL_640 = data[10]
            n = 640
            e640_L1_MUSCL = (sum(abs(f(x_640)-data[10][0]))*(L/n))
            t_MUSCL_640 = data[2]



        'FOU'
        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_20_FOU_CFL09_11.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_20 = data[6]
            Sg_FOU_20 = data[7]
            z_FOU_20 = data[10]
            n = 20
            e20_L1_FOU = (sum(abs(f(x_20)- data[10][0]))*(L/n))
            t_FOU_20 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_40_FOU_CFL09_21.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_40 = data[6]
            Sg_FOU_40 = data[7]
            z_FOU_40 = data[10]
            n = 40
            e40_L1_FOU = (sum(abs(f(x_40)-data[10][0]))*(L/n))
            t_FOU_40 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_80_FOU_CFL09_42.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_80 = data[6]
            Sg_FOU_80 = data[7]
            z_FOU_80 = data[10]
            n = 80
            e80_L1_FOU = (sum(abs(f(x_80)-data[10][0]))*(L/n))
            t_FOU_80 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_160_FOU_CFL09_85.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_160 = data[6]
            Sg_FOU_160 = data[7]
            z_FOU_160 = data[10]
            n = 160
            e160_L1_FOU = (sum(abs(f(x_160)-data[10][0]))*(L/n))
            t_FOU_160 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_320_FOU_CFL09_173.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_320 = data[6]
            Sg_FOU_320 = data[7]
            z_FOU_320 = data[10]
            n = 320
            e320_L1_FOU = (sum(abs(f(x_320)-data[10][0]))*(L/n))
            t_FOU_320 = data[2]


        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_640_FOU_CFL09_357.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_640 = data[6]
            Sg_FOU_640 = data[7]
            z_FOU_640 = data[10]
            n = 640
            e640_L1_FOU = (sum(abs(f(x_640)-data[10][0]))*(L/n))
            t_FOU_640 = data[2]


        plt.figure(1)
        plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        plt.plot(x_128, z_FR2_128[1,:], '-ys')
        plt.plot(x_128, z_FR3_128[1,:], '-bv')
        plt.plot(x_128, z_FR4_128[1,:], '-ro')
        plt.plot(x_128, z_MUSCL_128[1,:], '-g*')
        plt.legend(('MOC', 'FR-P1', 'FR-P2', 'FR-P3','MUSCL'))
        plt.grid()
        plt.xlim((1, 1.14))
        plt.ylim((0, 0.65))
        plt.title('NVCM example with 128x1x1 mesh')
        plt.ylabel('zC1')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zC1_FR_MLPmod_128.png')


        plt.figure(2)
        plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        plt.plot(x_64, z_FR2_64[1,:], '-ys', mfc='none')
        plt.plot(x_64, z_FR3_64[1,:], '-bv', mfc='none')
        plt.plot(x_64, z_FR4_64[1,:], '-ro', mfc='none')
        plt.plot(x_64, z_MUSCL_64[1,:], '-g*', mfc='none')
        plt.xlim((1, 1.5))
        plt.ylim((0, 0.65))
        plt.legend(('MOC', 'FR-P1', 'FR-P2','FR-P3', 'MUSCL'))
        plt.grid()
        plt.title('NVCM example with 64x1x1 mesh')
        plt.ylabel('zC1')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zC1_FR_MLPmod_64.png')


        plt.figure(3)
        plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        #plt.plot(x_256, z_FR2_256[1,:], '-gv', mfc='none')
        plt.plot(x_256, z_FR2_256[1,:], '-ys', mfc='none')
        plt.plot(x_256, z_FR3_256[1,:], '-bv', mfc='none')
        plt.plot(x_256, z_FR4_256[1,:], '-ro', mfc='none')
        plt.plot(x_256, z_MUSCL_256[1,:], '-g*', mfc='none')
        plt.legend(('MOC','FR-P1', 'FR-P2', 'FR-P3', 'MUSCL'))
        plt.grid()
        plt.xlim((1, 1.14))
        plt.ylim((0, 0.65))
        plt.title('NVCM example with 256x1x1 mesh')
        plt.ylabel('zC1')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/FR/5k_zC1_FR_MLPmod_256.png')

        plt.figure(4)
        e_FR2 = np.log10(np.array([e8_L1_FR2, e16_L1_FR2, e32_L1_FR2, e64_L1_FR2, e128_L1_FR2, e256_L1_FR2]))
        e_FR3 = np.log10(np.array([e8_L1_FR3, e16_L1_FR3, e32_L1_FR3, e64_L1_FR3, e128_L1_FR3, e256_L1_FR3]))
        e_FR4 = np.log10(np.array([e8_L1_FR4, e16_L1_FR4, e32_L1_FR4, e64_L1_FR4, e128_L1_FR4, e256_L1_FR4]))
        e_MUSCL = np.log10(np.array([e8_L1_MUSCL, e16_L1_MUSCL, e32_L1_MUSCL, e64_L1_MUSCL, e128_L1_MUSCL, e256_L1_MUSCL]))
        e_FOU = np.log10(np.array([e8_L1_FOU, e16_L1_FOU, e32_L1_FOU, e64_L1_FOU, e128_L1_FOU, e256_L1_FOU]))

        x = np.log10(np.array([8,16,32,64,128,256]))
        plt.plot(x, e_FOU, '-m<', mfc='none')
        plt.plot(x, e_MUSCL, '-g*', mfc='none')
        plt.plot(x, e_FR2, '-ys', mfc='none')
        plt.plot(x, e_FR3, '-bv', mfc='none')
        plt.plot(x, e_FR4, '-ro', mfc='none')
        #plt.plot(x_256, z_FR2_256[1,:], '-gv', mfc='none')
        plt.legend(('FOU', 'MUSCL', 'FR-$\mathcal{P}_1$', 'FR-$\mathcal{P}_2$', 'FR-$\mathcal{P}_3$' ))
        plt.grid()
        plt.ylabel('log(|$E_{L1}$|)')
        plt.xlabel('log($n_b$)')
        plt.savefig('results/compositional/FR/5k_erro_zC1_FR_MLPmod.png')


        plt.figure(5)
        x = np.array([8,16,32,64,128,256])
        t_FR2 = np.array([t_FR2_8, t_FR2_16, t_FR2_32, t_FR2_64, t_FR2_128, t_FR2_256])
        t_FR3 = np.array([t_FR3_8, t_FR3_16, t_FR3_32, t_FR3_64, t_FR3_128, t_FR3_256])
        t_FR4 = np.array([t_FR4_8, t_FR4_16, t_FR4_32, t_FR4_64, t_FR4_128, t_FR4_256])
        t_MUSCL = np.array([t_MUSCL_8, t_MUSCL_16, t_MUSCL_32, t_MUSCL_64, t_MUSCL_128, t_MUSCL_256])
        t_FOU = np.array([t_FOU_8, t_FOU_16, t_FOU_32, t_FOU_64, t_FOU_128, t_FOU_256, t_FOU_512, t_FOU_1024, t_FOU_2048])

        e_FR2 = np.array([e8_L1_FR2, e16_L1_FR2, e32_L1_FR2, e64_L1_FR2, e128_L1_FR2, e256_L1_FR2])
        e_FR3 = np.array([e8_L1_FR3, e16_L1_FR3, e32_L1_FR3, e64_L1_FR3, e128_L1_FR3, e256_L1_FR3])
        e_FR4 = np.array([e8_L1_FR4, e16_L1_FR4, e32_L1_FR4, e64_L1_FR4, e128_L1_FR4, e256_L1_FR4])
        e_MUSCL = np.array([e8_L1_MUSCL, e16_L1_MUSCL, e32_L1_MUSCL, e64_L1_MUSCL, e128_L1_MUSCL, e256_L1_MUSCL])
        e_FOU = np.array([e8_L1_FOU, e16_L1_FOU, e32_L1_FOU, e64_L1_FOU, e128_L1_FOU, e256_L1_FOU, e512_L1_FOU, e1024_L1_FOU, e2048_L1_FOU])

        plt.plot(t_FR2, e_FR2, '-ys', mfc='none')
        plt.plot(t_FR3, e_FR3, '-bo', mfc='none')
        plt.plot(t_FR4, e_FR4, '-rv', mfc='none')
        plt.plot(t_MUSCL, e_MUSCL, '-g*', mfc='none')
        plt.plot(t_FOU, e_FOU, '-m>', mfc='none')

        plt.legend(('FR-P1', 'FR-P2', 'FR-P3', 'MUSCL', 'FOU'))
        plt.grid()
        plt.xscale('log')
        plt.yscale('log')
        plt.title('NVCM example')
        plt.ylabel('$|E_{L1}|$')
        plt.xlabel('time[s]')
        plt.savefig('results/compositional/FR/5k_error_time_comp_MLPmod.png')

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
        plt.savefig('results/compositional/FR/5k_zC1_FR_MLPmod_32.png')
        import pdb; pdb.set_trace()
