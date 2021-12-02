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
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_4000 = data[6]
            Sg_FOU_4000 = data[7]
            z_FOU_4000 = data[10]
            xkj_FOU_4000 = data[13]

            xkj_FOU_4000[:,1,Sg_FOU_4000==0] = 0
            xkj_FOU_4000[:,0,So_FOU_4000==0] = 0
            x_4000 = np.linspace(0,1.5,len(So_FOU_4000))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_200_MUSCL_LLF_t_1609.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_LLF_200 = data[6]
            Sg_MUSCL_LLF_200 = data[7]
            z_MUSCL_LLF_200 = data[10]
            xkj_MUSCL_LLF_200 = data[13]

            xkj_MUSCL_LLF_200[:,1,Sg_MUSCL_LLF_200==0] = 0
            xkj_MUSCL_LLF_200[:,0,So_MUSCL_LLF_200==0] = 0
            t_MUSCL_LLF_200 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_200_MUSCL_MDW_1440.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_MDW_200 = data[6]
            Sg_MUSCL_MDW_200 = data[7]
            z_MUSCL_MDW_200 = data[10]
            xkj_MUSCL_MDW_200 = data[13]

            xkj_MUSCL_MDW_200[:,1,Sg_MUSCL_MDW_200==0] = 0
            xkj_MUSCL_MDW_200[:,0,So_MUSCL_MDW_200==0] = 0
            t_MUSCL_MDW_200 = data[2]

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_200_MUSCL_DW_1331.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_DW_200 = data[6]
            Sg_MUSCL_DW_200 = data[7]
            z_MUSCL_DW_200 = data[10]
            xkj_MUSCL_DW_200 = data[13]

            xkj_MUSCL_DW_200[:,1,Sg_MUSCL_DW_200==0] = 0
            xkj_MUSCL_DW_200[:,0,So_MUSCL_DW_200==0] = 0
            t_MUSCL_DW_200 = data[2]


        '--------------------------------LLF-----------------------------------'



        x_200 = np.linspace(0+1.5/400,1.5-1.5/400,200)

        plt.figure(1)
        plt.plot(x_MOC, Sg_MOC, 'k')
        plt.plot(x_200, Sg_MUSCL_LLF_200, 'r')
        plt.plot(x_200, Sg_MUSCL_DW_200, 'g')
        plt.plot(x_200, Sg_MUSCL_MDW_200, 'm')
        #plt.plot(x_200, Sg_MUSCL_200, 'g')
        plt.legend(('MOC', 'MUSCL LLF-200', 'MUSCL DW-200', 'MUSCL MDW-200'))
        plt.grid()
        plt.xlim((0.5, 1.2))
        plt.ylim((0.4, 1))
        plt.title('NVCM example with 200x1x1 mesh')
        plt.ylabel('Sg')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/TCC2/5k_Sg_MUSCL_200.png')

        plt.figure(2)
        plt.plot(x_zC1_MOC, zC1_MOC, 'k')
        plt.plot(x_200, z_MUSCL_LLF_200[1,:], 'r')
        plt.plot(x_200, z_MUSCL_DW_200[1,:], 'g')
        plt.plot(x_200, z_MUSCL_MDW_200[1,:], 'm')
        #plt.plot(x_200, Sg_MUSCL_200, 'g')
        plt.legend(('MOC', 'MUSCL LLF-200', 'MUSCL DW-200', 'MUSCL MDW-200'))
        plt.grid()
        plt.xlim((1, 1.14))
        plt.ylim((0, 0.65))
        plt.title('NVCM example with 200x1x1 mesh')
        plt.ylabel('zC1')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/TCC2/5k_zC1_MUSCL_200.png')


        import pdb; pdb.set_trace()
