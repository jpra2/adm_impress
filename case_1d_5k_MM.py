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
        x_MOC = np.array([0, 0.0615711, 0.121019, 0.181529, 0.214437, 0.214968, \
            0.396497, 0.544586, 0.725584, 0.726115, 0.835987, 0.835987, 1.0552, \
            1.05573, 1.09395, 1.09395, 1.5])
        Sg_MOC = np.array([1, 1, 1, 1, 1, 0.947941, 0.947941, 0.947941, 0.947552, \
        0.877622, 0.877622, 0.772727, 0.734266, 0.563326, 0.563326, 0, 0])

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_500_FOU_954.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_500 = data[6]
            Sg_FOU_500 = data[7]
            z_FOU_500 = data[10]
            xkj_FOU_500 = data[13]

            xkj_FOU_500[:,1,Sg_FOU_500==0] = 0
            xkj_FOU_500[:,0,So_FOU_500==0] = 0
            x_500 = np.linspace(0,1.5,len(So_FOU_500))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_200_FOU_395.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_200 = data[6]
            Sg_FOU_200 = data[7]
            z_FOU_200 = data[10]
            xkj_FOU_200 = data[13]

            xkj_FOU_200[:,1,Sg_FOU_200==0] = 0
            xkj_FOU_200[:,0,So_FOU_200==0] = 0
            x_200 = np.linspace(0,1.5,len(So_FOU_200))

        datas = np.load('flying/results_case1_Moshiri_Manzari_5k_100_FOU_204.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_100 = data[6]
            Sg_FOU_100 = data[7]
            z_FOU_100 = data[10]
            xkj_FOU_100 = data[13]

            xkj_FOU_100[:,1,Sg_FOU_100==0] = 0
            xkj_FOU_100[:,0,So_FOU_100==0] = 0
            x_100 = np.linspace(0,1.5,len(So_FOU_100))

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
        import pdb; pdb.set_trace()
