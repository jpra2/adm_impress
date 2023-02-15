import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline
from packs.utils.utils_old_test import get_box

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
n = 40

x_ref = np.array([0, 22.7064, 30.7339, 35.0917, 42.6606, 50])
y_ref = np.array([1, 0.95, 0.65, 0.35, 0.05, 0])

for arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/Heterogeneo_plot/results_Hoteit_Firoo_2k_ex5a_40x20_FI_VPI0_43_PLOT_655.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            zC1_FI = data[10][0]
            x40 = np.linspace(0, 50, n)
            centroids = data[11]
            p0 = [0,2,-1.0]
            p1 = [50,3,0]
            ind_ans = get_box(centroids,np.array([p0,p1]))

            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)

            zC1_FI_1 = zC1_FI[ind_ans]
            zC1_FI_2 = zC1_FI_1[ind_ans_sort]


        #import pdb; pdb.set_trace()
        plt.figure(1)
        #plt.title('BL Sw - mesh refinement')
        plt.plot(x_ref, y_ref, '+b-')
        plt.plot(x40, zC1_FI_1, '--r')

        plt.legend(('Solução de referência', 'FI 40x20 CVs'))
        plt.ylabel('Fração molar global do metano em y=2,50m')
        plt.xlabel('Distância (m)')
        plt.grid()
        plt.savefig('results/caso5_40x20_' + '.png')

        import pdb; pdb.set_trace()
