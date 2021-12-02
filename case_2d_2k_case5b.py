import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy import interpolate

flying = 'flying'
name = 'results_'
arquivos = os.listdir(flying)


for  arq in arquivos:
    if  arq.startswith(name):


        datas = np.load('flying/results_Hoteit_Firoo_2k_ex5b_IMPEC_FOU_4176.npy', allow_pickle=True)
        for data in datas[15:]:
            zC1 = data[10][0]
            vpi = data[1]
            print('vpi: ', vpi)
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[11]
            p0 = [0,1.,-1.0]
            p1 = [50.0,2.0,0.0]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            zC1_1 = zC1[ind_ans]

            p0 = [0,2.0,-1.0]
            p1 = [50.0,3,0.0]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            #ind_ans_sort = np.argsort(cent_mix)
            zC1_2 = zC1[ind_ans]
            #zC1_2 = zC1_2[ind_ans_sort]

            zC1_f_40 = 0.5*(zC1_1 + zC1_2)
            n=40
            x1_40 = np.linspace(0+50/(2*n),50-50/(2*n),n)
            #import pdb; pdb.set_trace()
            plt.figure(1)
            plt.title('FOU 40X20')
            plt.plot(x1_40, zC1_f_40, 'r')
            plt.grid()
            plt.ylabel('Mole Fraction of C1')
            plt.xlabel('Distance in X - Direction (m)')
            plt.savefig('results/compositional/TCC2/zC1_ex5b_FOU'  + '.png')
            import pdb; pdb.set_trace()
