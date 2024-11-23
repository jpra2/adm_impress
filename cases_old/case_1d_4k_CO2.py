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


        datas = np.load('flying/results_4kCO2_finescale_10001.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            oil_prod_rate = data[14]
            cum_oil_prod = data[16]
            time = data[17]

        plt.figure(1)
        plt.plot(time,cum_oil_prod, '-bo', mfc='none',markersize=3)
        plt.grid()
        plt.title('Cum oil production')
        plt.ylabel('m3')
        plt.xlabel('time[s]')
        plt.savefig('results/compositional/oil_prod_4k_CO2.png')

        plt.figure(2)
        plt.plot(np.array(time)/86400, np.array(oil_prod_rate)*86400, 'sr', mfc='none',markersize=3)
        plt.grid()
        plt.title('Oil production rate')
        plt.ylabel('m3/day')
        plt.xlabel('time[day]')
        #plt.ylim((1200,1400))
        plt.savefig('results/compositional/oil_prod_rate_4k_CO2.png')

        import pdb; pdb.set_trace()
