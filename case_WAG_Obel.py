import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
################################################################################
''' Dados Comparativos'''
tabela_q00012 = pd.read_csv("gem_q11m3d.txt" , sep = "\t", header = None).to_numpy()
time_gem = tabela_q00012[:,0]
q_gem = tabela_q00012[:,1]

################################################################################
'''simulation'''
for  arq in arquivos:
    if  arq.startswith(name):

        datas_wag = np.load('flying/results_4kCO2_finescale_10001.npy', allow_pickle=True)

        for data_wag in datas_wag[datas_wag.shape[0]-1:]:
            time_4k_wag = data_wag[22]
            time_4k_wag = np.asarray(time_4k_wag)/86400
            oil_Production_rate_wag = data_wag[16]
            oil_Production_rate_wag = np.asarray(oil_Production_rate_wag)*86400

################################################################################
        ''' Plotagem '''
        plt.figure()
        plt.plot(time_gem, q_gem, 'k')
        plt.plot(time_4k_wag, oil_Production_rate_wag,'r')
        #plt.plot(time_co2, oil_Production_rate_co2, 'b')
        #plt.plot(time_w, oil_Production_rate_w, 'g')
        plt.grid()
        plt.title('Oil Production Rate')
        plt.ylabel('[m3/d]')
        plt.xlabel('time[day]')
        #plt.legend(['WAG', 'CO2', 'Water'])
        # plt.xlim(0,1200)
        # plt.ylim(0,15)

        plt.savefig('results_wag_compare/WAG4k_compare_q02.png')
