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
tabela_gem = pd.read_csv("gem_cum_rf_rfold.txt" , sep = "\t", header = None).to_numpy()
time_cum = tabela_gem[:,0]
cum_oil = tabela_gem[:,1]
rf = tabela_gem[:,2]
rf_old = tabela_gem[:,3]


#### Dados de OriginPro
with open('gem_So_c3.txt', 'r') as gem_so:
    Tab_gem_so = [[float(entry.replace(",",".")) for entry in line.split()] for line in gem_so.readlines()]
x_so = []
y_so = []
for i in Tab_gem_so:
    x_s = i[0]
    y_s = i[1]
    x_so.append(x_s)
    y_so.append(y_s)
################################################################################
'''simulation'''
for  arq in arquivos:
    if  arq.startswith(name):

        datas_wag = np.load('flying/results_c2_obel_finescale_9001.npy', allow_pickle=True)

        for data_wag in datas_wag[datas_wag.shape[0]-1:]:
            time_sim = data_wag[23]
            time_sim = np.asarray(time_sim)/86400
            vpi = data_wag[1]
            FRoil = data_wag[24]
            FRoil = np.asarray(FRoil)*100
            # pressure_well = data_wag[25]
            # pressure_well = np.asarray(pressure_well)/1000
            # time_simples = np.asarray(time_simples)/86400
            oil_Production_rate_RC = data_wag[16]
            oil_Production_rate_RC = np.asarray(oil_Production_rate_RC)*86400
            # oil_Production_rate_SC = data_wag[26]
            # oil_Production_rate_SC = np.asarray(oil_Production_rate_SC)*86400
            oil_Production = data_wag[18]
            oil_Production = np.asarray(oil_Production)
            #oil_Production_SC = data_wag[22]
            # oil_Production_SC = np.asarray(oil_Production_SC)

            Sw = data_wag[5]
            So = data_wag[6]
            Sg = data_wag[7]

            n = 200
            x_200 = np.linspace(2731.2*(1/(2*n)),2731.2*(1-1/(2*n)),n)

################################################################################
        ''' Plotagem '''
        x_200_pl = x_200
        plt.figure()
        #plt.plot(x_200_pl, Sw,'k')
        plt.plot(time_sim, oil_Production_rate_RC, 'k')
        #plt.plot(time_cum, rf, 'r')
        plt.grid()
        #plt.title('Saturação')
        #plt.title('Saturação de óleo')
        #plt.title('Oil Production')
        plt.ylabel('FR')
        plt.xlabel('tempo[dia]')
        #plt.legend(['PADMEC', 'GEM'])
        # plt.xlim(0,1200)
        # plt.ylim(0,15)

        plt.savefig('Results_Ed_C3/teste_c2_obel.png')
