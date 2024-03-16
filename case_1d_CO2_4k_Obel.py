import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
################################################################################
''' Dados Comparativos'''
# with open('Qref_4kCO2_SC_CMG.txt', 'r') as q_cmg:
#     Tab_Qcmg = [[float(entry.replace(",",".")) for entry in line.split()] for line in q_cmg.readlines()]
# x_cmg = []
# y_cmg = []
# for i in Tab_Qcmg:
#     xc = i[0]
#     yc = i[1]
#     x_cmg.append(xc)
#     y_cmg.append(yc)
#
# with open('Qref_4kCO2_t1200.txt', 'r') as q_ref:
#     Tab_Qref = [[float(entry.replace(",",".")) for entry in line.split()] for line in q_ref.readlines()]
# x_ref = []
# y_ref = []
# for i in Tab_Qref:
#     xr = i[0]
#     yr = i[1]
#     x_ref.append(xr)
#     y_ref.append(yr)

###########simulation########################
for  arq in arquivos:
    if  arq.startswith(name):

        # datas = np.load('flying/results_4kCO2_finescale_q00001co2.npy', allow_pickle=True)
        # for data in datas[datas.shape[0]-1:]:
        #     time_co2 = data[22]
        #     time_co2 = np.asarray(time_co2)/86400
        #     oil_Production_rate_co2 = data[16]
        #     oil_Production_rate_co2 = np.asarray(oil_Production_rate_co2)*86400
        #
        # datas_w = np.load('flying/results_4kCO2_finescale_q00001w.npy', allow_pickle=True)
        # for data_w in datas_w[datas_w.shape[0]-1:]:
        #     time_w = data_w[22]
        #     time_w = np.asarray(time_w)/86400
        #     oil_Production_rate_w = data_w[16]
        #     oil_Production_rate_w = np.asarray(oil_Production_rate_w)*86400

#### WAG:flying/results_4kCO2_finescale_10001.npy
        datas_wag = np.load('flying/results_4kCO2_finescale_10001.npy', allow_pickle=True)
        #datas = np.load('flying/results_3k4ph_finescale_426.npy', allow_pickle=True)

        for data_wag in datas_wag[datas_wag.shape[0]-1:]:
            time_4k_wag = data_wag[22]
            time_4k_wag = np.asarray(time_4k_wag)/86400
            oil_Production_rate_wag = data_wag[16]
            oil_Production_rate_wag = np.asarray(oil_Production_rate_wag)*86400
################################################################################
        # plt.figure()
        # plt.plot(x_conv, y_conv, '.r')
        # plt.plot(x_ref, y_ref, 'r')
        # plt.grid()
        # plt.ylabel('Pressure')
        # plt.xlabel('Dimensionless distance')
        # plt.xlim(0,1800)
        # plt.ylim(0,15)
################################################################################
        plt.figure()
        #plt.plot(time_4k, oil_Production_rate)
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

        plt.savefig('resultados_WAG03/WAG4k_qw00012_tete2_w_wag.png')
