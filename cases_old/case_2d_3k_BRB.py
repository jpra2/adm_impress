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
        t_BRB = np.array([0, 0.73583, 1.83959, 6.43856, 9.3819, 13.7969, 38.0759, \
        157.837, 247.057, 301.876, 311.442, 330.758, 336.829, 344.371, 353.753, \
        368.469, 395.143, 413.723, 456.218, 498.528, 500])
        oil_prod_BRB = np.array([0, 325.084, 815.162, 1084.5, 1101, 1111.26, 1115.72, \
        1121.52, 1124.19, 1127.31, 1132.22, 1130.43, 1048.83, 914.158, 824.526, \
        748.718, 675.139, 643.478, 576.589, 528.428, 526.198])

        t_BRBg = np.array([0, 2.72331, 6.17284, 13.7981, 305.374, 330.065, 335.875, \
        337.863, 343.863, 350.036, 358.206, 373.275, 384.35, 409.768, 443.537, \
        479.121, 494.19, 500])
        gas_prod_BRB = np.array([0, 15.6495, 21.4085, 22.41, 22.6604, 24.4131, 39.3114, \
        86.385, 170.892, 186.291, 201.941, 224.726, 236.745, 261.158, 285.696, \
        306.729, 315.869, 320]) * 1e3

        datas = np.load('flying/results_case1_2d_BRB_IMPECn_1460.npy', allow_pickle=True)
        #datas = np.load('flying/results_case1_2d_BRB_IMPECf_1427.npy', allow_pickle=True)
        #datas = np.load('flying/results_Hoteit_Firoo_3k_500_upw2_3016.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            t = data[3]
            z = data[10]
            So = data[6]
            Sg = data[7]
            oil_prod_IMPEC_20x20 = np.array(data[14])*(24*60*60)
            gas_prod_IMPEC_20x20 = np.array(data[15])*(24*60*60)
            t_prod_IMPEC_20x20 = np.array(data[16])/(24*60*60)
            #t_prod_20x20 = np.linspace(0,500,147)

        datas = np.load('flying/results_case1_2d_BRB_40x40_IMPEC_6437.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            t = data[3]
            z = data[10]
            So = data[6]
            Sg = data[7]
            oil_prod_IMPEC_40x40 = np.array(data[14])*(24*60*60)
            gas_prod_IMPEC_40x40 = np.array(data[15])*(24*60*60)
            t_prod_IMPEC_40x40 = np.array(data[16])/(24*60*60)
            import pdb; pdb.set_trace()


        datas = np.load('flying/results_case1_2d_BRB_IMPSAT_1463.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            t = data[3]
            z = data[10]
            So = data[6]
            Sg = data[7]
            oil_prod_IMPSAT_20x20 = np.array(data[14])*(24*60*60)
            gas_prod_IMPSAT_20x20 = np.array(data[15])*(24*60*60)
            t_prod_IMPSAT_20x20 = np.array(data[16])/(24*60*60)
            #t_prod_20x20 = np.linspace(0,500,147)

        plt.figure(1)
        plt.plot(t_prod_IMPEC_20x20, oil_prod_IMPEC_20x20, 'r')
        plt.plot(t_prod_IMPEC_40x40, oil_prod_IMPEC_40x40, 'g')
        plt.plot(t_prod_IMPSAT_20x20, oil_prod_IMPSAT_20x20, 'b')
        plt.plot(t_BRB, oil_prod_BRB, 'k')
        plt.ylabel('Oil production rate [m$^3$/day]')
        plt.legend(( 'IMPEC-20x20', 'IMPEC 40x40', 'IMPSAT-20x20', 'Reference Fernandes 80x80'))
        plt.xlabel('time [days]')
        plt.title('Case1 - Fernandes\'s Dissertation')
        plt.savefig('results/compositional/case1BRB_oil_recov_IMPEC.png')

        plt.figure(2)
        plt.plot(t_prod_IMPEC_20x20, gas_prod_IMPEC_20x20/1e3, 'r')
        plt.plot(t_prod_IMPEC_40x40, gas_prod_IMPEC_40x40/1e3, 'r')
        plt.plot(t_prod_IMPSAT_20x20, gas_prod_IMPSAT_20x20/1e3, 'b')
        plt.plot(t_BRBg, gas_prod_BRB/1e3, 'k')
        plt.ylabel('Gas production rate [10$^3$m$^3$/day]')
        plt.legend(( 'IMPEC-20x20', 'IMPEC 40x40', 'IMPSAT-20x20', 'Reference Fernandes 80x80'))
        plt.xlabel('time [days]')
        plt.title('Case1 - Fernandes\'s Dissertation')
        plt.savefig('results/compositional/case1BRB_gas_recov_IMPEC.png')
        import pdb; pdb.set_trace()
