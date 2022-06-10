import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
n = 8

for arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_water_inj_6k_8_FOU_modified_caset2_21.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[2:]:
            Sw_FOU = data[5]
            So_FOU = data[6]
            Sg_FOU = data[7]
            Oil_p_FOU = data[8]
            Gas_p_FOU = data[9]
            pressure_FOU = data[4]/1e3
            time_FOU = data[3]
            x1 = np.linspace(0.54624/n,2731.2/n*(n-1),n)

        datas2 = np.load('flying/results_water_inj_6k_8_FI_modified_caset2_69.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas2[2:]:
            Sw_FI = data[5]
            So_FI = data[6]
            Sg_FI = data[7]
            Oil_p_FI = data[8]
            Gas_p_FI = data[9]
            pressure_FI = data[4]/1e3

        datas3 = np.load('flying/results_water_inj_6k_8_FI_modified_caset2_CFL_0,5_NewSC_and_limits_21.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[2:]:
            Sw_FI_CFL05 = data[5]
            So_FI_CFL05 = data[6]
            Sg_FI_CFL05 = data[7]
            Oil_p_FI_CFL05 = data[8]
            Gas_p_FI_CFL05 = data[9]
            pressure_FI_CFL05 = data[4]/1e3


        #import pdb; pdb.set_trace()
        plt.figure(1)
        plt.title('t = 20 days - 8x1x1 mesh')
        plt.plot(x1, pressure_FOU, 'k')
        plt.plot(x1, pressure_FI, '-bo')
        plt.plot(x1, pressure_FI_CFL05, 'r')
        plt.ylabel('Pressure (kPa)')
        plt.xlabel('Distance')
        plt.legend(('FOU', 'FI', 'FI CFL 05'))
        plt.grid()
        plt.savefig('results/pressure_6k_inj' + '.png')

        plt.figure(2)
        plt.title('t = 20 days - 8x1x1 mesh')
        plt.plot(x1, So_FOU, 'k')
        plt.plot(x1, So_FI, '-bo')
        plt.plot(x1, So_FI_CFL05, 'r')
        plt.legend(('FOU', 'FI', 'FI CFL 05'))
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/saturation_oil_6k_inj' + '.png')

        plt.figure(3)
        plt.title('t = 20 days - 8x1x1 mesh')
        plt.plot(x1, Sw_FOU, 'k')
        plt.plot(x1, Sw_FI, '-bo')
        plt.plot(x1, Sw_FI_CFL05, 'r')
        plt.legend(('FOU', 'FI', 'FI CFL 05'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/saturation_water_6k_FOU_inj' + '.png')

        plt.figure(4)
        plt.title('t = 20 days - 8x1x1 mesh')
        plt.plot(x1, Sg_FOU, 'k')
        plt.plot(x1, Sg_FI, '-bo')
        plt.plot(x1, Sg_FI_CFL05, 'r')
        plt.legend(('FOU', 'FI', 'FI CFL 05'))
        plt.ylabel('Gas saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/saturation_gas_6k_FOU_inj' + '.png')

        print('Done')
        import pdb; pdb.set_trace()
