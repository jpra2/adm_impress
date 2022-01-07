import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)

for arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_water_inj_6k_128_FOU_modified_case_1914.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[2:]:
            Sw = data[5]
            So = data[6]
            Sg = data[7]
            Oil_p = data[8]
            Gas_p = data[9]
            pressure = data[4]/1e3
            time = data[3]
            x1 = np.linspace(0.0,2725.7376,128)
        
        plt.figure(1)
        plt.title('t = 200 days - 128x1x1 mesh')
        plt.plot(x1, pressure)
        plt.ylabel('Pressure (kPa)')
        plt.xlabel('Distance')
        #plt.legend(('PADMEC-FOUM'))
        plt.grid()
        plt.savefig('results/pressure_6k_FOU_inj' + '.png')

        plt.figure(2)
        plt.title('t = 200 days - 128x1x1 mesh')
        plt.plot(x1, So)
        #plt.legend(('PADMEC-FOUM'))
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/saturation_oil_6k_FOU_inj' + '.png')

        plt.figure(3)
        plt.title('t = 200 days - 128x1x1 mesh')
        plt.plot(x1, Sw)
        #plt.legend(('PADMEC-FOUM'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/saturation_water_6k_FOU_inj' + '.png')

        plt.figure(4)
        plt.title('t = 200 days - 128x1x1 mesh')
        plt.plot(x1, Sg)
        #plt.legend(('PADMEC-FOUM'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/saturation_gas_6k_FOU_inj' + '.png')

        print('Done')
        import pdb; pdb.set_trace()
