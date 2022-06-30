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
        fx = open('x1024.txt','r')
        x_CMG = [float(line.rstrip('\n\r')) for line in fx]

        fP = open('P1024.txt','r')
        P_CMG = [float(line.rstrip('\n\r')) for line in fP]

        fSw = open('Sw1024.txt','r')
        Sw_CMG = [float(line.rstrip('\n\r')) for line in fSw]

        fSo = open('So1024.txt','r')
        So_CMG = [float(line.rstrip('\n\r')) for line in fSo]

        fSg = open('Sg1024.txt','r')
        Sg_CMG = [float(line.rstrip('\n\r')) for line in fSg]
        n=128

        #datas = np.load('flying/results_water_inj_6k_modified_case_upw_4326.npy', allow_pickle=True)
        datas = np.load('flying/results_water_inj_6k_128_case_4520.npy', allow_pickle=True)
        import pdb; pdb.set_trace()
        for data in datas[2:]:
            Sw = data[5]
            So = data[6]
            Sg = data[7]
            Oil_p = data[8]
            Gas_p = data[9]
            pressure = data[4]/1e3
            time = data[3]
            #x1 = np.linspace(0.54624,2725.7376,500)
        datas = np.load('flying/results_water_inj_6k_128_MUSCL_case_8773.npy', allow_pickle=True)
        for data in datas[2:]:
            Sw_MUSCL = data[5]
            So_MUSCL = data[6]
            Sg_MUSCL = data[7]
            Oil_p_MUSCL = data[8]
            Gas_p_MUSCL = data[9]
            pressure_MUSCL = data[4]/1e3
            time_MUSCL = data[3]

        x = np.linspace(0.54624/n,2731.2/n*(n-1),n)
        x1 = x
        plt.figure(1)
        plt.title('t = 200 days - 8x1x1 mesh')
        plt.plot(x1, pressure, 'k', x_CMG, P_CMG, 'r', x, pressure_MUSCL, 'y')
        plt.ylabel('Pressure (kPa)')
        plt.xlabel('Distance')
        plt.legend(('PADMEC-FOUM', 'CMG', 'PADMEC-MUSCL'))
        plt.grid()
        plt.savefig('results/compositional/MUSCL_tests/pressure_6k_MUSCL_inj' + '{}'.format(n) + '.png')

        plt.figure(2)
        plt.title('t = 200 days')
        plt.plot(x1, So, 'k', x_CMG, So_CMG, 'r', x, So_MUSCL, 'y')
        plt.legend(('PADMEC-FOUM', 'CMG', 'PADMEC-MUSCL'))
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/compositional/MUSCL_tests/saturation_oil_6k_MUSCL_inj' + '{}'.format(n) + '.png')

        plt.figure(3)
        plt.title('t = 200 days')
        plt.plot(x1, Sw, 'k', x_CMG, Sw_CMG, 'r', x, Sw_MUSCL, 'y')
        plt.legend(('PADMEC-FOUM', 'CMG', 'PADMEC-MUSCL'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/compositional/MUSCL_tests/saturation_water_6k_MUSCL_inj' + '{}'.format(n) + '.png')

        plt.figure(4)
        plt.title('Perfil de Saturação de Gás')
        plt.plot(x1, Sg, 'k', x_CMG, Sg_CMG, 'r', x, Sg_MUSCL, 'y')
        plt.legend(('FOU-128', 'CMG', 'MUSCL-128'))
        plt.ylabel('Saturação de gás ')
        plt.xlabel('Distância')
        plt.grid()
        plt.savefig('results/compositional/TCC2/saturation_gas_6k_MUSCL_inj' +'{}'.format(n) + '.png')

        '''x = np.linspace(0.54624,2725.7376,500)
        plt.figure(1)
        plt.title('t = 200 days')
        plt.plot(x, pressure, 'k', x_CMG, P_CMG, 'r', x, pressure_MUSCL, 'y')
        plt.ylabel('Pressure (kPa)')
        plt.xlabel('Distance')
        plt.legend(('PADMEC-FOUM', 'CMG', 'PADMEC-MUSCL'))
        plt.grid()
        plt.savefig('results/compositional/pressure_6k_MUSCL_inj' + '.png')

        plt.figure(2)
        plt.title('t = 200 days')
        plt.plot(x, So, 'k', x_CMG, So_CMG, 'r', x, So_MUSCL, 'y')
        plt.legend(('PADMEC-FOUM', 'CMG', 'PADMEC-MUSCL'))
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/compositional/saturation_oil_6k_MUSCL_inj' + '.png')

        plt.figure(3)
        plt.title('t = 200 days')
        plt.plot(x, Sw, 'k', x_CMG, Sw_CMG, 'r', x, Sw_MUSCL, 'y')
        plt.legend(('PADMEC-FOUM', 'CMG', 'PADMEC-MUSCL'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/compositional/saturation_water_6k_MUSCL_inj' + '.png')'''
        '''plt.figure(1)
        plt.title('t = 200 days')
        plt.plot(x, pressure, 'k', x_CMG, P_CMG, 'r')
        plt.ylabel('Pressure (kPa)')
        plt.xlabel('Distance')
        plt.legend(('PADMEC', 'CMG'))
        plt.grid()
        plt.savefig('results/compositional/pressure_6k_inj' + '.png')

        plt.figure(2)
        plt.title('t = 200 days')
        plt.plot(x, So, 'k', x_CMG, So_CMG, 'r')
        plt.legend(('PADMEC', 'CMG'))
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/compositional/saturation_oil_6k_inj' + '.png')

        plt.figure(3)
        plt.title('t = 200 days')
        plt.plot(x, Sw, 'k', x_CMG, Sw_CMG, 'r')
        plt.legend(('PADMEC', 'CMG'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/compositional/saturation_water_6k_inj' + '.png')'''

        import pdb; pdb.set_trace()
