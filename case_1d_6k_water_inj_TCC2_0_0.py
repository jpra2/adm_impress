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
        fx = open('x_incomp.txt','r')
        x_CMG = [float(line.rstrip('\n\r')) for line in fx]

        fP = open('P_incomp_low.txt','r')
        P_CMG = [float(line.rstrip('\n\r')) for line in fP]

        fSw = open('Sw_incomp_low.txt','r')
        Sw_CMG = [float(line.rstrip('\n\r')) for line in fSw]

        fSo = open('So_incomp_low.txt','r')
        So_CMG = [float(line.rstrip('\n\r')) for line in fSo]

        fSg = open('Sg_incomp_low.txt','r')
        Sg_CMG = [float(line.rstrip('\n\r')) for line in fSg]

        #datas = np.load('flying/results_water_inj_6k_modified_case_upw_4326.npy', allow_pickle=True)
        #datas = np.load('flying/results_water_inj_6k_200_FOU_modified_case_17445.npy', allow_pickle=True)
        datas = np.load('flying/results_water_inj_6k_200_upw_teste_777.npy', allow_pickle=True)
        
        for data in datas[2:]:
            Sw = data[5]
            So = data[6]
            Sg = data[7]
            Oil_p = data[8]
            Gas_p = data[9]
            pressure = data[4]/1e3
            time = data[3]
            n = 200

        datas = np.load('flying/results_water_inj_6k_200_MUSCL_MDW_modified_case_9412.npy', allow_pickle=True)

        for data in datas[2:]:
            Sw_MDW = data[5]
            So_MDW = data[6]
            Sg_MDW = data[7]
            Oil_p = data[8]
            Gas_p = data[9]
            pressure = data[4]/1e3
            time = data[3]
            n = 200
            x = np.linspace(2731.2*(1/(2*n)),2731.2*(1-1/(2*n)),n)

        x1=x
        sizeletter = 12
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300

        plt.rcParams.update({'font.size': sizeletter})
        plt.figure(1)
        #plt.title('t = 200 days')
        plt.plot(x_CMG, So_CMG, 'k', x1, So, '-ro', x, So_MDW, '-bs', mfc='none', markersize=5)
        plt.legend(('CMG', 'FOU-200', 'MUSCL-200'), prop={'size': sizeletter})
        plt.ylabel('Saturação de óleo')
        plt.xlabel('Distância [m]')
        plt.grid()
        plt.savefig('results/compositional/TCC2/saturation_oil_6k_MUSCL_inj' + '{}'.format(n) + '.png')

        plt.figure(2)
        #plt.title('t = 200 days')
        plt.plot(x_CMG, Sw_CMG, 'k', x1, Sw, '-ro', x, Sw_MDW, '-bs', mfc='none', markersize=5)
        plt.legend(('CMG', 'FOU-200', 'MUSCL-200'), prop={'size': sizeletter})
        plt.ylabel('Saturaço de água')
        plt.xlabel('Distância [m]')
        plt.grid()
        plt.savefig('results/compositional/TCC2/saturation_water_6k_MUSCL_inj' + '{}'.format(n) + '.png')

        plt.figure(3)
        #plt.title('Perfil de Saturação de Gás')
        plt.plot(x_CMG, Sg_CMG, 'k', x1, Sg, '-ro', x, Sg_MDW, '-bs', mfc='none', markersize=5)
        plt.legend(('CMG', 'FOU-200', 'MUSCL-200'), prop={'size': sizeletter})
        plt.ylabel('Saturação de gás')
        plt.xlabel('Distância [m]')
        plt.grid()
        plt.savefig('results/compositional/TCC2/saturation_gas_6k_MUSCL_inj' +'{}'.format(n) + '.png')

        plt.figure(4)
        #plt.title('t = 200 days')
        plt.plot(x_CMG, So_CMG, 'k', x1, So, '-ro', x, So_MDW, '-bs', mfc='none', markersize=5)
        plt.legend(('CMG', 'FOU-200', 'MUSCL-200'), prop={'size': sizeletter})
        plt.ylabel('Saturação de óleo')
        plt.xlabel('Distância [m]')
        plt.xlim(0,400)
        plt.ylim(0., 0.2)
        plt.grid()
        plt.savefig('results/compositional/TCC2/saturation_oil_6k_MUSCL_inj_zoom' + '{}'.format(n) + '.png')

        plt.figure(5)
        #plt.title('t = 200 days')
        plt.plot(x_CMG, Sw_CMG, 'k', x1, Sw, '-ro', x, Sw_MDW, '-bs', mfc='none', markersize=5)
        plt.legend(('CMG', 'FOU-200', 'MUSCL-200'), prop={'size': sizeletter})
        plt.ylabel('Saturaço de água')
        plt.xlabel('Distância [m]')
        plt.xlim(0,600)
        plt.ylim(0.2,1)
        plt.grid()
        plt.savefig('results/compositional/TCC2/saturation_water_6k_MUSCL_inj_zoom' + '{}'.format(n) + '.png')

        plt.figure(6)
        #plt.title('Perfil de Saturação de Gás')
        plt.plot(x_CMG, Sg_CMG, 'k', x1, Sg, '-ro', x, Sg_MDW, '-bs', mfc='none', markersize=5)
        plt.legend(('CMG', 'FOU-200', 'MUSCL-200'),prop={'size': sizeletter})
        plt.ylabel('Saturação de gás')
        plt.xlabel('Distância [m]')
        plt.xlim(0,1500)
        plt.grid()
        plt.savefig('results/compositional/TCC2/saturation_gas_6k_MUSCL_inj_zoom' +'{}'.format(n) + '.png')

        import pdb; pdb.set_trace()
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
