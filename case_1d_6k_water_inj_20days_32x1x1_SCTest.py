import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
n = 32

fx = open('x_points_CMG_SchmallNEW.txt','r')
x_CMG = [float(line.rstrip('\n\r')) for line in fx]

#fP = open('P1024.txt','r')
#P_CMG = [float(line.rstrip('\n\r')) for line in fP]

fSw = open('Sw_points_CMG_SchmallNEW.txt','r')
Sw_CMG = [float(line.rstrip('\n\r')) for line in fSw]

fSo = open('So_points_CMG_SchmallNEW.txt','r')
So_CMG = [float(line.rstrip('\n\r')) for line in fSo]

fSg = open('Sg_points_CMG_SchmallNEW.txt','r')
Sg_CMG = [float(line.rstrip('\n\r')) for line in fSg]


for arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/Teste_criterio_de_parada/32_CV/results_6k_SchmallNEW_32_FI_20days_ADIMENSIONALIZADO_e6_202.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FI_e6 = data[5]
            So_FI_e6 = data[6]
            Sg_FI_e6 = data[7]
            Oil_p_FI_e6 = data[8]
            Gas_p_FI_e6 = data[9]
            pressure_FI_e6 = data[4]/1e3
            time_FI_e6 = data[3]
            zC1_FI_e6 = data[10][0]
            zC3_FI_e6 = data[10][1]
            zC6_FI_e6 = data[10][2]
            zC10_FI_e6 = data[10][3]
            zC15_FI_e6 = data[10][4]
            zC20_FI_e6 = data[10][5]
            x1 = np.linspace(0, 50, n)

        datas2 = np.load('flying/Teste_criterio_de_parada/32_CV/results_6k_SchmallNEW_32_FI_20days_ADIMENSIONALIZADO_e3_NEWDT_50.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas2[1:]:
            Sw_FI_e3 = data[5]
            So_FI_e3 = data[6]
            Sg_FI_e3 = data[7]
            Oil_p_FI_e3 = data[8]
            Gas_p_FI_e3 = data[9]
            pressure_FI_e3 = data[4]/1e3
            time_FI_e3 = data[3]
            zC1_FI_e3 = data[10][0]
            zC3_FI_e3 = data[10][1]
            zC6_FI_e3 = data[10][2]
            zC10_FI_e3 = data[10][3]
            zC15_FI_e3 = data[10][4]
            zC20_FI_e3 = data[10][5]

        datas3 = np.load('flying/Teste_criterio_de_parada/32_CV/results_6k_SchmallNEW_32_FI_20days_ADIMENSIONALIZADO_e2_2002.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI_e2 = data[5]
            So_FI_e2 = data[6]
            Sg_FI_e2 = data[7]
            Oil_p_FI_e2 = data[8]
            Gas_p_FI_e2 = data[9]
            pressure_FI_e2 = data[4]/1e3
            time_FI_e2 = data[3]
            zC1_FI_e2 = data[10][0]
            zC3_FI_e2 = data[10][1]
            zC6_FI_e2 = data[10][2]
            zC10_FI_e2 = data[10][3]
            zC15_FI_e2 = data[10][4]
            zC20_FI_e2 = data[10][5]

        """datas4 = np.load('flying/Teste_criterio_de_parada/32_CV/.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas4[1:]:
            Sw_FI_e1 = data[5]
            So_FI_e1 = data[6]
            Sg_FI_e1 = data[7]
            Oil_p_FI_e1 = data[8]
            Gas_p_FI_e1 = data[9]
            pressure_FI_e1 = data[4]/1e3
            time_FI_e1 = data[3]
            zC1_FI_e1 = data[10][0]
            zC3_FI_e1 = data[10][1]
            zC6_FI_e1 = data[10][2]
            zC10_FI_e1 = data[10][3]
            zC15_FI_e1 = data[10][4]
            zC20_FI_e1 = data[10][5]

        datas5 = np.load('flying/Teste_criterio_de_parada/32_CV/.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas5[1:]:
            Sw_FI_e0 = data[5]
            So_FI_e0 = data[6]
            Sg_FI_e0 = data[7]
            Oil_p_FI_e0 = data[8]
            Gas_p_FI_e0 = data[9]
            pressure_FI_e0 = data[4]/1e3
            time_FI_e0 = data[3]
            zC1_FI_e0 = data[10][0]
            zC3_FI_e0 = data[10][1]
            zC6_FI_e0 = data[10][2]
            zC10_FI_e0 = data[10][3]
            zC15_FI_e0 = data[10][4]
            zC20_FI_e0 = data[10][5]"""

        
        plt.figure(1)
        plt.title('t = 20 days - 32x1x1 mesh')
        plt.plot(x1, pressure_FI_e6, '+y')
        plt.plot(x1, pressure_FI_e3, '-r')
        plt.plot(x1, pressure_FI_e2, ':b')
        #plt.plot(x1, pressure_FI_e1, '+g')
        #plt.plot(x1, pressure_FI_e0, '*c')
        plt.ylabel('Pressure (kPa)')
        plt.xlabel('Distance')
        plt.legend(('1e-6', '1e-3', '1e-2'))
        plt.grid()
        plt.savefig('results/Teste_criterio_de_parada/6k/32_CV/pressure_20days' + '.png')

        plt.figure(2)
        plt.title('t = 20 days - 32x1x1 mesh')
        plt.plot(x1, So_FI_e6, '+y')
        plt.plot(x1, So_FI_e3, '-r')
        plt.plot(x1, So_FI_e2, ':b')
        #plt.plot(x1, So_FI_e1, '+g')
        #plt.plot(x1, So_FI_e0, '*c')
        plt.plot(x_CMG, So_CMG, 'k')
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance')
        plt.legend(('1e-6', '1e-3', '1e-2', 'CMG'))
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/Teste_criterio_de_parada/6k/32_CV/saturation_oil_6k_20days' + '.png')

        plt.figure(3)
        plt.title('t = 20 days - 32x1x1 mesh')
        plt.plot(x1, Sw_FI_e6, '+y')
        plt.plot(x1, Sw_FI_e3, '-r')
        plt.plot(x1, Sw_FI_e2, ':b')
        #plt.plot(x1, Sw_FI_e1, '+g')
        #plt.plot(x1, Sw_FI_e0, '*c')
        plt.plot(x_CMG, Sw_CMG, 'k')
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.legend(('1e-6', '1e-3', '1e-2', 'CMG'))
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/Teste_criterio_de_parada/6k/32_CV/saturation_water_6k_20days' + '.png')

        plt.figure(4)
        plt.title('t = 20 days - 32x1x1 mesh')
        plt.plot(x1, Sg_FI_e6, '+y')
        plt.plot(x1, Sg_FI_e3, '-r')
        plt.plot(x1, Sg_FI_e2, ':b')
        #plt.plot(x1, Sg_FI_e1, '+g')
        #plt.plot(x1, Sg_FI_e0, '*c')
        plt.plot(x_CMG, Sg_CMG, 'k')
        plt.legend(('1e-6', '1e-3', '1e-2', 'CMG'))
        plt.ylabel('Gas saturation')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/Teste_criterio_de_parada/6k/32_CV/saturation_gas_6k_20days' + '.png')


        print('Done')
        import pdb; pdb.set_trace()
