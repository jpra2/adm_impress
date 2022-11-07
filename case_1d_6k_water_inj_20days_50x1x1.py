import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
n = 50

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

        datas = np.load('flying/6k/results_6k_SchmallNEW_50_IMPEC_20days_401.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FOU = data[5]
            So_FOU = data[6]
            Sg_FOU = data[7]
            Oil_p_FOU = data[8]
            Gas_p_FOU = data[9]
            pressure_FOU = data[4]/1e3
            time_FOU = data[3]
            zC1_FOU = data[10][0]
            zC3_FOU = data[10][1]
            zC6_FOU = data[10][2]
            zC10_FOU = data[10][3]
            zC15_FOU = data[10][4]
            zC20_FOU = data[10][5]
            x1 = np.linspace(0, 50, n)


        datas3 = np.load('flying/6k/results_6k_SchmallNEW_50_FI_20days_POCONEW_203.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI = data[5]
            So_FI = data[6]
            Sg_FI = data[7]
            Oil_p_FI = data[8]
            Gas_p_FI = data[9]
            pressure_FI = data[4]/1e3
            zC1_FI = data[10][0]
            zC3_FI = data[10][1]
            zC6_FI = data[10][2]
            zC10_FI = data[10][3]
            zC15_FI = data[10][4]
            zC20_FI = data[10][5]

        datas3 = np.load('flying/6k/results_6k_SchmallNEW_50_FI_20days_2DTeste_203.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI_TesteW = data[5]
            So_FI_TesteW = data[6]
            Sg_FI_TesteW = data[7]
            Oil_p_FI_TesteW = data[8]
            Gas_p_FI_TesteW = data[9]
            pressure_FI_TesteW = data[4]/1e3
            zC1_FI_TesteW = data[10][0]
            zC3_FI_TesteW = data[10][1]
            zC6_FI_TesteW = data[10][2]
            zC10_FI_TesteW = data[10][3]
            zC15_FI_TesteW = data[10][4]
            zC20_FI_TesteW = data[10][5]


        
        plt.figure(1)
        plt.title('t = 20 days - 50x1x1 mesh')
        plt.plot(x1, pressure_FOU, 'k')
        #plt.plot(x1, pressure_MUSCL, 'g')
        #plt.plot(x_CMG, pressure_CMG, '-bo')
        plt.plot(x1, pressure_FI, 'r')
        plt.plot(x1, pressure_FI_TesteW, '--b')
        plt.ylabel('Pressure (kPa)')
        plt.xlabel('Distance')
        plt.legend(('IMPEC', 'FI', 'FI teste final'))
        plt.grid()
        plt.savefig('results/6k/6k_50/20days/pressure_6k_20days' + '.png')

        plt.figure(2)
        plt.title('t = 20 days - 50x1x1 mesh')
        #plt.plot(x1, So_FOU, 'k')
        #plt.plot(x1, So_MUSCL, 'g')
        plt.plot(x_CMG, So_CMG, 'k')
        plt.plot(x1, So_FI, 'r')
        plt.plot(x1, So_FI_TesteW, '--b')
        plt.legend(('CMG', 'FI', 'FI teste final'))
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k/6k_50/20days/saturation_oil_6k_20days' + '.png')

        plt.figure(3)
        plt.title('t = 20 days - 50x1x1 mesh')
        #plt.plot(x1, Sw_FOU, 'k')
        #plt.plot(x1, Sw_MUSCL, 'g')
        plt.plot(x_CMG, Sw_CMG, 'k')
        plt.plot(x1, Sw_FI, 'r')
        plt.plot(x1, Sw_FI_TesteW, '--b')
        plt.legend(('CMG', 'FI', 'FI teste final'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k/6k_50/20days/saturation_water_6k_20days' + '.png')

        plt.figure(4)
        plt.title('t = 20 days - 50x1x1 mesh')
        #plt.plot(x1, Sg_FOU, 'k')
        #plt.plot(x1, Sg_MUSCL, 'g')
        plt.plot(x_CMG, Sg_CMG, 'k')
        plt.plot(x1, Sg_FI, 'r')
        plt.plot(x1, Sg_FI_TesteW, '--b')
        plt.legend(('CMG', 'FI', 'FI teste final'))
        plt.ylabel('Gas saturation')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k/6k_50/20days/saturation_gas_6k_20days' + '.png')

        print('Done')
        import pdb; pdb.set_trace()
