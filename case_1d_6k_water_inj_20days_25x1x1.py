import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
n_200 = 200
n_100 = 100
n_50 = 50
n_25 = 25

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


        datas3 = np.load('flying/6k/25_0_0/results_6k_SchmallNEW_25_FI_20days_POCONEW_40.npy', allow_pickle=True)
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
            x1_25 = np.linspace(0, 50, n_25)


        plt.figure(2)
        plt.title('t = 20 days - 25x1x1 mesh')
        #plt.plot(x1, So_FOU, 'k')
        #plt.plot(x1, So_MUSCL, 'g')
        plt.plot(x_CMG, So_CMG, 'k')
        plt.plot(x1_25, So_FI, '--r')
        #plt.plot(x1_100, So_FI_2, '-b')
        #plt.legend(('CMG', 'FI novo poço 200 CVs', 'FI novo poço 100 CVs dt menor'))
        plt.legend(('CMG - 5000 elements', 'FI'))
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance (m)')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k/6k_25/20days/saturation_oil_6k_20days' + '.png')

        plt.figure(3)
        plt.title('t = 20 days - 25x1x1 mesh')
        #plt.plot(x1, Sw_FOU, 'k')
        #plt.plot(x1, Sw_MUSCL, 'g')
        plt.plot(x_CMG, Sw_CMG, 'k')
        plt.plot(x1_25, Sw_FI, '--b')
        #plt.plot(x1_100, Sw_FI_2, '-b')
        plt.legend(('CMG - 5000 elements', 'FI'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance (m)')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k/6k_25/20days/saturation_water_6k_20days' + '.png')

        plt.figure(4)
        plt.title('t = 20 days - 25x1x1 mesh')
        #plt.plot(x1, Sg_FOU, 'k')
        #plt.plot(x1, Sg_MUSCL, 'g')
        plt.plot(x_CMG, Sg_CMG, 'k')
        plt.plot(x1_25, Sg_FI, '--g')
        #plt.plot(x1_100, Sg_FI_2, '-b')
        plt.legend(('CMG - 5000 elements', 'FI'))
        plt.ylabel('Gas saturation')
        plt.xlabel('Distance (m)')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k/6k_25/20days/saturation_gas_6k_20days' + '.png')



        print('Done')
        import pdb; pdb.set_trace()
