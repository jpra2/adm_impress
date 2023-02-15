import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
n= 400

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

fzC1 = open('ZC1_6k.txt','r')
zC1_CMG = [float(line.rstrip('\n\r')) for line in fzC1]

fzC20 = open('ZC20_6k.txt','r')
zC20_CMG = [float(line.rstrip('\n\r')) for line in fzC20]


for arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_6k_SchmallNEW_400_IMPEC_20days_1568.npy', allow_pickle=True)
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

        datas = np.load('flying/results_6k_SchmallNEW_400_FI_20days_CFL3_0_558.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FI = data[5]
            So_FI = data[6]
            Sg_FI = data[7]
            Oil_p_FI = data[8]
            Gas_p_FI = data[9]
            pressure_FI = data[4]/1e3
            time_FI = data[3]
            zC1_FI = data[10][0]
            zC3_FI = data[10][1]
            zC6_FI = data[10][2]
            zC10_FI = data[10][3]
            zC15_FI = data[10][4]
            zC20_FI = data[10][5]



        #import pdb; pdb.set_trace()
        plt.figure(1)
        #plt.title('t = 20 days - 200x1x1 mesh')
        plt.plot(x_CMG, So_CMG, 'k')
        #plt.plot(x1_200, So_FI, '--r')
        #plt.plot(x1, So_FOU, 'r')
        plt.plot(x1, So_FI, '--r')
        plt.legend(('CMG - 5000 CVs', 'FI 400 CVs'))
        plt.ylabel('Saturação de óleo')
        plt.xlabel('Distância (m)')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k/6k_400/saturation_oil_6k_20days_' + '.png')

        plt.figure(2)
        #plt.title('t = 20 days - 200x1x1 mesh')
        plt.plot(x_CMG, Sw_CMG, 'k')
        #plt.plot(x1_200, Sw_FI, '--r')
        #plt.plot(x1, Sw_FOU, 'r')
        plt.plot(x1, Sw_FI, '--b')
        plt.legend(('CMG - 5000 CVs', 'FI 400 CVs'))
        plt.ylabel('Saturação de água')
        plt.xlabel('Distância (m)')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k/6k_400/saturation_water_6k_20days_' + '.png')

        plt.figure(3)
        #plt.title('t = 20 days - 200x1x1 mesh')
        plt.plot(x_CMG, Sg_CMG, 'k')
        #plt.plot(x1_200, Sg_FI, '--r')
        #plt.plot(x1, Sg_FOU, 'r')
        plt.plot(x1, Sg_FI, '--g')
        plt.legend(('CMG - 5000 CVs', 'FI 400 CVs'))
        plt.ylabel('Saturação de gás')
        plt.xlabel('Distância (m)')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k/6k_400/saturation_gas_6k_20days_' + '.png')

        plt.figure(5)
        #plt.title('t = 10 days - 50x1x1 mesh')
        #plt.plot(x1, zC1_FOU, 'k')
        plt.plot(x_CMG, zC1_CMG, 'k')
        plt.plot(x1, zC1_FI, '--b')
        plt.legend(('CMG - 5000 CVs', 'FI 200 CVs'), loc=1)
        plt.ylabel('Fração molar global do C1')
        plt.xlabel('Distância (m)')
        plt.grid()
        plt.ylim((0.37, 0.45))
        plt.savefig('results/6k/6k_400/zC1_6k_FI_20days_lim' + '.png')

        plt.figure(10)
        #plt.title('t = 10 days - 50x1x1 mesh')
        #plt.plot(x1, zC20_FOU, 'k')
        plt.plot(x_CMG, zC20_CMG, 'k')
        plt.plot(x1, zC20_FI, '--b')
        plt.legend(('CMG - 5000 CVs', 'FI 200 CVs'), loc=4)
        #plt.legend(('IMPEC - 5000 CVs', 'FI 200 CVs', 'CMG - 5000 CVs'), loc=1)
        plt.ylabel('Fração molar global do C20')
        plt.xlabel('Distância (m)')
        plt.grid()
        plt.ylim((0.055, 0.064))
        plt.savefig('results/6k/6k_400/zC20_6k_FI_20days_lim' + '.png')

        print('Done')
        import pdb; pdb.set_trace()
