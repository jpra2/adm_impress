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
            x1_50 = np.linspace(0, 50, n_50)

        datas2 = np.load('flying/6k/100_0_0/results_6k_SchmallNEW_100_FI_20days_POCONEW_dtmenor_2003.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas2[1:]:
            Sw_FI_2 = data[5]
            So_FI_2 = data[6]
            Sg_FI_2 = data[7]
            Oil_p_FI_2 = data[8]
            Gas_p_FI_2 = data[9]
            pressure_FI_2 = data[4]/1e3
            zC1_FI_2 = data[10][0]
            zC3_FI_2 = data[10][1]
            zC6_FI_2 = data[10][2]
            zC10_FI_2 = data[10][3]
            zC15_FI_2 = data[10][4]
            zC20_FI_2 = data[10][5]
            x1_100 = np.linspace(0, 50, n_100)

        datas3 = np.load('flying/6k/200_0_0/results_6k_SchmallNEW_200_FI_20days_POCONEW_2003.npy', allow_pickle=True)
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
            x1_200 = np.linspace(0, 50, n_200)


        
        plt.figure(1)
        plt.title('t = 20 days - 200x1x1 mesh')
        plt.plot(x1_50, pressure_FOU, 'k')
        #plt.plot(x1, pressure_MUSCL, 'g')
        #plt.plot(x_CMG, pressure_CMG, '-bo')
        plt.plot(x1_200, pressure_FI, 'r')
        plt.plot(x1_100, pressure_FI_2, '-b')
        plt.ylabel('Pressure (kPa)')
        plt.xlabel('Distance')
        plt.legend(('IMPEC 50 CVs', 'FI novo poço 200 CVs', 'FI novo poço 100 CVs dt menor'))
        #plt.legend(('IMPEC', 'FI', 'FI Teste'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.grid()
        plt.savefig('results/6k_200/20days/pressure_6k_20days' + '.png')

        plt.figure(2)
        plt.title('t = 20 days - 200x1x1 mesh')
        #plt.plot(x1, So_FOU, 'k')
        #plt.plot(x1, So_MUSCL, 'g')
        plt.plot(x_CMG, So_CMG, 'k')
        plt.plot(x1_200, So_FI, 'r')
        #plt.plot(x1_100, So_FI_2, '-b')
        #plt.legend(('CMG', 'FI novo poço 200 CVs', 'FI novo poço 100 CVs dt menor'))
        plt.legend(('CMG', 'FI novo poço 200 CVs'))
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_200/20days/saturation_oil_6k_20days' + '.png')

        plt.figure(3)
        plt.title('t = 20 days - 200x1x1 mesh')
        #plt.plot(x1, Sw_FOU, 'k')
        #plt.plot(x1, Sw_MUSCL, 'g')
        plt.plot(x_CMG, Sw_CMG, 'k')
        plt.plot(x1_200, Sw_FI, 'r')
        #plt.plot(x1_100, Sw_FI_2, '-b')
        plt.legend(('CMG', 'FI novo poço 200 CVs', 'FI novo poço 100 CVs dt menor'))
        plt.legend(('CMG', 'FI novo poço 200 CVs'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_200/20days/saturation_water_6k_20days' + '.png')

        plt.figure(4)
        plt.title('t = 20 days - 200x1x1 mesh')
        #plt.plot(x1, Sg_FOU, 'k')
        #plt.plot(x1, Sg_MUSCL, 'g')
        plt.plot(x_CMG, Sg_CMG, 'k')
        plt.plot(x1_200, Sg_FI, 'r')
        #plt.plot(x1_100, Sg_FI_2, '-b')
        plt.legend(('CMG', 'FI novo poço 200 CVs', 'FI novo poço 100 CVs dt menor'))
        plt.legend(('CMG', 'FI novo poço 200 CVs'))
        plt.ylabel('Gas saturation')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_200/20days/saturation_gas_6k_20days' + '.png')

        """plt.figure(5)
        plt.title('t = 10 days - 50x1x1 mesh')
        plt.plot(x1, zC1_FOU, 'k')
        #plt.plot(x1, zC1_MUSCL, 'g')
        plt.plot(x1, zC1_FI, '-bo')
        #plt.plot(x1, zC1_FI2, 'r')
        plt.legend(('IMPEC', 'FI'))
        #plt.legend(('IMPEC', 'FI', 'FI Teste'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.ylabel('zC1')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_50/10days/zC1_6k_FOU_10days' + '.png')

        plt.figure(6)
        plt.title('t = 10 days - 50x1x1 mesh')
        plt.plot(x1, zC3_FOU, 'k')
        #plt.plot(x1, zC3_MUSCL, 'g')
        plt.plot(x1, zC3_FI, '-bo')
        #plt.plot(x1, zC3_FI2, 'r')
        plt.legend(('IMPEC', 'FI'))
        #plt.legend(('IMPEC', 'FI', 'FI Teste'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.ylabel('zC3')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_50/10days/zC3_6k_FOU_10days' + '.png')

        plt.figure(7)
        plt.title('t = 10 days - 50x1x1 mesh')
        plt.plot(x1, zC6_FOU, 'k')
        #plt.plot(x1, zC6_MUSCL, 'g')
        plt.plot(x1, zC6_FI, '-bo')
        #plt.plot(x1, zC6_FI2, 'r')
        plt.legend(('IMPEC', 'FI'))
        #plt.legend(('IMPEC', 'FI', 'FI Teste'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.ylabel('zC6')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_50/10days/zC6_6k_FOU_10days' + '.png')

        plt.figure(8)
        plt.title('t = 10 days - 50x1x1 mesh')
        plt.plot(x1, zC10_FOU, 'k')
        #plt.plot(x1, zC10_MUSCL, 'g')
        plt.plot(x1, zC10_FI, '-bo')
        #plt.plot(x1, zC10_FI2, 'r')
        plt.legend(('IMPEC', 'FI'))
        #plt.legend(('IMPEC', 'FI', 'FI Teste'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.ylabel('zC10')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_50/10days/zC10_6k_FOU_10days' + '.png')

        plt.figure(9)
        plt.title('t = 10 days - 50x1x1 mesh')
        plt.plot(x1, zC15_FOU, 'k')
        #plt.plot(x1, zC15_MUSCL, 'g')
        plt.plot(x1, zC15_FI, '-bo')
        #plt.plot(x1, zC15_FI2, 'r')
        plt.legend(('IMPEC', 'FI'))
        #plt.legend(('IMPEC', 'FI', 'FI Teste'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.ylabel('zC15')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_50/10days/zC15_6k_FOU_10days' + '.png')

        plt.figure(10)
        plt.title('t = 10 days - 50x1x1 mesh')
        plt.plot(x1, zC20_FOU, 'k')
        #plt.plot(x1, zC20_MUSCL, 'g')
        plt.plot(x1, zC20_FI, '-bo')
        #plt.plot(x1, zC20_FI2, 'r')
        plt.legend(('IMPEC', 'FI'))
        #plt.legend(('IMPEC', 'FI', 'FI Teste'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.ylabel('zC20')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_50/10days/zC20_6k_FOU_10days' + '.png')"""

        print('Done')
        import pdb; pdb.set_trace()
