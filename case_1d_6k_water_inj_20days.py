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

        datas = np.load('flying/6k/results_6k_SchmallNEW_32_IMPEC_FOU_20days_Stones_1002.npy', allow_pickle=True)
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

        """datas2 = np.load('flying/6k/.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas2[1:]:
            Sw_MUSCL = data[5]
            So_MUSCL = data[6]
            Sg_MUSCL = data[7]
            Oil_p_MUSCL = data[8]
            Gas_p_MUSCL = data[9]
            pressure_MUSCL = data[4]/1e3
            time_MUSCL = data[3]
            zC1_MUSCL = data[10][0]
            zC3_MUSCL = data[10][1]
            zC6_MUSCL = data[10][2]
            zC10_MUSCL = data[10][3]
            zC15_MUSCL = data[10][4]
            zC20_MUSCL = data[10][5]"""

        datas3 = np.load('flying/6k/results_6k_SchmallNEW_32_FI_20days_KRNOVO_1001.npy', allow_pickle=True)
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


        #import pdb; pdb.set_trace()
        plt.figure(1)
        plt.title('t = 20 days - 32x1x1 mesh')
        plt.plot(x1, pressure_FOU, '-k+')
        #plt.plot(x1, pressure_MUSCL, 'g')
        plt.plot(x1, pressure_FI, '-bo')
        #plt.plot(x1, pressure_FI2, 'r')
        plt.ylabel('Pressure (kPa)')
        plt.xlabel('Distance')
        plt.legend(('IMPEC', 'FI'))
        #plt.legend(('IMPEC', 'FI', 'FI Teste'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.grid()
        plt.savefig('results/6k_32/20days/pressure_6k_20days_2' + '.png')

        plt.figure(2)
        plt.title('t = 20 days - 32x1x1 mesh')
        plt.plot(x1, So_FOU, '-k+')
        #plt.plot(x1, So_MUSCL, 'g')
        plt.plot(x1, So_FI, '-bo')
        plt.plot(x_CMG, So_CMG, 'r')
        #plt.legend(('IMPEC', 'FI'))
        plt.legend(('IMPEC', 'FI', 'CMG'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_32/20days/saturation_oil_6k_20days_2' + '.png')

        plt.figure(3)
        plt.title('t = 20 days - 32x1x1 mesh')
        plt.plot(x1, Sw_FOU, '-k+')
        #plt.plot(x1, Sw_MUSCL, 'g')
        plt.plot(x1, Sw_FI, '-bo')
        plt.plot(x_CMG, Sw_CMG, 'r')
        plt.legend(('IMPEC', 'FI'))
        plt.legend(('IMPEC', 'FI', 'CMG'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_32/20days/saturation_water_6k_20days_2' + '.png')

        plt.figure(4)
        plt.title('t = 20 days - 32x1x1 mesh')
        plt.plot(x1, Sg_FOU, '-k+')
        #plt.plot(x1, Sg_MUSCL, 'g')
        plt.plot(x1, Sg_FI, '-bo')
        plt.plot(x_CMG, Sg_CMG, 'r')
        #plt.legend(('IMPEC', 'FI'))
        plt.legend(('IMPEC', 'FI', 'CMG'))
        #plt.legend(('IMPEC', 'MUSCL', 'FI dt 864', 'FI dt 1728'))
        plt.ylabel('Gas saturation')
        plt.xlabel('Distance')
        plt.grid()
        #plt.ylim((0,1))
        plt.savefig('results/6k_32/20days/saturation_gas_6k_20days_2' + '.png')

        """plt.figure(5)
        plt.title('t = 20 days - 32x1x1 mesh')
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
        plt.savefig('results/6k_32/20days/zC1_6k_FOU_20days' + '.png')

        plt.figure(6)
        plt.title('t = 20 days - 32x1x1 mesh')
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
        plt.savefig('results/6k_32/20days/zC3_6k_FOU_20days' + '.png')

        plt.figure(7)
        plt.title('t = 20 days - 32x1x1 mesh')
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
        plt.savefig('results/6k_32/20days/zC6_6k_FOU_20days' + '.png')

        plt.figure(8)
        plt.title('t = 20 days - 32x1x1 mesh')
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
        plt.savefig('results/6k_32/20days/zC10_6k_FOU_20days' + '.png')

        plt.figure(9)
        plt.title('t = 20 days - 32x1x1 mesh')
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
        plt.savefig('results/6k_32/20days/zC15_6k_FOU_20days' + '.png')

        plt.figure(10)
        plt.title('t = 20 days - 32x1x1 mesh')
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
        plt.savefig('results/6k_32/20days/zC20_6k_FOU_20days' + '.png')"""

        print('Done')
        import pdb; pdb.set_trace()
