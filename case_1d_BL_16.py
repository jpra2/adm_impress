import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
xD = np.loadtxt('case_plots/x_BL_semi_analytical.txt')
SwD = np.loadtxt('case_plots/Sw_BL_semi_analytical.txt')

for arq in arquivos:
    if  arq.startswith(name):
        n=16

        """datas = np.load('flying/BL_Teste2/results_Buckley_Leverett_case_32_IMPEC_TESTE2_60.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw = data[5]
            So = data[6]
            Sg = data[7]
            Oil_p = data[8]
            Gas_p = data[9]
            pressure = data[4]/1e3
            time = data[3]
            x1 = np.linspace(0.0, 0.6096, n)
            x1 = x1 / 0.6096
            #import pdb; pdb.set_trace()"""

        datas2 = np.load('flying/BL_Teste2/results_Buckley_Leverett_case_16_FI_36.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas2[1:]:
            Sw_FI = data[5]
            So_FI = data[6]
            Sg_FI = data[7]
            Oil_p_FI = data[8]
            Gas_p_FI = data[9]
            pressure_FI = data[4]/1e3
            time_FI = data[3]
            x1 = np.linspace(0.0, 0.6096, n)
            x1 = x1 / 0.6096

        datas3 = np.load('flying/results_Buckley_Leverett_case_16_FI_NEW_36.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI_new = data[5]
            So_FI_new = data[6]
            Sg_FI_new = data[7]
            Oil_p_FI_new = data[8]
            Gas_p_FI_new = data[9]
            pressure_FI_new = data[4]/1e3
            time_FI_new = data[3]

        """plt.figure(1)
        plt.title('BL Pressure - 32x1x1 mesh')
        plt.plot(x1, pressure, 'k')
        plt.plot(x1, pressure_FI, '-bo')
        plt.plot(x1, pressure_FI_new, '-r+')
        plt.ylabel('Pressure (kPa)')
        plt.xlabel('Distance')
        plt.legend(('IMPEC', 'Fully Implicit - back', 'Fully Implicit - new'))
        plt.grid()
        plt.savefig('results/BL_Pressure_32_1_2' + '{}'.format(n) + '.png')"""

        plt.figure(2)
        #plt.title('BL Sw - 32x1x1 mesh')
        #plt.plot(x1, Sw, 'k')
        #plt.plot(x1, Sw_FI, '-bo')
        plt.plot(xD, SwD, 'b')
        plt.plot(x1, Sw_FI_new, 'r')
        plt.legend(('Analytical Solution', 'Fully Implicit'))
        #plt.legend(('IMPEC', 'Fully Implicit - back', 'Analytical Solution', 'Fully Implicit - new'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance (m)')
        plt.grid()
        plt.savefig('results/BL_Sw_16_' + '{}'.format(n) + '.png')

        """plt.figure(3)
        plt.title('BL So - 32x1x1 mesh')
        plt.plot(x1, So, 'k')
        plt.plot(x1, So_FI, '-bo')
        plt.plot(x1, So_FI_new, '-r+')
        plt.legend(('IMPEC', 'Fully Implicit - back', 'Fully Implicit - new'))
        plt.ylabel('Oil saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/BL_So_32_1_2' + '{}'.format(n) + '.png')

        plt.figure(4)
        plt.title('BL Sg - 32x1x1 mesh')
        plt.plot(x1, Sg, 'k')
        plt.plot(x1, Sg_FI, '-bo')
        plt.plot(x1, Sg_FI_new, '-r+')
        plt.legend(('IMPEC', 'Fully Implicit - back', 'Fully Implicit - new'))
        plt.ylabel('Saturação de gás ')
        plt.xlabel('Distância')
        plt.grid()
        plt.savefig('results/BL_Sg_32_1_2' +'{}'.format(n) + '.png')"""

        import pdb; pdb.set_trace()