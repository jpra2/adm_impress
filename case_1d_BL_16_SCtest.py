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

        datas = np.load('flying/Teste_criterio_de_parada/BL/16_CV/results_Buckley_Leverett_case_16_FI_ADIMENSIONALIZADO_e9dt20_36.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FI_e9 = data[5]
            So_FI_e9 = data[6]
            Sg_FI_e9 = data[7]
            Oil_p_FI_e9 = data[8]
            Gas_p_FI_e9 = data[9]
            pressure_FI_e9 = data[4]/1e3
            time_FI_e9 = data[3]
            x1 = np.linspace(0.0, 0.6096, n)
            x1 = x1 / 0.6096
            #import pdb; pdb.set_trace()

        datas = np.load('flying/Teste_criterio_de_parada/BL/16_CV/results_Buckley_Leverett_case_16_FI_ADIMENSIONALIZADO_e6dt20_36.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FI_e6 = data[5]
            So_FI_e6 = data[6]
            Sg_FI_e6 = data[7]
            Oil_p_FI_e6 = data[8]
            Gas_p_FI_e6 = data[9]
            pressure_FI_e6 = data[4]/1e3
            time_FI_e6 = data[3]

        datas = np.load('flying/Teste_criterio_de_parada/BL/16_CV/results_Buckley_Leverett_case_16_FI_ADIMENSIONALIZADO_e3dt20_36.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FI_e3 = data[5]
            So_FI_e3 = data[6]
            Sg_FI_e3 = data[7]
            Oil_p_FI_e3 = data[8]
            Gas_p_FI_e3 = data[9]
            pressure_FI_e3 = data[4]/1e3
            time_FI_e3 = data[3]

        datas = np.load('flying/Teste_criterio_de_parada/BL/16_CV/results_Buckley_Leverett_case_16_FI_ADIMENSIONALIZADO_e2dt20_36.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FI_e2 = data[5]
            So_FI_e2 = data[6]
            Sg_FI_e2 = data[7]
            Oil_p_FI_e2 = data[8]
            Gas_p_FI_e2 = data[9]
            pressure_FI_e2 = data[4]/1e3
            time_FI_e2 = data[3]


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
        plt.title('BL Sw - 8x1x1 mesh')
        #plt.plot(x1, Sw, 'k')
        #plt.plot(x1, Sw_FI, '-bo')
        plt.plot(xD, SwD, 'k')
        plt.plot(x1, Sw_FI_e9, 'r')
        plt.plot(x1, Sw_FI_e6, '-bo')
        plt.plot(x1, Sw_FI_e3, '+y')
        plt.plot(x1, Sw_FI_e2, ':g')
        #plt.legend(('Analytical Solution', '1e-9', '1e-9 dt maior'))
        plt.legend(('Analytical Solution', '1e-9', '1e-6', '1e-3', '1e-2'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance (m)')
        plt.grid()
        plt.savefig('results/Teste_criterio_de_parada/BL/16_CV/BL_Sw_16_' + '{}'.format(n) + '.png')

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
