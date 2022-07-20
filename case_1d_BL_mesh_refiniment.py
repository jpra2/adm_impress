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

        datas = np.load('flying/BL_Teste2/results_Buckley_Leverett_case_128_FI_TESTE_139.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FI_128 = data[5]
            So_FI_128 = data[6]
            Sg_FI_128 = data[7]
            Oil_p_FI_128 = data[8]
            Gas_p_FI_128 = data[9]
            pressure_FI_128 = data[4]/1e3
            time_FI_128 = data[3]
            x128 = np.linspace(0.0, 0.6096, 128)
            x128 = x128/0.6096

        datas2 = np.load('flying/BL_Teste2/results_Buckley_Leverett_case_64_FI_TESTE_70.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas2[1:]:
            Sw_FI_64 = data[5]
            So_FI_64 = data[6]
            Sg_FI_64 = data[7]
            Oil_p_FI_64 = data[8]
            Gas_p_FI_64 = data[9]
            pressure_FI_64 = data[4]/1e3
            time_FI_64 = data[3]
            x64 = np.linspace(0.0, 0.6096, 64)
            x64 = x64/0.6096

        datas3 = np.load('flying/BL_Teste2/results_Buckley_Leverett_32_FI_BACKTOBACK_43.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI_32 = data[5]
            So_FI_32 = data[6]
            Sg_FI_32 = data[7]
            Oil_p_FI_32 = data[8]
            Gas_p_FI_32 = data[9]
            pressure_FI_32 = data[4]/1e3
            time_FI_32 = data[3]
            x32 = np.linspace(0.0, 0.6096, 32)
            x32 = x32/0.6096


        plt.figure(1)
        plt.title('BL Sw - mesh refinement')
        plt.plot(xD, SwD, 'b')
        plt.plot(x32, Sw_FI_32, ':g')
        plt.plot(x64, Sw_FI_64, '.-k')
        plt.plot(x128, Sw_FI_128, '--r')
        plt.legend(('Analytical Solution', '32 CV', '64 CV', '128 CV'))
        #plt.legend(('IMPEC', 'Fully Implicit - back', 'Analytical Solution', 'Fully Implicit - new'))
        plt.ylabel('Water saturation')
        plt.xlabel('Distance')
        plt.grid()
        plt.savefig('results/BL_Sw_refinement' + '.png')

        import pdb; pdb.set_trace()
