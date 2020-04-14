import numpy as np
import matplotlib.pyplot as plt
import os

flying = 'flying'
name = 'all_compositional_'
arquivos = os.listdir(flying)
i=1
for arq in arquivos:
    if  arq.startswith(name):
        datas = np.load('flying/all_compositional_m2_results_21.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure = (data[4] - 13.78951458E6 * np.ones(100))/6894.76
            time = data[3]
            loop = data[0]
            i=loop
            #flux_vols = data[5]
            x = np.linspace(0,1,100)
            plt.figure(i)
            plt.plot(x, pressure)
            plt.grid()
            plt.ylabel('Pressure Drop')
            plt.xlabel('Dimensionless distance')
            plt.savefig('results/compositional/results_' + str(loop) + '.png')
