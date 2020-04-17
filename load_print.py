import numpy as np
import matplotlib.pyplot as plt
import os

flying = 'flying'
name = 'all_compositional_'
arquivos = os.listdir(flying)
i=1
for arq in arquivos:
    if  arq.startswith(name):
        datas = np.load('flying/all_compositional_mono_comp_results_22.npy', allow_pickle=True)

        for data in datas[3:]:
            #pressure = data[4] / 6894.76
            pressure = (data[4] - 13.78951458E6 * np.ones(100))/6894.76
            time = data[3]
            loop = data[0]
            flux = data[5]
            flux_vector = data[6]
            i=loop
            #flux_vols = data[5]
            x = np.linspace(0,1,100)
            plt.figure(1)
            plt.plot(x, pressure)
            plt.grid()
            plt.ylabel('Pressure (psi)')
            plt.xlabel('Dimensionless distance')
            plt.savefig('results/compositional/pressure_mc' + str(loop) + '.png')

            import pdb; pdb.set_trace()


'''
        p1 = (datas[1][4] - 13.78951458E6 * np.ones(100))/6894.76
        p2 = (datas[2][4] - 13.78951458E6 * np.ones(100))/6894.76
        p3 = (datas[3][4] - 13.78951458E6 * np.ones(100))/6894.76
        t1 = datas[1][3]
        t2 = datas[2][3]
        t3 = datas[3][3]
        x = np.linspace(0,1,100)
        plt.figure(1)
        plt.plot(x, p1, 'r', x, p2, 'b', x, p3, 'g')
        plt.grid()
        plt.legend(('%f seg' %t1, '%f seg' %t2, '%f seg' %t3, ))
        plt.ylabel('Pressure Drop (psi)')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/pressure_mc' + '.png')
        import pdb; pdb.set_trace()'''
