import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import InterpolatedUnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)

for  arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_Buckley_Leverett_case_2D_20x20_FOU_403.npy', allow_pickle=True)
        for data in datas[-1:]:
            Sw_20x20_FOU = data[5]
            centroid = data[11]
            centroid_x = centroid[:,0]
            centroid_y = centroid[:,1]
            ind = abs(centroid_x - (0.6096-centroid_y))<1e-15
            Sw_diagonal_20x20_FOU = Sw_20x20_FOU[ind]
            x_diagonal_20x20_FOU = centroid_x[ind]

        datas = np.load('flying/results_Buckley_Leverett_case_2D_40x40_FOU_1586.npy', allow_pickle=True)
        for data in datas[-1:]:
            Sw_40x40_FOU = data[5]
            centroid = data[11]
            centroid_x = centroid[:,0]
            centroid_y = centroid[:,1]
            ind = abs(centroid_x - (0.6096-centroid_y))<1e-15
            Sw_diagonal_40x40_FOU = Sw_40x40_FOU[ind]
            x_diagonal_40x40_FOU = centroid_x[ind]

        datas = np.load('flying/results_Buckley_Leverett_case_2D_20x20_MUSCL_407.npy', allow_pickle=True)
        for data in datas[-1:]:
            Sw_20x20_MUSCL = data[5]
            centroid = data[11]
            centroid_x = centroid[:,0]
            centroid_y = centroid[:,1]
            ind = abs(centroid_x - (0.6096-centroid_y))<1e-15
            Sw_diagonal_20x20_MUSCL = Sw_20x20_MUSCL[ind]
            x_diagonal_20x20_MUSCL = centroid_x[ind]

        datas = np.load('flying/results_Buckley_Leverett_case_2D_40x40_MUSCL_1586.npy', allow_pickle=True)
        for data in datas[-1:]:
            Sw_40x40_MUSCL = data[5]
            centroid = data[11]
            centroid_x = centroid[:,0]
            centroid_y = centroid[:,1]
            ind = abs(centroid_x - (0.6096-centroid_y))<1e-15
            Sw_diagonal_40x40_MUSCL = Sw_40x40_MUSCL[ind]
            x_diagonal_40x40_MUSCL = centroid_x[ind]
            import pdb; pdb.set_trace()


        plt.figure(1)
        plt.plot(x_diagonal_20x20_FOU, Sw_diagonal_20x20_FOU, 'm')
        plt.plot(x_diagonal_20x20_MUSCL, Sw_diagonal_20x20_MUSCL, 'b')
        plt.grid()
        plt.legend(('FOU 20x20', 'MUSCL 20x20'))
        plt.title('Water saturation profile')
        plt.ylabel('Water saturation')
        plt.xlabel('Diagonal distance')
        plt.savefig('results/compositional/Swdiag_buckley2D_20x20.png')

        plt.figure(2)
        plt.plot(x_diagonal_40x40_FOU, Sw_diagonal_40x40_FOU, 'm')
        plt.plot(x_diagonal_40x40_MUSCL, Sw_diagonal_40x40_MUSCL, 'b')
        plt.grid()
        plt.legend(('FOU 40x40', 'MUSCL 40x40'))
        plt.title('Water saturation profile')
        plt.ylabel('Water saturation')
        plt.xlabel('Diagonal distance')
        plt.savefig('results/compositional/Swdiag_buckley2D_40x40.png')

        import pdb; pdb.set_trace()
