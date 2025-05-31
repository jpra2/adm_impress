import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import InterpolatedUnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)

for  arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_Buckley_Leverett_case_80x80_IMPEC_MUSCL_6328.npy', allow_pickle=True)
        for data in datas[-1:]:
            Sw_80x80_MUSCL = data[5]
            centroid = data[13]
            centroid_x = centroid[:,0]
            centroid_y = centroid[:,1]
            ind = abs(centroid_x - (0.6096-centroid_y))<1e-15
            Sw_diagonal_80x80_MUSCL = Sw_80x80_MUSCL[ind]
            x_diagonal_80x80_MUSCL = centroid_x[ind]
        
        datas = np.load('flying/results_Buckley_Leverett_case_80x80_IMPEC_FOU_2583.npy', allow_pickle=True)
        datas = np.load('flying/results_Buckley_Leverett_case_80x80_2D_IMPEC_FOU_2583.npy', allow_pickle=True)
        for data in datas[-1:]:
            Sw_80x80_FOU = data[5]
            centroid = data[13]
            centroid_x = centroid[:,0]
            centroid_y = centroid[:,1]
            ind = abs(centroid_x - (0.6096-centroid_y))<1e-15
            Sw_diagonal_80x80_FOU = Sw_80x80_FOU[ind]
            x_diagonal_80x80_FOU = centroid_x[ind]

        datas = np.load('results/flying/results_Buckley_Leverett_case_2D_80x80_IMPEC_663.npy', allow_pickle=True)
        for data in datas[-1:]:
            Sw_80x80_FOU_ref = data[5]
            centroid = data[11]
            centroid_x = centroid[:,0]
            centroid_y = centroid[:,1]
            ind = abs(centroid_x - (0.6096-centroid_y))<1e-15
            Sw_diagonal_80x80_FOU_ref = Sw_80x80_FOU_ref[ind]
            x_diagonal_80x80_FOU_ref = centroid_x[ind]

        
        plt.figure(1)
        plt.plot(x_diagonal_80x80_FOU_ref, Sw_diagonal_80x80_FOU_ref, 'm')
        plt.plot(x_diagonal_80x80_FOU, Sw_diagonal_80x80_FOU, 'b')
        plt.plot(x_diagonal_80x80_MUSCL, Sw_diagonal_80x80_MUSCL, 'g')
        plt.grid()
        plt.legend(('FOU 80x80 ref', 'FOU 80x80', 'MUSCL 80x80'))
        plt.title('Water saturation profile')
        plt.ylabel('Water saturation')
        plt.xlabel('Diagonal distance')
        plt.savefig('results/compositional/Swdiag_buckley2D_80x80.png')

        import pdb; pdb.set_trace()
