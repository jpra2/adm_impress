import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)

for  arq in arquivos:
    if  arq.startswith(name):
        '''-------------------------MUSCL LLF RESULTS------------------------'''

        zCO2_x = np.array([0,0,0.0251938, 0.0377907, 0.0494186, 0.0818798, 0.120155, \
                            0.150194, 0.202035, 0.248062, 0.300872, 0.374031, \
                            0.436047, 0.506783, 0.560087, 0.619186, 0.682655, \
                            0.750484, 0.82655, 0.895833, 0.95688, 1.02374, \
                            1.0843, 1.16037, 1.23062, 1.297, 1.297, 1.5])
        zCO2 = np.array([0.995, 0.78667, 0.741667, 0.707222, 0.686111, 0.639444, \
                        0.604444, 0.583333, 0.558333, 0.541667, 0.526111, 0.508333, \
                        0.496667, 0.486111, 0.480556, 0.472222, 0.466667, 0.459444, \
                        0.452778, 0.447778, 0.443333, 0.438889, 0.436667, 0.432222,  \
                        0.427778, 0.426111, 0., 0.])


        datas = np.load('flying/results_Orr_2k_C4_IMPEC_FOU_525.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            zCO2_200_FOU = data[10][0]
            n=200
            x_200 = np.linspace(0+1.5/(2*n),1.5*(1-1/(2*n)),n)


        plt.figure(1)
        plt.plot(x_200, zCO2_200_FOU, '-bo', mfc='none')
        plt.plot(zCO2_x, zCO2, '-k')
        #plt.plot(x_axis_xCH4_LLF_50, xCH4_LLF_50, '-sk', mfc='none')
        plt.grid()
        plt.legend('FOU - 200', 'MOC')
        plt.ylabel('CO$_2$ global molar fraction ')
        plt.title('HIgh volatile intermediate component case - Orr')
        plt.xlabel('Distance')
        plt.savefig('results/compositional/2k_CO2_Orr_NVCM.png')
        import pdb; pdb.set_trace()
