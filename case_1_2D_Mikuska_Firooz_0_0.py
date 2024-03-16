import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import InterpolatedUnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)

for  arq in arquivos:
    if  arq.startswith(name):
        x_diagonal_40x40_CMG = np.array([0, 1.767766953, 3.535533905, 5.303300858, \
        7.07106781, 8.838834763, 10.60660172, 12.37436867, 14.14213562, 15.90990257, \
        17.67766953, 19.44543648, 21.21320343, 22.98097038, 24.74873734, 26.51650429, \
        28.28427124, 30.05203819, 31.81980515, 33.58757401, 35.35533905, 37.12310791, \
        38.89087296, 40.65864182, 42.42640686, 44.19417572, 45.96194077, 47.72970963, \
        49.49747467, 51.26524353, 53.03300858, 54.80077744, 56.56854248, 58.33631134, \
        60.10407639, 61.87184525, 63.63961029, 65.40737915, 67.17514801, 68.94290924])
        zC1_diagonal_40x40_CMG = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0.99999994, \
        1, 1, 1, 1, 1, 1, 1, 1, 0.99999994, 0.999996662, 0.999931812, 0.999212563, \
        0.994407058, 0.973902345, 0.915378213, 0.798873782, 0.629756272, 0.442055434, \
        0.275489837, 0.152939081, 0.076220013, 0.034544773, 0.014491771, 0.005745118, \
        0.002203954, 0.000840748, 0.000329051, 0.000136752, 6.23E-05, 1.28E-05])

        x_diagonal_50x50_CMG = np.array([0, 1.414213538, 2.828427076, 4.242640495, \
        5.656854153, 7.07106781, 8.485280991, 9.899495125, 11.31370831, 12.72792244, \
        14.14213562, 15.5563488, 16.97056198, 18.38477707, 19.79899025, 21.21320343, \
        22.62741661, 24.04162979, 25.45584488, 26.87005806, 28.28427124, 29.69848442, \
        31.1126976, 32.52691269, 33.94112396, 35.35533905, 36.76955414, 38.18376541, \
        39.5979805, 41.01219177, 42.42640686, 43.84062195, 45.25483322, 46.66904831, \
        48.08325958, 49.49747467, 50.91168976, 52.32590103, 53.74011612, 55.15432739, \
        56.56854248, 57.98275757, 59.39696884, 60.81118393, 62.2253952, 63.63961029, \
        65.05382538, 66.46804047, 67.88224793, 69.29646301])
        zC1_diagonal_50x50_CMG = np.array([1, 1, 1, 1, 1, 1, 0.999999881, 1, 0.999999762, \
        0.999988198, 0.999999821, 0.999993682, 1, 1, 0.999999642, 0.999990106, 1, 1, \
        0.99999994, 0.999998808, 1, 0.99999994, 1, 1, 0.999999881, 0.999996603, 0.999948263, \
        0.999491334, 0.996616483, 0.984175324, 0.946150005, 0.862054944, 0.723623037, \
        0.54798609, 0.369852424, 0.221005201, 0.116587371, 0.054480899, 0.022775447, \
        0.008617969, 0.002986676, 0.0009609, 0.000291741, 8.53E-05, 2.46E-05, 7.17E-06, \
        2.20E-06, 7.37E-07, 2.82E-07, 2.53E-08])

        x_diagonal_80x80_CMG = np.array([0, 0.883883476, 1.767766953, 2.651650429, \
        3.535533905, 4.419417381, 5.303300858, 6.187184334, 7.07106781, 7.954951286, \
        8.838834763, 9.722718239, 10.60660172, 11.49048519, 12.37436867, 13.25825214, \
        14.14213562, 15.0260191, 15.90990257, 16.793787, 17.67766953, 18.56155396, \
        19.44543648, 20.32932091, 21.21320343, 22.09708786, 22.98097038, 23.86485481, \
        24.74873734, 25.63262177, 26.51650429, 27.40038872, 28.28427124, 29.16815567, \
        30.05203819, 30.93592262, 31.81980515, 32.70368958, 33.58757401, 34.47145462, \
        34.91339874, 35.35533905, 36.23922348, 38.44893265, 40.65864182, 41.54252243, \
        42.42640686, 43.31029129, 44.19417572, 45.07805634, 45.96194077, 46.8458252, \
        47.72970963, 48.61359024, 49.49747467, 50.3813591, 51.26524353, 52.14912415, \
        53.03300858, 53.91689301, 54.80077744, 55.68465805, 56.56854248, 57.45242691, \
        58.33631134, 59.22019196, 60.10407639, 60.98796082, 61.87184525, 62.75572586, \
        63.63961029, 64.52349091, 65.40737915, 66.29125977, 67.17514801, 68.05902863, \
        68.94290924, 69.82679749, 70.7106781, 71.59455872, 72.47844696, 73.36232758])

        zC1_diagonal_80x80_CMG = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
        1, 0.999999762, 0.99999845, 0.999996006, 0.999990761, 0.999953985, 0.999996006, \
        0.999806643, 0.999304116, 0.997831643, 0.994088888, 0.985766292, 0.969451785, \
        0.941059291, 0.896887422, 0.8350088, 0.756337404, 0.664775074, 0.566333175, \
        0.467674226, 0.37471348, 0.291727841, 0.221075371, 0.163385257, 0.117995597, \
        0.083448485, 0.05792382, 0.039560262, 0.026657011, 0.01777532, 0.011767492, \
        0.00776055, 0.005116371, 0.003383839, 0.002252761, 0.001514557, 0.001031375, \
        0.000713245, 0.00039574, 0.000284534, 0.000231793, 0.000185183, 0.000146023, \
        0.000112076, 7.91E-05])

        datas = np.load('flying/results_Hoteit_Firoo_2k_ex1_IMPEC_FOU_5596.npy', allow_pickle=True)
        for data in datas[-1:]:
            zC1_40x40_FOU = data[10][0]
            centroid = data[11]
            centroid_x = centroid[:,0]
            centroid_y = centroid[:,1]
            ind = abs(centroid_x - (centroid_y))<1e-15
            zC1_diagonal_40x40_FOU = zC1_40x40_FOU[ind]
            x_diagonal_40x40_FOU = (((centroid_x[ind])**2)*2)**(1/2)
            
        datas = np.load('flying/results_Hoteit_Firoo_2k_ex1_IMPEC_MUSCL_upw_5560.npy', allow_pickle=True)
        for data in datas[-1:]:
            zC1_40x40_MUSCL_upw = data[10][0]
            centroid = data[11]
            centroid_x = centroid[:,0]
            centroid_y = centroid[:,1]
            ind = abs(centroid_x - (centroid_y))<1e-15
            zC1_diagonal_40x40_MUSCL_upw = zC1_40x40_MUSCL_upw[ind]
            x_diagonal_40x40_MUSCL_upw = (((centroid_x[ind])**2)*2)**(1/2)

        datas = np.load('flying/results_Hoteit_Firoo_2k_ex1_IMPEC_MUSCL_LLF_5619.npy', allow_pickle=True)
        for data in datas[-1:]:
            zC1_40x40_MUSCL_LLF = data[10][0]
            centroid = data[11]
            centroid_x = centroid[:,0]
            centroid_y = centroid[:,1]
            ind = abs(centroid_x - (centroid_y))<1e-15
            zC1_diagonal_40x40_MUSCL_LLF = zC1_40x40_MUSCL_LLF[ind]
            x_diagonal_40x40_MUSCL_LLF = (((centroid_x[ind])**2)*2)**(1/2)

        plt.figure(1)
        plt.plot(x_diagonal_40x40_FOU, zC1_diagonal_40x40_FOU, 'm')
        plt.plot(x_diagonal_40x40_MUSCL_upw, zC1_diagonal_40x40_MUSCL_upw, 'g')
        plt.plot(x_diagonal_50x50_CMG, zC1_diagonal_50x50_CMG, 'k')
        plt.plot(x_diagonal_40x40_MUSCL_LLF, zC1_diagonal_40x40_MUSCL_LLF, 'b')

        #plt.plot(x_diagonal_80x80_CMG, zC1_diagonal_80x80_CMG, 'b')
        plt.grid()
        plt.legend(('FOU 40x40', 'MUSCL+upw 40x40', 'CMG 50x50', 'MUSCL+LLF 40x40'))
        plt.title('Methane global composition profile')
        plt.ylabel('Methane global composition')
        plt.xlabel('Diagonal distance')
        plt.savefig('results/compositional/TCC2/zC1_Firoo_ex1_40x40.png')
        import pdb; pdb.set_trace()
