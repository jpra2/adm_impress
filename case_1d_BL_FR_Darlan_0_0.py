import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline
import quadpy

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)


xD = np.loadtxt('case_plots/x_BL_Darlan_semi_analytical.txt')
SwD = np.loadtxt('case_plots/Sw_BL_Darlan_semi_analytical.txt')
L = 0.6096
xD*=L
f = interp1d(xD,SwD)
ksi_w = 999.55/18.015e-3

for  arq in arquivos:
    if  arq.startswith(name):

        '''--------------FLUX RECONSTRUCTION RESULTS (2nd order)-------------'''
        # MLPu1 mod

        #datas = np.load('flying/results_BL_Darlan_32_FR2v_191.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_32_FR2_373.npy', allow_pickle=True) #159
        for data in datas[datas.shape[0]-1:]:
            Nk32_FR2 = data[12][-1].flatten()
            n = 32
            Sw32_FR2 = Nk32_FR2 / ksi_w /(0.6096/n*0.03048**2)
            x32_2 = np.empty((n,2))
            GL = quadpy.c1.gauss_lobatto(2)
            points = GL.points
            for i in range(len(points)):
                x32_2[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)*L
            x32_2 = x32_2.flatten()
            Sw32_ans_2 = np.empty_like(Nk32_FR2)
            for i in range(n*2):
                Sw32_ans_2[i] = f(x32_2[i])

            e32_L1_FR2_t = (sum(abs(Sw32_ans_2-Sw32_FR2))*(L/n))

            Sw32_FR_u1mod = data[5]
            Nk32_FR = data[12][0]
            n = 32
            x32 = np.linspace(0+1/64, 1-1/64, 32)*L
            t_32_FR2 = data[2]
            e32_L1_FR2 = (sum(abs(f(x32)-Sw32_FR_u1mod))*(L/n))
            e32_L2_FR2 = np.sqrt(np.sum((f(x32)-Sw32_FR_u1mod)**2) * L/n)
            t_32_acum_FR2 = data[16]
            t_step_32_FR2 = np.copy(t_32_acum_FR2)
            t_step_32_FR2[0] = 0
            for i in range(1,len(t_32_acum_FR2)):
                t_step_32_FR2[i] = t_32_acum_FR2[i] - t_32_acum_FR2[i-1]

        #datas = np.load('flying/results_BL_Darlan_64_FR2v_378.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_64_FR2_77.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw64_FR_u1mod = data[5]
            Nk64_FR = data[12][0]
            n = 64
            x64 = np.linspace(0+1/128, 1-1/128, 64)*L
            t_64_FR2 = data[2]
            e64_L1_FR2 = (sum(abs(f(x64)-Sw64_FR_u1mod))*(L/n))
            e64_L2_FR2 = np.sqrt(np.sum((f(x64)-Sw64_FR_u1mod)**2) * L/n)

            t_64_acum_FR2 = data[16]
            t_step_64_FR2 = np.copy(t_64_acum_FR2)
            t_step_64_FR2[0] = 0
            for i in range(1,len(t_64_acum_FR2)):
                t_step_64_FR2[i] = t_64_acum_FR2[i] - t_64_acum_FR2[i-1]

        #datas = np.load('flying/results_BL_Darlan_128_FR2v_746.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_128_FR2_1422.npy', allow_pickle=True) #590
        for data in datas[datas.shape[0]-1:]:
            Sw128_FR_u1mod = data[5]
            Nk128_FR = data[12][0]
            n = 128
            x128 = np.linspace(0+1/256, 1-1/256, 128)*L
            t_128_FR2 = data[2]
            e128_L1_FR2 = (sum(abs(f(x128)-Sw128_FR_u1mod))*(L/n))
            e128_L2_FR2 = np.sqrt(np.sum((f(x128)-Sw128_FR_u1mod)**2) * L/n)

            t_128_acum_FR2 = data[16]
            t_step_128_FR2 = np.copy(t_128_acum_FR2)
            t_step_128_FR2[0] = 0
            for i in range(1,len(t_128_acum_FR2)):
                t_step_128_FR2[i] = t_128_acum_FR2[i] - t_128_acum_FR2[i-1]

        datas = np.load('flying/results_BL_Darlan_256_FR2_2832.npy', allow_pickle=True) #1172
        #datas = np.load('flying/results_BL_Darlan_256_FR2v_1486.npy', allow_pickle=True) #1172
        for data in datas[datas.shape[0]-1:]:
            Sw256_FR_u1mod = data[5]
            Nk256_FR = data[12][0]
            n = 256
            x256 = np.linspace(0+1/512, 1-1/512, 256)*L
            t_256_FR2 = data[2]
            e256_L1_FR2 = (sum(abs(f(x256)-Sw256_FR_u1mod))*(L/n))
            e256_L2_FR2 = np.sqrt(np.sum((f(x256)-Sw256_FR_u1mod)**2) * L/n)

            t_256_acum_FR2 = data[16]
            t_step_256_FR2 = np.copy(t_256_acum_FR2)
            t_step_256_FR2[0] = 0
            for i in range(1,len(t_256_acum_FR2)):
                t_step_256_FR2[i] = t_256_acum_FR2[i] - t_256_acum_FR2[i-1]

        '''--------------FLUX RECONSTRUCTION RESULTS (3rd order)-------------'''
        # MLPu1 mod

        #datas = np.load('flying/results_BL_Darlan_32_FR3v_285.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_32_FR3_626.npy', allow_pickle=True) #273
        for data in datas[datas.shape[0]-1:]:
            Nk32_FR3 = data[12][-1].flatten()
            n = 32
            Sw32_FR3 = Nk32_FR3 / ksi_w /(0.6096/n*0.03048**2)
            x32_3 = np.empty((n,3))
            GL = quadpy.c1.gauss_lobatto(3)
            points = GL.points
            for i in range(len(points)):
                x32_3[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)*L
            x32_3 = x32_3.flatten()
            Sw32_ans_3 = np.empty_like(Nk32_FR3)
            for i in range(n*3):
                Sw32_ans_3[i] = f(x32_3[i])

            e32_L1_FR3_t = (sum(abs(Sw32_ans_3-Sw32_FR3))*(L/n))
            Sw32_FR3_u1mod = data[5]
            t_32_FR3 = data[2]
            t_32_FR3 = data[2]
            n=32
            e32_L1_FR3 = (sum(abs(f(x32)-Sw32_FR3_u1mod))*(L/n))
            e32_L2_FR3 = np.sqrt(np.sum((f(x32)-Sw32_FR3_u1mod)**2) * L / n)

            t_32_acum_FR3 = data[16]
            t_step_32_FR3 = np.copy(t_32_acum_FR3)
            t_step_32_FR3[0] = 0
            for i in range(1,len(t_32_acum_FR3)):
                t_step_32_FR3[i] = t_32_acum_FR3[i] - t_32_acum_FR3[i-1]

        #datas = np.load('flying/results_BL_Darlan_64_FR3v_571.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_64_FR3_1171.npy', allow_pickle=True) #555
        for data in datas[datas.shape[0]-1:]:
            Sw64_FR3_u1mod = data[5]
            t_64_FR3 = data[2]
            n=64
            e64_L1_FR3 = (sum(abs(f(x64)-Sw64_FR3_u1mod))*(L/64))
            e64_L2_FR3 = np.sqrt(np.sum((f(x64)-Sw64_FR3_u1mod)**2) * L / 64)

            t_64_acum_FR3 = data[16]
            t_step_64_FR3 = np.copy(t_64_acum_FR3)
            t_step_64_FR3[0] = 0
            for i in range(1,len(t_64_acum_FR3)):
                t_step_64_FR3[i] = t_64_acum_FR3[i] - t_64_acum_FR3[i-1]

        #datas = np.load('flying/results_BL_Darlan_128_FR3v_1092.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_128_FR3_3571.npy', allow_pickle=True) #1127
        for data in datas[datas.shape[0]-1:]:
            Sw128_FR3_u1mod = data[5]
            t_128_FR3 = data[2]
            e128_L1_FR3 = (sum(abs(f(x128)-Sw128_FR3_u1mod))*(L/128))
            e128_L2_FR3 = np.sqrt(np.sum((f(x128)-Sw128_FR3_u1mod)**2) * L / 128)

            t_128_acum_FR3 = data[16]
            t_step_128_FR3 = np.copy(t_128_acum_FR3)
            t_step_128_FR3[0] = 0
            for i in range(1,len(t_128_acum_FR3)):
                t_step_128_FR3[i] = t_128_acum_FR3[i] - t_128_acum_FR3[i-1]

        datas = np.load('flying/results_BL_Darlan_256_FR3_4714.npy', allow_pickle=True) #2282
        for data in datas[datas.shape[0]-1:]:
            Sw256_FR3_u1mod = data[5]
            t_256_FR3 = data[2]
            e256_L1_FR3 = (sum(abs(f(x256)-Sw256_FR3_u1mod))*(L/256))
            e256_L2_FR3 = np.sqrt(np.sum((f(x256)-Sw256_FR3_u1mod)**2) * L / 256)

            t_256_acum_FR3 = data[16]
            t_step_256_FR3 = np.copy(t_256_acum_FR3)
            t_step_256_FR3[0] = 0
            for i in range(1,len(t_256_acum_FR3)):
                t_step_256_FR3[i] = t_256_acum_FR3[i] - t_256_acum_FR3[i-1]

        '''--------------FLUX RECONSTRUCTION RESULTS (4th order)-------------'''
        # MLPu1 mod

        #datas = np.load('flying/results_BL_Darlan_32_FR4v_374.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_32_FR4_835.npy', allow_pickle=True) #376
        for data in datas[datas.shape[0]-1:]:
            Nk32_FR4 = data[12][-1].flatten()
            n = 32
            Sw32_FR4 = Nk32_FR4 / ksi_w /(0.6096/n*0.03048**2)
            x32_4 = np.empty((n,4))
            GL = quadpy.c1.gauss_lobatto(4)
            points = GL.points
            for i in range(len(points)):
                x32_4[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x32_4 = x32_4.flatten()*L
            Sw32_ans_4 = np.empty_like(Nk32_FR4)
            for i in range(n*4):
                Sw32_ans_4[i] = f(x32_4[i])

            e32_L1_FR4_t = (sum(abs(Sw32_ans_4-Sw32_FR4))*(L/n))

            Sw32_FR4_u1mod = data[5]
            t_32_FR4 = data[2]
            e32_L1_FR4 = (sum(abs(f(x32)-Sw32_FR4_u1mod))*(L/n))
            e32_L2_FR4 = np.sqrt(np.sum((f(x32)-Sw32_FR4_u1mod)**2) * L / n)

            t_32_acum_FR4 = data[16]
            t_step_32_FR4 = np.copy(t_32_acum_FR4)
            t_step_32_FR4[0] = 0
            for i in range(1,len(t_32_acum_FR4)):
                t_step_32_FR4[i] = t_32_acum_FR4[i] - t_32_acum_FR4[i-1]

        #datas = np.load('flying/results_BL_Darlan_64_FR4v_728.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_64_FR4_1626.npy', allow_pickle=True) #995
        for data in datas[datas.shape[0]-1:]:
            n = 64

            Sw64_FR4 = data[12][-1].flatten() / ksi_w /(0.6096/n*0.03048**2)
            x64_4 = np.empty((n,4))
            GL = quadpy.c1.gauss_lobatto(4)
            points = GL.points
            for i in range(len(points)):
                x64_4[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x64_4 = x64_4.flatten()
            Sw64_ans_4 = np.empty_like(x64_4)
            for i in range(n*4):
                Sw64_ans_4[i] = f(x64_4[i])

            e64_L1_FR4_t = (sum(abs(Sw64_ans_4-Sw64_FR4))*(L/n))

            Sw64_FR4_u1mod = data[5]
            t_64_FR4 = data[2]
            e64_L1_FR4 = (sum(abs(f(x64)-Sw64_FR4_u1mod))*(L/n))
            e64_L2_FR4 = np.sqrt(np.sum((f(x64)-Sw64_FR4_u1mod)**2) * L / n)
            t_64_acum_FR4 = data[16]
            t_step_64_FR4 = np.copy(t_64_acum_FR4)
            t_step_64_FR4[0] = 0
            for i in range(1,len(t_64_acum_FR4)):
                t_step_64_FR4[i] = t_64_acum_FR4[i] - t_64_acum_FR4[i-1]

        #datas = np.load('flying/results_BL_Darlan_128_FR4v_1488.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_128_FR4_4936.npy', allow_pickle=True) #1566
        for data in datas[datas.shape[0]-1:]:
            Sw128_FR4_u1mod = data[5]
            t_128_FR4 = data[2]
            e128_L1_FR4 = (sum(abs(f(x128)-Sw128_FR4_u1mod))*(L/128))
            e128_L2_FR4 = np.sqrt(np.sum((f(x128)-Sw128_FR4_u1mod)**2) * L / 128)
            t_128_acum_FR4 = data[16]
            t_step_128_FR4 = np.copy(t_128_acum_FR4)
            t_step_128_FR4[0] = 0
            for i in range(1,len(t_128_acum_FR4)):
                t_step_128_FR4[i] = t_128_acum_FR4[i] - t_128_acum_FR4[i-1]

        datas = np.load('flying/results_BL_Darlan_256_FR4_9933.npy', allow_pickle=True) #3156
        for data in datas[datas.shape[0]-1:]:
            Sw256_FR4_u1mod = data[5]
            t_256_FR4 = data[2]
            e256_L1_FR4 = (sum(abs(f(x256)-Sw256_FR4_u1mod))*(L/256))
            e256_L2_FR4 = np.sqrt(np.sum((f(x256)-Sw256_FR4_u1mod)**2) * L / 256)
            t_256_acum_FR4 = data[16]
            t_step_256_FR4 = np.copy(t_256_acum_FR4)
            t_step_256_FR4[0] = 0
            for i in range(1,len(t_256_acum_FR4)):
                t_step_256_FR4[i] = t_256_acum_FR4[i] - t_256_acum_FR4[i-1]

        '''----------------------FIRST ORDER UPWIND -------------------------'''

        #datas = np.load('flying/results_BL_Darlan_32_FOU_146.npy', allow_pickle=True) #73
        datas = np.load('flying/results_BL_Darlan_32_FOU_132.npy', allow_pickle=True) #73
        for data in datas[datas.shape[0]-1:]:
            Sw32_FOU = data[5]
            t_32_FOU = data[2]
            e32_L1_FOU = (sum(abs(f(x32)-Sw32_FOU))*(L/32))
            e32_L2_FOU = np.sqrt(np.sum((f(x32)-Sw32_FOU)**2) * L/ 32)
            t_32_acum_FOU = data[16]
            t_step_32_FOU = np.copy(t_32_acum_FOU)
            t_step_32_FOU[0] = 0
            for i in range(1,len(t_32_acum_FOU)):
                t_step_32_FOU[i] = t_32_acum_FOU[i] - t_32_acum_FOU[i-1]

        #datas = np.load('flying/results_BL_Darlan_64_FOU_279.npy', allow_pickle=True) #139
        datas = np.load('flying/results_BL_Darlan_64_FOU_73.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw64_FOU = data[5]
            t_64_FOU = data[2]
            e64_L1_FOU = (sum(abs(f(x64)-Sw64_FOU))*(L/64))
            e64_L2_FOU = np.sqrt(np.sum((f(x64)-Sw64_FOU)**2) * L / 64)
            R64_L1_FOU = math.log(e32_L1_FOU/e64_L1_FOU,2)
            t_64_acum_FOU = data[16]
            t_step_64_FOU = np.copy(t_64_acum_FOU)
            t_step_64_FOU[0] = 0
            for i in range(1,len(t_64_acum_FOU)):
                t_step_64_FOU[i] = t_64_acum_FOU[i] - t_64_acum_FOU[i-1]

        #datas = np.load('flying/results_BL_Darlan_128_FOU_542.npy', allow_pickle=True) #270
        datas = np.load('flying/results_BL_Darlan_128_FOU_145.npy', allow_pickle=True) #270
        for data in datas[datas.shape[0]-1:]:
            Sw128_FOU = data[5]
            t_128_FOU = data[2]
            e128_L1_FOU = (sum(abs(f(x128)-Sw128_FOU))*(L/128))
            e128_L2_FOU = np.sqrt(np.sum((f(x128)-Sw128_FOU)**2) * L / 128)
            R128_L1_FOU = math.log(e64_L1_FOU/e128_L1_FOU,2)

            t_128_acum_FOU = data[16]
            t_step_128_FOU = np.copy(t_128_acum_FOU)
            t_step_128_FOU[0] = 0
            for i in range(1,len(t_128_acum_FOU)):
                t_step_128_FOU[i] = t_128_acum_FOU[i] - t_128_acum_FOU[i-1]


        #datas = np.load('flying/results_BL_Darlan_256_FOU_1066.npy', allow_pickle=True) #532
        datas = np.load('flying/results_BL_Darlan_256_FOU_290.npy', allow_pickle=True) #532
        for data in datas[datas.shape[0]-1:]:
            Sw256_FOU = data[5]
            t_256_FOU = data[2]
            e256_L1_FOU = (sum(abs(f(x256)-Sw256_FOU))*(L/256))
            e256_L2_FOU = np.sqrt(np.sum((f(x256)-Sw256_FOU)**2) * L / 256)
            R256_L1_FOU = math.log(e128_L1_FOU/e256_L1_FOU,2)

            t_256_acum_FOU = data[16]
            t_step_256_FOU = np.copy(t_256_acum_FOU)
            t_step_256_FOU[0] = 0
            for i in range(1,len(t_256_acum_FOU)):
                t_step_256_FOU[i] = t_256_acum_FOU[i] - t_256_acum_FOU[i-1]


        '''----------------------------- MUSCL ------------------------------'''

        datas = np.load('flying/results_BL_Darlan_32_MUSCL_401.npy', allow_pickle=True) #85
        for data in datas[datas.shape[0]-1:]:
            Sw32_MUSCL = data[5]
            t_32_MUSCL = data[2]
            e32_L1_MUSCL = (sum(abs(f(x32)-Sw32_MUSCL))*(L/32))
            e32_L2_MUSCL = np.sqrt(np.sum((f(x32)-Sw32_MUSCL)**2) * L / 32)
            t_32_acum_MUSCL = data[16]
            t_step_32_MUSCL = np.copy(t_32_acum_MUSCL)
            t_step_32_MUSCL[0] = 0
            for i in range(1,len(t_32_acum_MUSCL)):
                t_step_32_MUSCL[i] = t_32_acum_MUSCL[i] - t_32_acum_MUSCL[i-1]

        datas = np.load('flying/results_BL_Darlan_64_MUSCL_225.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw64_MUSCL = data[5]
            t_64_MUSCL = data[2]
            e64_L1_MUSCL = (sum(abs(f(x64)-Sw64_MUSCL))*(L/64))
            e64_L2_MUSCL = np.sqrt(np.sum((f(x64)-Sw64_MUSCL)**2) * L / 64)
            t_64_acum_MUSCL = data[16]
            t_step_64_MUSCL = np.copy(t_64_acum_MUSCL)
            t_step_64_MUSCL[0] = 0
            for i in range(1,len(t_64_acum_MUSCL)):
                t_step_64_MUSCL[i] = t_64_acum_MUSCL[i] - t_64_acum_MUSCL[i-1]

        datas = np.load('flying/results_BL_Darlan_128_MUSCL_454.npy', allow_pickle=True)#677

        for data in datas[datas.shape[0]-1:]:
            Sw128_MUSCL = data[5]
            t_128_MUSCL = data[2]
            e128_L1_MUSCL = (sum(abs(f(x128)-Sw128_MUSCL))*(L/128))
            e128_L2_MUSCL = np.sqrt(np.sum((f(x128)-Sw128_MUSCL)**2) * L / 128)
            t_128_acum_MUSCL = data[16]
            t_step_128_MUSCL = np.copy(t_128_acum_MUSCL)
            t_step_128_MUSCL[0] = 0
            for i in range(1,len(t_step_128_MUSCL)):
                t_step_128_MUSCL[i] = t_128_acum_MUSCL[i] - t_128_acum_MUSCL[i-1]

        datas = np.load('flying/results_BL_Darlan_256_MUSCL_911.npy', allow_pickle=True) #1351
        for data in datas[datas.shape[0]-1:]:
            Sw256_MUSCL = data[5]
            t_256_MUSCL = data[2]
            e256_L1_MUSCL = (sum(abs(f(x256)-Sw256_MUSCL))*(L/256))
            e256_L2_MUSCL = np.sqrt(np.sum((f(x256)-Sw256_MUSCL)**2) * L / 256)
            t_256_acum_MUSCL = data[16]
            t_step_256_MUSCL = np.copy(t_256_acum_MUSCL)
            t_step_256_MUSCL[0] = 0
            for i in range(1,len(t_256_acum_MUSCL)):
                t_step_256_MUSCL[i] = t_256_acum_MUSCL[i] - t_256_acum_MUSCL[i-1]


        plt.figure(1)
        plt.plot(x64, Sw64_FOU, '-rs', \
            x64, Sw64_FR_u1mod, '-mP', \
            x64, Sw64_FR3_u1mod, '-yv', \
            x64, Sw64_FR4_u1mod, '-b<', \
            x64, Sw64_MUSCL, '-co', mfc='none')
        plt.plot( xD, SwD, '-k')
        plt.legend(('FOU', 'FR P1', 'FR P2', \
            'FR P3', 'MUSCL', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - 64 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_64FR.png', format='png')

        plt.figure(2)
        plt.plot(x32, Sw32_FOU, '-rs',\
            x32, Sw32_FR_u1mod, '-mP', \
            x32, Sw32_MUSCL, '-co',  mfc='none')
        plt.plot( xD, SwD, '-k')
        plt.legend(('FOU', 'FR P1', 'MUSCL', 'Analytical Solution'))
        plt.xlim(0.3, 0.5)
        plt.ylim(0.6, 0.75)
        plt.title('Buckley-Leverett Solution Example - 32 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_32FR.png', format='png')

        plt.figure(3)
        plt.plot(x128, Sw128_FOU, '-rs',\
            x128, Sw128_FR_u1mod, '-mP', \
            x128, Sw128_FR3_u1mod, '-yv', \
            x128, Sw128_FR4_u1mod, '-b<', \
            x128, Sw128_MUSCL, '-co',  mfc='none')
        plt.plot( xD, SwD, '-k')
        plt.legend(('FOU', 'FR P1', 'FR P2', \
            'FR P3', 'MUSCL', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - 128 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_128FR.png', format='png')

        plt.figure(4)
        plt.plot(x32, Sw32_FR_u1mod, '-mP', \
            x64, Sw64_FR_u1mod, '-yv', \
            x128, Sw128_FR_u1mod, '-b<',
            x256, Sw256_FR_u1mod, '-co', mfc='none')
        plt.plot( xD, SwD, '-k')
        plt.legend(('32 cells', '63 cells', '128 cells', '256 cells', 'Analytical Solution'))
        #plt.xlim(0.5, 0.8)
        #plt.ylim(0.6, 0.75)
        plt.title('Buckley-Leverett Solution Example - FR P1')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_FR2_conv.png', format='png')

        plt.figure(5)
        plt.plot(x32, Sw32_FR3_u1mod, '-mP', \
            x64, Sw64_FR3_u1mod, '-yv', \
            x128, Sw128_FR3_u1mod, '-b<',
            x256, Sw256_FR3_u1mod, '-co', mfc='none')
        plt.plot( xD, SwD, '-k')
        plt.legend(('32 cells', '63 cells', '128 cells', '256 cells','Analytical Solution'))
        #plt.xlim(0.5, 0.8)
        #plt.ylim(0.6, 0.75)
        plt.title('Buckley-Leverett Solution Example - FR P2')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_FR3_conv.png', format='png')

        plt.figure(6)
        plt.plot(x32, Sw32_FR4_u1mod, '-mP', \
            x64, Sw64_FR4_u1mod, '-yv', \
            x128, Sw128_FR4_u1mod, '-b<',
            x256, Sw256_FR4_u1mod, '-co', mfc='none')
        plt.plot( xD, SwD, '-k')
        plt.legend(('32 cells', '63 cells', '128 cells', '256 cells','Analytical Solution'))
        #plt.xlim(0.5, 0.8)
        #plt.ylim(0.6, 0.75)
        plt.title('Buckley-Leverett Solution Example - FR P3')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_FR4_conv.png', format='png')

        plt.figure(7)
        plt.plot(t_32_acum_FOU, t_step_32_FOU, '-y')
        plt.plot(t_32_acum_MUSCL, t_step_32_MUSCL, '-b')
        plt.plot(t_32_acum_FR2, t_step_32_FR2, '-g')
        plt.plot(t_32_acum_FR3, t_step_32_FR3, '-k')
        plt.plot(t_32_acum_FR4, t_step_32_FR4, '-c')
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.title('Time-step for Buckley-Leverett Example 32 cells')
        plt.ylabel('Time-step [s]')
        plt.xlabel('Time [s]')
        plt.savefig('results/compositional/FR/time_step_BL_Darlan_32comp.png', format='png')

        plt.figure(8)
        plt.plot(t_64_acum_FOU, t_step_64_FOU, '-y')
        plt.plot(t_64_acum_MUSCL, t_step_64_MUSCL, '-b')
        plt.plot(t_64_acum_FR2, t_step_64_FR2, '-g')
        plt.plot(t_64_acum_FR3, t_step_64_FR3, '-k')
        plt.plot(t_64_acum_FR4, t_step_64_FR4, '-c')
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.title('Time-step for Buckley-Leverett Example 64 cells')
        plt.ylabel('Time-step [s]')
        plt.xlabel('Time [s]')
        plt.savefig('results/compositional/FR/time_step_BL_Darlan_64comp.png', format='png')

        plt.figure(9)
        plt.plot(t_128_acum_FOU, t_step_128_FOU, '-y')
        plt.plot(t_128_acum_MUSCL, t_step_128_MUSCL, '-b')
        plt.plot(t_128_acum_FR2, t_step_128_FR2, '-g')
        plt.plot(t_128_acum_FR3, t_step_128_FR3, '-k')
        plt.plot(t_128_acum_FR4, t_step_128_FR4, '-c')
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.title('Time-step for Buckley-Leverett Example 128 cells')
        plt.ylabel('Time-step [s]')
        plt.xlabel('Time [s]')
        plt.savefig('results/compositional/FR/time_step_BL_Darlan_128comp.png', format='png')

        plt.figure(10)
        plt.plot(t_256_acum_FOU, t_step_256_FOU, '-y')
        plt.plot(t_256_acum_MUSCL, t_step_256_MUSCL, '-b')
        plt.plot(t_256_acum_FR2, t_step_256_FR2, '-g')
        plt.plot(t_256_acum_FR3, t_step_256_FR3, '-k')
        plt.plot(t_256_acum_FR4, t_step_256_FR4, '-c')
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.title('Time-step for Buckley-Leverett Example 256 cells')
        plt.ylabel('Time-step [s]')
        plt.xlabel('Time [s]')
        plt.savefig('results/compositional/FR/time_step_BL_Darlan_256comp.png', format='png')

        plt.figure(11)
        x = np.log10(np.array([32,64,128,256]))
        y_FR2 = np.log10(np.array([e32_L1_FR2, e64_L1_FR2, e128_L1_FR2,
            e256_L1_FR2]))
        y_FR3 = np.log10(np.array([e32_L1_FR3, e64_L1_FR3, e128_L1_FR3,
            e256_L1_FR3]))
        y_FR4 = np.log10(np.array([e32_L1_FR4, e64_L1_FR4, e128_L1_FR4,
            e256_L1_FR4]))
        y_MUSCL = np.log10(np.array([e32_L1_MUSCL, e64_L1_MUSCL,
            e128_L1_MUSCL, e256_L1_MUSCL]))
        y_FOU = np.log10(np.array([e32_L1_FOU, e64_L1_FOU, e128_L1_FOU,
            e256_L1_FOU]))

        y_ref = -x-0.5

        plt.plot(x, y_FOU, '-rs', x, y_FR2, '-mP', x, y_FR3, '-yv', x, y_FR4, '-b<', \
            x, y_MUSCL, '-co', mfc='none')
        plt.plot(x, y_ref, 'k')
        plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{10}({E}_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.legend(('FOU', 'FR-P1', 'FR-P2', 'FR-P3', 'MUSCL', 'Linear Convergence'))
        plt.grid()
        plt.savefig('results/compositional/FR/BL_Darlan_L1_convergence_FRs.png')

        plt.figure(12)
        plt.plot(x32, Sw32_FOU, '-r', x64, Sw64_FOU, '-g', x128, Sw128_FOU, '-b', \
            x256, Sw256_FOU, '-y', xD, SwD, 'k')
        plt.legend(('FOU 32', 'FOU 64', 'FOU 128', 'FOU 256', 'Reference'))
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_FOU')

        plt.figure(13)
        x = np.log10(np.array([32,64,128,256]))
        y_FR2 = np.log10(np.array([e32_L1_FR2, e64_L1_FR2, e128_L1_FR2,
            e256_L1_FR2]))
        y_FR3 = np.log10(np.array([e32_L1_FR3, e64_L1_FR3, e128_L1_FR3,
            e256_L1_FR3]))
        y_FR4 = np.log10(np.array([e32_L1_FR4, e64_L1_FR4, e128_L1_FR4,
            e256_L1_FR4]))

        y_ref = -x-0.8

        plt.plot(x, y_FR2, '-mP', x, y_FR3, '-yv', x, y_FR4, '-b<', \
                mfc='none')
        plt.plot(x, y_ref, 'k')
        plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{10}({E}_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.legend(('FR-P1', 'FR-P2', 'FR-P3','Linear Convergence'))
        plt.grid()
        plt.savefig('results/compositional/FR/BL_Darlan_L1_convergence_FRs_only.png')

        import pdb; pdb.set_trace()
