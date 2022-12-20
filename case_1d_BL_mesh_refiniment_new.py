import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
xD = np.loadtxt('case_plots/x_BL_semi_analytical.txt')
xD[-1] = 0.0
SwD = np.loadtxt('case_plots/Sw_BL_semi_analytical.txt')
f = interp1d(xD,SwD)

for arq in arquivos:
    if  arq.startswith(name):

        datas5 = np.load('flying/results_Buckley_Leverett_case_8_FI_NEW_36.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas5[1:]:
            Sw_FI_8 = data[5]
            So_FI_8 = data[6]
            Sg_FI_8 = data[7]
            Oil_p_FI_8 = data[8]
            Gas_p_FI_8 = data[9]
            pressure_FI_8 = data[4]/1e3
            time_FI_8 = data[3]
            x8 = np.linspace(0.0, 0.6096, 8)
            x8 = x8 / 0.6096

            e8_L1_FI = (sum(abs(f(x8)-Sw_FI_8))*(1/8))
            #R8_L1_FI = math.log(e8_L1_FI/e16_L1_FI,2)
            e8_L2_FI = np.sqrt(np.sum((f(x8)-Sw_FI_8)**2) * 1 / 8)
            #R8_L2_FI = math.log(e8_L2_FI/e16_L2_FI,2)

        datas4 = np.load('flying/results_Buckley_Leverett_case_16_FI_NEW_36.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas4[1:]:
            Sw_FI_16 = data[5]
            So_FI_16 = data[6]
            Sg_FI_16 = data[7]
            Oil_p_FI_16 = data[8]
            Gas_p_FI_16 = data[9]
            pressure_FI_16 = data[4]/1e3
            time_FI_16 = data[3]
            x16 = np.linspace(0.0, 0.6096, 16)
            x16 = x16 / 0.6096

            e16_L1_FI = (sum(abs(f(x16)-Sw_FI_16))*(1/16))
            R16_L1_FI = math.log(e8_L1_FI/e16_L1_FI,2)
            e16_L2_FI = np.sqrt(np.sum((f(x16)-Sw_FI_16)**2) * 1 / 16)
            R16_L2_FI = math.log(e8_L2_FI/e16_L2_FI,2)

        datas3 = np.load('flying/results_Buckley_Leverett_case_32_FI_NEW_44.npy', allow_pickle=True)
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

            e32_L1_FI = (sum(abs(f(x32)-Sw_FI_32))*(1/32))
            R32_L1_FI = math.log(e16_L1_FI/e32_L1_FI,2)
            e32_L2_FI = np.sqrt(np.sum((f(x32)-Sw_FI_32)**2) * 1 / 32)
            R32_L2_FI = math.log(e16_L2_FI/e32_L2_FI,2)

        datas2 = np.load('flying/BL_Teste2/results_Buckley_Leverett_case_64_FI_CFLMAIOR_70.npy', allow_pickle=True)
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

            e64_L1_FI = (sum(abs(f(x64)-Sw_FI_64))*(1/64))
            R64_L1_FI = math.log(e32_L1_FI/e64_L1_FI,2)
            e64_L2_FI = np.sqrt(np.sum((f(x64)-Sw_FI_64)**2) * 1 / 64)
            R64_L2_FI = math.log(e32_L2_FI/e64_L2_FI,2)

        datas = np.load('flying/BL_Teste2/results_Buckley_Leverett_case_128_FI_NOVO_139.npy', allow_pickle=True)
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
            #import pdb; pdb.set_trace()
            e128_L1_FI = (sum(abs(f(x128)-Sw_FI_128))*(1/128))
            R128_L1_FI = math.log(e64_L1_FI/e128_L1_FI,2)
            e128_L2_FI = np.sqrt(np.sum((f(x128)-Sw_FI_128)**2) * 1 / 128)
            R128_L2_FI = math.log(e64_L2_FI/e128_L2_FI,2)

        datas = np.load('flying/results_Buckley_Leverett_case_256_FI_1384.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FI_256 = data[5]
            So_FI_256 = data[6]
            Sg_FI_256 = data[7]
            Oil_p_FI_256 = data[8]
            Gas_p_FI_256 = data[9]
            pressure_FI_256 = data[4]/1e3
            time_FI_256 = data[3]
            x256 = np.linspace(0.0, 0.6096, 256)
            x256 = x256/0.6096
            #import pdb; pdb.set_trace()
            e256_L1_FI = (sum(abs(f(x256)-Sw_FI_256))*(1/256))
            R256_L1_FI = math.log(e128_L1_FI/e256_L1_FI,2)
            e256_L2_FI = np.sqrt(np.sum((f(x256)-Sw_FI_256)**2) * 1 / 256)
            R256_L2_FI = math.log(e128_L2_FI/e256_L2_FI,2)

        datas = np.load('flying/results_Buckley_Leverett_case_512_FI_693.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw_FI_512 = data[5]
            So_FI_512 = data[6]
            Sg_FI_512 = data[7]
            Oil_p_FI_512 = data[8]
            Gas_p_FI_512 = data[9]
            pressure_FI_512 = data[4]/1e3
            time_FI_512 = data[3]
            x512 = np.linspace(0.0, 0.6096, 512)
            x512 = x512/0.6096

            e512_L1_FI = (sum(abs(f(x512)-Sw_FI_512))*(1/512))
            R512_L1_FI = math.log(e256_L1_FI/e512_L1_FI,2)
            e512_L2_FI = np.sqrt(np.sum((f(x512)-Sw_FI_512)**2) * 1 / 512)
            R512_L2_FI = math.log(e256_L2_FI/e512_L2_FI,2)

        #import pdb; pdb.set_trace()
        plt.figure(1)
        #plt.title('BL Sw - mesh refinement')
        plt.plot(xD, SwD, 'b')
        plt.plot(x64, Sw_FI_64, '.-c')
        plt.plot(x128, Sw_FI_128, '-.k')
        plt.plot(x256, Sw_FI_256, ':g')
        plt.plot(x512, Sw_FI_512, '--r')
        plt.legend(('Solução analítica', 'FI 64 CVs', 'FI 128 CVs', 'FI 256 CVs', 'FI 512 CVs'))
        #plt.legend(('IMPEC', 'Fully Implicit - back', 'Analytical Solution', 'Fully Implicit - new'))
        plt.ylabel('Saturação de água')
        plt.xlabel('Distância (m)')
        plt.grid()
        plt.savefig('results/BL_Sw_refinement_3' + '.png')


        plt.figure(2)
        x = np.log2(np.array([8, 16, 32, 64, 128, 256, 512]))
        y_FI = np.log2(np.array([e8_L1_FI, e16_L1_FI, e32_L1_FI, e64_L1_FI, e128_L1_FI, e256_L1_FI, e512_L1_FI]))

        y_ref = -x-2.0

        plt.plot(x, y_FI, '-mP', mfc='none')
        plt.plot(x, y_ref, 'k')
        #plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{2}({E}_{L_1})$')
        plt.xlabel('$log_{2}(N)$')
        plt.legend(('FI', 'Convergência linear'))
        plt.grid()
        plt.savefig('results/BL_L1_convergence_FI_2.png')

        plt.figure(3)
        x = np.log2(np.array([8,16,32,64,128,256,512]))
        y = np.log2(np.array([e8_L2_FI, e16_L2_FI, e32_L2_FI, e64_L2_FI, e128_L2_FI, e256_L2_FI, e512_L2_FI]))
        plt.plot(x, y,'-bo')

        ref_line = x[0:]/2
        #plt.plot(x[0:2], -np.array([ref_line[1],ref_line[1]])-5,'-k')
        #plt.plot(np.array([x[0],x[0]]),-ref_line-5,'-k')
        plt.plot(x[0:],-ref_line-3.5,'-k')
        plt.plot()
        plt.legend(('FI','2$nd$ order'))
        plt.ylabel('$log_2({E}_{L_2})$')
        plt.xlabel('$log_2(N)$')
        plt.grid()
        plt.savefig('results/BL_L2_convergence_FI.png')

        import pdb; pdb.set_trace()
