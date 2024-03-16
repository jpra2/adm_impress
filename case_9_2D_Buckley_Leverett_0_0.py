import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import InterpolatedUnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)

for  arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_Buckley_Leverett_case_2D_20x20_IMPEC_66.npy', allow_pickle=True)
        for data in datas[-1:]:
            t_20x20_IMPEC = data[2]
            t_vec_20x20_IMPEC = np.array(data[16])/86400
            time_step_size_20x20_IMPEC = np.copy(t_vec_20x20_IMPEC)
            for i in range(len(data[16])):
                if i==0:
                    time_step_size_20x20_IMPEC[i] = t_vec_20x20_IMPEC[i]
                else: time_step_size_20x20_IMPEC[i] = t_vec_20x20_IMPEC[i] - t_vec_20x20_IMPEC[i-1]
            time_step_size_20x20_IMPEC = time_step_size_20x20_IMPEC*1e3

        datas = np.load('flying/results_Buckley_Leverett_case_2D_20x20_IMPSAT_30.npy', allow_pickle=True)
        for data in datas[-1:]:
            t_20x20_IMPSAT = data[2]
            t_vec_20x20_IMPSAT = np.array(data[16])/86400
            time_step_size_20x20_IMPSAT = np.copy(t_vec_20x20_IMPSAT)

            for i in range(len(data[16])):
                if i==0:
                    time_step_size_20x20_IMPSAT[i] = t_vec_20x20_IMPSAT[i]
                else: time_step_size_20x20_IMPSAT[i] = t_vec_20x20_IMPSAT[i] - t_vec_20x20_IMPSAT[i-1]
            time_step_size_20x20_IMPSAT = time_step_size_20x20_IMPSAT*1e3


        datas = np.load('flying/results_Buckley_Leverett_case_2D_40x40_IMPEC_132.npy', allow_pickle=True)
        for data in datas[-1:]:
            t_40x40_IMPEC = data[2]
            t_vec_40x40_IMPEC = np.array(data[16])/86400
            time_step_size_40x40_IMPEC = np.copy(t_vec_40x40_IMPEC)
            for i in range(len(data[16])):
                if i==0:
                    time_step_size_40x40_IMPEC[i] = t_vec_40x40_IMPEC[i]
                else: time_step_size_40x40_IMPEC[i] = t_vec_40x40_IMPEC[i] - t_vec_40x40_IMPEC[i-1]
            time_step_size_40x40_IMPEC = time_step_size_40x40_IMPEC*1e3

        datas = np.load('flying/results_Buckley_Leverett_case_2D_40x40_IMPSAT_64.npy', allow_pickle=True)
        for data in datas[-1:]:
            t_40x40_IMPSAT = data[2]
            t_vec_40x40_IMPSAT = np.array(data[16])/86400
            time_step_size_40x40_IMPSAT = np.copy(t_vec_40x40_IMPSAT)

            for i in range(len(data[16])):
                if i==0:
                    time_step_size_40x40_IMPSAT[i] = t_vec_40x40_IMPSAT[i]
                else: time_step_size_40x40_IMPSAT[i] = t_vec_40x40_IMPSAT[i] - t_vec_40x40_IMPSAT[i-1]
            time_step_size_40x40_IMPSAT = time_step_size_40x40_IMPSAT*1e3


        datas = np.load('flying/results_Buckley_Leverett_case_2D_80x80_IMPEC_663.npy', allow_pickle=True)
        for data in datas[-1:]:
            t_80x80_IMPEC = data[2]
            t_vec_80x80_IMPEC = np.array(data[16])/86400
            time_step_size_80x80_IMPEC = np.copy(t_vec_80x80_IMPEC)
            for i in range(len(data[16])):
                if i==0:
                    time_step_size_80x80_IMPEC[i] = t_vec_80x80_IMPEC[i]
                else: time_step_size_80x80_IMPEC[i] = t_vec_80x80_IMPEC[i] - t_vec_80x80_IMPEC[i-1]
            time_step_size_80x80_IMPEC = time_step_size_80x80_IMPEC*1e3

        datas = np.load('flying/results_Buckley_Leverett_case_2D_80x80_IMPSAT_205.npy', allow_pickle=True)
        for data in datas[-1:]:
            t_80x80_IMPSAT = data[2]
            t_vec_80x80_IMPSAT = np.array(data[16])/86400
            time_step_size_80x80_IMPSAT = np.copy(t_vec_80x80_IMPSAT)

            for i in range(len(data[16])):
                if i==0:
                    time_step_size_80x80_IMPSAT[i] = t_vec_80x80_IMPSAT[i]
                else: time_step_size_80x80_IMPSAT[i] = t_vec_80x80_IMPSAT[i] - t_vec_80x80_IMPSAT[i-1]
            time_step_size_80x80_IMPSAT = time_step_size_80x80_IMPSAT*1e3

        plt.figure(1)
        plt.plot(t_vec_20x20_IMPEC, time_step_size_20x20_IMPEC, 'm')
        plt.plot(t_vec_20x20_IMPSAT, time_step_size_20x20_IMPSAT, 'b')
        plt.grid()
        plt.legend(('IMPEC 20x20', 'IMPSAT 20x20'))
        plt.title('Análise do passo de tempo')
        plt.ylabel('Tamanho do passo de tempo [10$^{-3}$ dia]')
        plt.xlabel('Tempo [dia]')
        plt.savefig('results/compositional/time_stepsBL_IMPSAT20.png')

        plt.figure(2)
        plt.plot(t_vec_40x40_IMPEC, time_step_size_40x40_IMPEC, 'm')
        plt.plot(t_vec_40x40_IMPSAT, time_step_size_40x40_IMPSAT, 'b')
        plt.grid()
        plt.legend(('IMPEC 40x40', 'IMPSAT 40x40'))
        plt.title('Análise do passo de tempo')
        plt.ylabel('Tamanho do passo de tempo [10$^{-3}$ dia]')
        plt.xlabel('Tempo [dia]')
        plt.savefig('results/compositional/time_stepsBL_IMPSAT40.png')

        plt.figure(3)
        plt.plot(t_vec_80x80_IMPEC, time_step_size_80x80_IMPEC, 'm')
        plt.plot(t_vec_80x80_IMPSAT, time_step_size_80x80_IMPSAT, 'b')
        plt.grid()
        plt.legend(('IMPEC 80x80', 'IMPSAT 80x80'))
        plt.title('Análise do passo de tempo')
        plt.ylabel('Tamanho do passo de tempo [10$^{-3}$ dia]')
        plt.xlabel('Tempo [dia]')
        plt.savefig('results/compositional/time_stepsBL_IMPSAT80.png')
        import pdb; pdb.set_trace()
