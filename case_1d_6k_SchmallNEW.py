import numpy as np
import matplotlib.pyplot as plt
import os
import math
import csv

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)

for arq in arquivos:
    if  arq.startswith(name):
        fx = open('x_points_CMG_SchmallNEW.txt','r')
        x_CMG = [float(line.rstrip('\n\r')) for line in fx]

        #fP = open('P1024.txt','r')
        #P_CMG = [float(line.rstrip('\n\r')) for line in fP]

        fSw = open('Sw_points_CMG_SchmallNEW.txt','r')
        Sw_CMG = [float(line.rstrip('\n\r')) for line in fSw]

        fSo = open('So_points_CMG_SchmallNEW.txt','r')
        So_CMG = [float(line.rstrip('\n\r')) for line in fSo]

        fSg = open('Sg_points_CMG_SchmallNEW.txt','r')
        Sg_CMG = [float(line.rstrip('\n\r')) for line in fSg]

        #datas = np.load('flying/results_water_inj_6k_modified_case_upw_4326.npy', allow_pickle=True)
        '''--------------------------------FOU-------------------------------'''
        datas = np.load('flying/results_6k_SchmallNEW_8_IMPEC_FOU_561.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_FOU_8 = data[5]
            So_FOU_8 = data[6]
            Sg_FOU_8 = data[7]
            z_FOU_8 = data[10]
            Nk_FOU_8 = data[12]
            p_FOU_8 = data[4]/1e3
            t_FOU_8 = data[3]
            n=8
            x_8 = np.linspace(50/(2*n),50*(1-1/(2*n)),n)

        '''-------------------------------MUSCL------------------------------'''
        datas = np.load('flying/results_6kSchmallNEW_8_IMPEC_MUSCL_minmod_RK3_553.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_MUSCL_8 = data[5]
            So_MUSCL_8 = data[6]
            Sg_MUSCL_8 = data[7]
            z_MUSCL_8 = data[10]
            Nk_MUSCL_8 = data[12]
            p_MUSCL_8 = data[4]/1e3
            t_MUSCL_8 = data[3]

        datas = np.load('flying/results_6k_SchmallNEW_16_IMPEC_MUSCL_VA_LLF_RK3_274.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_MUSCL_16 = data[5]
            So_MUSCL_16 = data[6]
            Sg_MUSCL_16 = data[7]
            z_MUSCL_16 = data[10]
            Nk_MUSCL_16 = data[12]
            p_MUSCL_16 = data[4]/1e3
            t_MUSCL_16 = data[3]

        '''--------------------------------FR2-------------------------------'''
        datas = np.load('flying/results_6k_SchmallNEW_8_IMPEC_FR2_LLF_RK3t_104.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_FR2_8 = data[5]
            So_FR2_8 = data[6]
            Sg_FR2_8 = data[7]
            z_FR2_8 = data[10]
            Nk_FR2_8 = data[12]
            p_FR2_8 = data[4]/1e3
            t_FR2_8 = data[3]

        datas = np.load('flying/results_6k_SchmallNEW_16_IMPEC_FR2_LLF_RK3_252.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_FR2_16 = data[5]
            So_FR2_16 = data[6]
            Sg_FR2_16 = data[7]
            z_FR2_16 = data[10]
            Nk_FR2_16 = data[12]
            p_FR2_16 = data[4]/1e3
            t_FR2_16 = data[3]


        '''--------------------------------FR3-------------------------------'''
        datas = np.load('flying/results_6k_SchmallNEW_8_IMPEC_FR3_LLF_RK3t_201.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_FR3_8 = data[5]
            So_FR3_8 = data[6]
            Sg_FR3_8 = data[7]
            z_FR3_8 = data[10]
            Nk_FR3_8 = data[12]
            p_FR3_8 = data[4]/1e3
            t_FR3_8 = data[3]

        datas = np.load('flying/results_6k_SchmallNEW_16_IMPEC_FR3_LLF_RK3_455.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_FR3_16 = data[5]
            So_FR3_16 = data[6]
            Sg_FR3_16 = data[7]
            z_FR3_16 = data[10]
            Nk_FR3_16 = data[12]
            p_FR3_16 = data[4]/1e3
            t_FR3_16 = data[3]

        '''--------------------------------FR4-------------------------------'''
        datas = np.load('flying/results_6k_SchmallNEW_8_IMPEC_FR4_LLF_RK3t_268.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_FR4_8 = data[5]
            So_FR4_8 = data[6]
            Sg_FR4_8 = data[7]
            z_FR4_8 = data[10]
            Nk_FR4_8 = data[12]
            p_FR4_8 = data[4]/1e3
            t_FR4_8 = data[3]

        datas = np.load('flying/results_6k_SchmallNEW_16_IMPEC_FR4_LLF_RK3t_602.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_FR4_16 = data[5]
            So_FR4_16 = data[6]
            Sg_FR4_16 = data[7]
            z_FR4_16 = data[10]
            Nk_FR4_16 = data[12]
            p_FR4_16 = data[4]/1e3
            t_FR4_16 = data[3]
            n=16
            x_16 = np.linspace(50/(2*n),50*(1-1/(2*n)),n)


        plt.figure(1)
        plt.title('t = 200 days - 8x1x1 mesh')
        plt.plot(x_CMG, Sw_CMG, 'k')
        plt.plot(x_8, Sw_FOU_8, '-bo', mfc='none')
        plt.plot(x_8, Sw_MUSCL_8, '-mv', mfc='none')
        plt.plot(x_8, Sw_FR2_8, '-gs', mfc='none')
        plt.plot(x_8, Sw_FR3_8, '-yv', mfc='none')
        plt.plot(x_8, Sw_FR4_8, '-r*', mfc='none')
        plt.ylabel('$S_w$')
        plt.xlabel('Distance')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.savefig('results/compositional/SPE_paper/6k_Sw_8.png')
        plt.close()

        plt.figure(1)
        plt.title('t = 200 days - 8x1x1 mesh')
        plt.plot(x_CMG, Sw_CMG, 'k')
        #plt.plot(x_8, Sw_FOU_8, '-bo', mfc='none')
        plt.plot(x_16, Sw_MUSCL_16, '-mv', mfc='none')
        plt.plot(x_16, Sw_FR2_16, '-gs', mfc='none')
        plt.plot(x_16, Sw_FR3_16, '-yv', mfc='none')
        plt.plot(x_16, Sw_FR4_16, '-r*', mfc='none')
        plt.ylabel('$S_w$')
        plt.xlabel('Distance')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.savefig('results/compositional/SPE_paper/6k_Sw_16.png')
        plt.close()

        plt.figure(1)
        plt.title('t = 200 days - 8x1x1 mesh')
        plt.plot(x_CMG, So_CMG, 'k')
        plt.plot(x_8, So_FOU_8, '-bo', mfc='none')
        plt.plot(x_8, So_MUSCL_8, '-mv', mfc='none')
        plt.plot(x_8, So_FR2_8, '-gs', mfc='none')
        plt.plot(x_8, So_FR3_8, '-yv', mfc='none')
        plt.plot(x_8, So_FR4_8, '-r*', mfc='none')
        plt.ylabel('$S_o$')
        plt.xlabel('Distance')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.savefig('results/compositional/SPE_paper/6k_So_8.png')
        plt.close()

        plt.figure(1)
        plt.title('t = 200 days - 8x1x1 mesh')
        plt.plot(x_CMG, So_CMG, 'k')
        #plt.plot(x_16, So_FOU_16, '-bo', mfc='none')
        plt.plot(x_16, So_MUSCL_16, '-mv', mfc='none')
        plt.plot(x_16, So_FR2_16, '-gs', mfc='none')
        plt.plot(x_16, So_FR3_16, '-yv', mfc='none')
        plt.plot(x_16, So_FR4_16, '-r*', mfc='none')
        plt.ylabel('$S_o$')
        plt.xlabel('Distance')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.savefig('results/compositional/SPE_paper/6k_So_16.png')
        plt.close()

        plt.figure(1)
        plt.title('t = 200 days - 8x1x1 mesh')
        #plt.plot(x_CMG, So_CMG, 'k')
        plt.plot(x_8, Nk_FOU_8[-1,:], '-bo', mfc='none')
        plt.plot(x_8, Nk_MUSCL_8[-1,:], '-mv', mfc='none')
        plt.plot(x_8, Nk_FR2_8[-1,:], '-gs', mfc='none')
        plt.ylabel('$N_w$')
        plt.xlabel('Distance')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1'))
        plt.grid()
        plt.savefig('results/compositional/SPE_paper/6k_Nw_8.png')
        plt.close()

        plt.figure(1)
        plt.title('t = 200 days - 8x1x1 mesh')
        #plt.plot(x_CMG, So_CMG, 'k')
        plt.plot(x_8, z_FOU_8[0,:], '-bo', mfc='none')
        plt.plot(x_8, z_MUSCL_8[0,:], '-mv', mfc='none')
        plt.plot(x_8, z_FR2_8[0,:], '-gs', mfc='none')
        plt.plot(x_8, z_FR3_8[0,:], '-y<', mfc='none')
        plt.plot(x_8, z_FR4_8[0,:], '-r*', mfc='none')
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.legend(( 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.savefig('results/compositional/SPE_paper/6k_zC1_8.png')
        plt.close()

        plt.figure(1)
        plt.title('t = 200 days - 8x1x1 mesh')
        #plt.plot(x_CMG, So_CMG, 'k')
        plt.plot(x_8, z_FOU_8[0,:], '-bo', mfc='none')
        plt.plot(x_16, z_MUSCL_16[0,:], '-mv', mfc='none')
        plt.plot(x_16, z_FR2_16[0,:], '-gs', mfc='none')
        plt.plot(x_16, z_FR3_16[0,:], '-y<', mfc='none')
        plt.plot(x_16, z_FR4_16[0,:], '-r*', mfc='none')
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.legend(( 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.savefig('results/compositional/SPE_paper/6k_zC1_16.png')
        plt.close()


        import pdb; pdb.set_trace()
