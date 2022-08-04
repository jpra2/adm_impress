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

        fzC1 = open('zC1_points_CMG_SchmallNEW.txt','r')
        zC1_CMG = [float(line.rstrip('\n\r')) for line in fzC1]

        #datas = np.load('flying/results_water_inj_6k_modified_case_upw_4326.npy', allow_pickle=True)
        '''--------------------------------FOU-------------------------------'''

        datas = np.load('flying/results_6k_SchmallNEW_50_IMPEC_FOU_132.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_FOU_50 = data[5]
            So_FOU_50 = data[6]
            Sg_FOU_50 = data[7]
            z_FOU_50 = data[10]
            Nk_FOU_50 = data[12]
            p_FOU_50 = data[4]/1e3
            t_FOU_50 = data[3]
            n=50
            x_50 = np.linspace(50/(2*n),50*(1-1/(2*n)),n)

        '''-------------------------------MUSCL------------------------------'''

        datas = np.load('flying/results_6k_SchmallNEW_50_IMPEC_MUSCL_LLF_VA_E_479.npy', allow_pickle=True)
        datas = np.load('flying/results_6k_SchmallNEW_50_IMPEC_MUSCL_MDW_VA_E_409.npy', allow_pickle=True)
        datas = np.load('flying/results_6k_SchmallNEW_50_IMPEC_MUSCL_MDW_VL_E_401.npy', allow_pickle=True)
        datas = np.load('flying/results_6k_SchmallNEW_50_IMPEC_MUSCL_MDW_VL_RK3_399.npy', allow_pickle=True)
        #datas = np.load('flying/results_6k_SchmallNEW_50_IMPEC_MUSCL_LLF_VL_E_477.npy', allow_pickle=True)
        #datas = np.load('flying/results_6k_SchmallNEW_50_IMPEC_MUSCL_LLF_VL_RK3_472.npy', allow_pickle=True)
        datas = np.load('flying/results_6k_SchmallNEW_50_IMPEC_MUSCL_ROE_VL_RK3_494.npy', allow_pickle=True)
        for data in datas[1:]:
            Sw_MUSCL_50 = data[5]
            So_MUSCL_50 = data[6]
            Sg_MUSCL_50 = data[7]
            z_MUSCL_50 = data[10]
            Nk_MUSCL_50 = data[12]
            p_MUSCL_50 = data[4]/1e3
            t_MUSCL_50 = data[3]
            n=50
            x_50 = np.linspace(50/(2*n),50*(1-1/(2*n)),n)

        plt.figure(1)
        plt.title('t = 200 days - 50x1x1 mesh')
        plt.plot(x_CMG, So_CMG, 'k')
        plt.plot(x_50, So_FOU_50, '-b', mfc='none')
        plt.plot(x_50, So_MUSCL_50, '-m', mfc='none')
        #plt.plot(x_16, z_FR2_16[0,:], '-gs', mfc='none')
        #plt.plot(x_16, z_FR3_16[0,:], '-y<', mfc='none')
        #plt.plot(x_16, z_FR4_16[0,:], '-r*', mfc='none')
        plt.ylabel('$S_o$')
        plt.xlabel('Distance')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.savefig('results/compositional/SPE_paper/6k_So_50.png')
        plt.close()

        plt.figure(1)
        plt.title('t = 200 days - 50x1x1 mesh')
        plt.plot(x_CMG, Sw_CMG, 'k')
        plt.plot(x_50, Sw_FOU_50, '-b', mfc='none')
        plt.plot(x_50, Sw_MUSCL_50, '-m', mfc='none')
        #plt.plot(x_16, z_FR2_16[0,:], '-gs', mfc='none')
        #plt.plot(x_16, z_FR3_16[0,:], '-y<', mfc='none')
        #plt.plot(x_16, z_FR4_16[0,:], '-r*', mfc='none')
        plt.ylabel('$S_w$')
        plt.xlabel('Distance')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.savefig('results/compositional/SPE_paper/6k_Sw_50.png')
        plt.close()

        plt.figure(1)
        plt.title('t = 200 days - 50x1x1 mesh')
        plt.plot(x_CMG, zC1_CMG, 'k')
        plt.plot(x_50, z_FOU_50[0,:], '-bo', mfc='none')
        plt.plot(x_50, z_MUSCL_50[0,:], '-mv', mfc='none')
        #plt.plot(x_16, z_FR2_16[0,:], '-gs', mfc='none')
        #plt.plot(x_16, z_FR3_16[0,:], '-y<', mfc='none')
        #plt.plot(x_16, z_FR4_16[0,:], '-r*', mfc='none')
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance')
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.grid()
        plt.savefig('results/compositional/SPE_paper/6k_zC1_50.png')
        plt.close()


        import pdb; pdb.set_trace()
