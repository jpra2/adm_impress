import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
x_CMG = np.loadtxt('xCH4_case2MM_CMG.txt')
xCH4_CMG = np.loadtxt('xxCH4_case2MM_CMG.txt')
for  arq in arquivos:
    if  arq.startswith(name):
        '''-------------------------MUSCL LLF RESULTS------------------------'''
        '''x_axis_xC30 = np.linspace(0, 21.4, 100)
        xC30 = np.zeros(len(x_axis_xC30))
        x_axis_xC3 = np.array([21.5484, 21.7646, 21.9171, 21.968, 22.1459, 22.5655, 24.8665, 26.2904,
            26.5446])
        xC3 = np.array([2.867384e-4, 0.00258065, 0.00659498, 0.0219355, 0.0840143, 0.317563, 0.308817, 0.303656,
            0.])'''

        x_axis_xCH4 = np.array([0, 23.1051, 23.1051, 23.1476, 23.9119, 24.0817, 24.4958, 24.9522,  25.7909,
                        26.8737, 27.7495, 27.7866, 27.8822, 27.9989, 28.0945, 28.259,  28.4289, 28.5881,
                        28.7686, 28.9756, 29.241, 29.6125, 29.8779, 30.483,50])
        xCH4 = np.array([0, 1.743896e-4, 0.196886, 0.332735, 0.332735, 0.331166, 0.329248, 0.325237, 0.321749, 0.318784, 0.318261,
                        0.307623, 0.239437, 0.179098, 0.139686, 0.0953911, 0.0636522, 0.0430742, 0.0273792,
                        0.0162182, 0.00767314, 0.00348779, 0.00156951, 0.0,0])

        #x_axis_xC31 = np.linspace(x_axis_xC3[-1], 50, 100)
        #xC31 = np.zeros(len(x_axis_xC31))

        '''x_axis_xC3 = np.concatenate((x_axis_xC30,x_axis_xC3))
        x_axis_xC3 = np.concatenate((x_axis_xC3,x_axis_xC31))
        xC3 = np.concatenate((xC30,xC3))
        xC3 = np.concatenate((xC3,xC31))'''


        #f = interp1d(x_axis_xCH4,xCH4)

        """---------------------- Convergence Study -------------------------"""

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_5000_FOU_UPW_b0_wa_5487.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_5000 = data[6]
            Sg_FOU_5000 = data[7]
            xkj_FOU_5000 = data[13]
            z_FOU_5000 = data[10]
            P_FOU_5000 = data[4]
            xkj_FOU_5000[:,1,Sg_FOU_5000==0] = 0
            xkj_FOU_5000[:,0,So_FOU_5000==0] = 0
            n = 5000
            x_5000 = np.linspace(0+50/(2*n),50-50/(2*n),n)
            t_FOU_5000 = data[2]
            f = interp1d(x_5000,xkj_FOU_5000[0,0])

        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FOU_LLF_b0_u_406.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FOU_LLF_b0_ha_502.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FOU_UPW_b0_ha_420.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FOU_MDW_b0_ha_417.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FOU_MDW_b0_hao_424.npy', allow_pickle=True)
        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FOU_UPW_b0_wa_415.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FOU_UPW_b0_u_366.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FOU_MDW_b0_u_427.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FOU_LLF_b0_wa_419.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FOU_200 = data[6]
            Sg_FOU_200 = data[7]
            xkj_FOU_200 = data[13]
            z_FOU_200 = data[10]
            P_FOU_200 = data[4]
            xkj_FOU_200[:,1,Sg_FOU_200==0] = 0
            xkj_FOU_200[:,0,So_FOU_200==0] = 0
            n = 200
            x_200 = np.linspace(0+50/(2*n),50-50/(2*n),n)
            t_FOU_200 = data[2]
            e200_L1_FOU = (sum(abs(f(x_200)-xkj_FOU_200[0,0,:]))*(1/n))


        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_MUSCL_MDW_mmod_b0_ha_894.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_MUSCL_MDW_mmod_b0_hao_907.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_MUSCL_MDW_mmod_b0_wa_892.npy', allow_pickle=True)
        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_MUSCL_MDW_RK3_b0_wa_953.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_MUSCL_LLF_RK3_b0_wa_1173.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_MUSCL_MDW_RK3_VA_b0_wa_885.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_200 = data[6]
            Sg_MUSCL_200 = data[7]
            xkj_MUSCL_200 = data[13]
            z_MUSCL_200 = data[10]
            P_MUSCL_200 = data[4]
            n = 200
            xkj_MUSCL_200[:,1,Sg_MUSCL_200==0] = 0
            xkj_MUSCL_200[:,0,So_MUSCL_200==0] = 0
            e200_L1_MUSCL = (sum(abs(f(x_200)-xkj_MUSCL_200[0,0,:]))*(1/n))


        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FR2_MDW_b0_wa_914.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FR2_LLF_b0_wa_1085.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR2_200 = data[6]
            Sg_FR2_200 = data[7]
            xkj_FR2_200 = data[13]
            z_FR2_200 = data[10]
            P_FR2_200 = data[4]
            n = 200
            xkj_FR2_200[:,1,Sg_FR2_200==0] = 0
            xkj_FR2_200[:,0,So_FR2_200==0] = 0
            e200_L1_FR2 = (sum(abs(f(x_200)-xkj_FR2_200[0,0,:]))*(1/n))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FR3_MDW_b0_wa_1350.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FR3a_MDW_b0_wa_1350.npy', allow_pickle=True)
        #datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FR3_MDW_b0_u_1418.npy', allow_pickle=True)
        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FR3_MDW_b0_wa_1852.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR3_200 = data[6]
            Sg_FR3_200 = data[7]
            xkj_FR3_200 = data[13]
            z_FR3_200 = data[10]
            P_FR3_200 = data[4]
            n = 200
            xkj_FR3_200[:,1,Sg_FR3_200==0] = 0
            xkj_FR3_200[:,0,So_FR3_200==0] = 0
            e200_L1_FR3 = (sum(abs(f(x_200)-xkj_FR3_200[0,0,:]))*(1/n))

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_100_FR4_MDW_wa_978.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR4_100 = data[6]
            Sg_FR4_100 = data[7]
            xkj_FR4_100 = data[13]
            n = 100
            xkj_FR4_100[:,1,Sg_FR4_100==0] = 0
            xkj_FR4_100[:,0,So_FR4_100==0] = 0

        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_200_FR4_MDW_b0_wa_1845.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_FR4_200 = data[6]
            Sg_FR4_200 = data[7]
            xkj_FR4_200 = data[13]
            z_FR4_200 = data[10]
            P_FR4_200 = data[4]
            n = 200
            xkj_FR4_200[:,1,Sg_FR4_200==0] = 0
            xkj_FR4_200[:,0,So_FR4_200==0] = 0
            e200_L1_FR4 = (sum(abs(f(x_200)-xkj_FR4_200[0,0,:]))*(1/n))



        datas = np.load('flying/results_case2_Moshiri_Manzari_3k_100_MUSCL_MDW_mmod_hao_478.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            So_MUSCL_100 = data[6]
            Sg_MUSCL_100 = data[7]
            xkj_MUSCL_100 = data[13]
            z_MUSCL_100 = data[10]
            xkj_MUSCL_100[:,1,Sg_MUSCL_100==0] = 0
            xkj_MUSCL_100[:,0,So_MUSCL_100==0] = 0
            n = 100
            x_100 = np.linspace(0+50/(2*n),50-50/(2*n),n)


        plt.figure(1)
        plt.plot(x_200, xkj_FOU_200[0,0], '-bo', mfc='none')
        plt.plot(x_200, xkj_MUSCL_200[0,0], '-m<', mfc='none')
        plt.plot(x_200, xkj_FR2_200[0,0], '-gs', mfc='none')
        plt.plot(x_200, xkj_FR3_200[0,0], '-yv', mfc='none')
        plt.plot(x_200, xkj_FR4_200[0,0], '-r*', mfc='none')
        #plt.plot(x_CMG, xCH4_CMG, 'k')
        #plt.plot(x_axis_xCH4, xCH4, '-k')
        plt.plot(x_5000, xkj_FOU_5000[0,0], '-k')

        plt.xlim(20,35)
        plt.ylim(-0.05, 0.4)
        plt.grid()
        plt.legend(( 'FOU (200 CV)', 'MUSCL (200 CV)', 'FR-P1 (200 CV)', 'FR-P2 (200 CV)', 'FR-P3 (200 CV)', 'FOU (5000 CV)'))
        plt.title('Results after 186 days')
        plt.ylabel('$x_{C1}$')
        plt.xlabel('Distance[m]')
        plt.savefig('results/compositional/SPE_paper/3k_methane_x_MM_FOU_200.png')
        plt.close()

        plt.figure(1)
        plt.plot(x_200, z_FOU_200[0,:], '-bo', mfc='none')
        plt.plot(x_200, z_MUSCL_200[0,:], '-m<', mfc='none')
        plt.plot(x_200, z_FR2_200[0,:], '-gs', mfc='none')
        plt.plot(x_200, z_FR3_200[0,:], '-yv', mfc='none')
        plt.plot(x_200, z_FR4_200[0,:], '-r*', mfc='none')
        #plt.plot(x_CMG, xCH4_CMG, 'k')
        #plt.plot(x_axis_xCH4, xCH4, '-k')
        plt.plot(x_5000, z_FOU_5000[0,:], '-k')
        plt.grid()
        plt.xlim(20,35)
        plt.legend(( 'FOU', 'MUSCL-200', 'FR-P1', 'FR-P2', 'FR-P3', 'FOU-5000'))
        plt.title('Results for 200x1x1 mesh')
        plt.ylabel('$z_{C1}$')
        plt.xlabel('Distance[m]')
        plt.savefig('results/compositional/SPE_paper/3k_methane_z_MM_FOU_200.png')
        plt.close()

        plt.figure(1)
        plt.plot(x_200, z_FOU_200[1,:], '-bo', mfc='none')
        plt.plot(x_200, z_MUSCL_200[1,:], '-m<', mfc='none')
        plt.plot(x_200, z_FR2_200[1,:], '-gs', mfc='none')
        plt.plot(x_200, z_FR3_200[1,:], '-yv', mfc='none')
        plt.plot(x_200, z_FR4_200[1,:], '-r*', mfc='none')
        #plt.plot(x_CMG, xCH4_CMG, 'k')
        #plt.plot(x_axis_xCH4, xCH4, '-k')
        plt.plot(x_5000, z_FOU_5000[1,:], '-k')
        plt.grid()
        plt.xlim(20,35)
        plt.legend(( 'FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3', 'FOU-5000'))
        plt.title('Results for 200x1x1 mesh')
        plt.ylabel('$z_{C2}$')
        plt.xlabel('Distance[m]')
        plt.savefig('results/compositional/SPE_paper/3k_ethane_z_MM_FOU_200.png')
        plt.close()

        plt.figure(1)
        plt.plot(x_200, P_FOU_200, '-bo', mfc='none')
        plt.plot(x_200, P_MUSCL_200, '-m<', mfc='none')
        plt.plot(x_200, P_FR2_200, '-gs', mfc='none')
        plt.plot(x_200, P_FR3_200, '-yv', mfc='none')
        plt.plot(x_200, P_FR4_200, '-r*', mfc='none')
        #plt.plot(x_CMG, xCH4_CMG, 'k')
        #plt.plot(x_axis_xCH4, xCH4, '-k')
        plt.plot(x_5000, P_FOU_5000, '-k')
        plt.grid()
        #plt.xlim(20,35)
        plt.legend(( 'FOU', 'MUSCL-200', 'FR-P1', 'FR-P2', 'FR-P3', 'FOU-5000'))
        plt.title('Results for 200x1x1 mesh')
        plt.ylabel('$p$')
        plt.xlabel('Distance[m]')
        plt.savefig('results/compositional/SPE_paper/3k_p_MM_FOU_200.png')
        plt.close()
        import pdb; pdb.set_trace()
