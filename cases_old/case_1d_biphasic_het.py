import numpy as np
import matplotlib.pyplot as plt #Biblioteca para plotar resultados
from scipy.interpolate import interp1d
import os


flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
L=6.096

i1 = 1
i2 = 2

for  arq in arquivos:
    if  arq.startswith(name):
        fx = open('x_biphasic_het_CMG','r')
        x_CMG = [float(line.rstrip('\n\r')) for line in fx]
        x_CMG.append(6.096)

        fSw = open('Sw_biphasic_het_CMG','r')
        Sw_CMG = [float(line.rstrip('\n\r')) for line in fSw]

        fSw10sec = open('Sw_biphasic_het_10sec_CMG','r')
        Sw_10sec_CMG = [float(line.rstrip('\n\r')) for line in fSw10sec]
        Sw_10sec_CMG.append(0.436)
        fSw_CMG = interp1d(x_CMG,Sw_10sec_CMG)

        """  FOU  """

        datas = np.load('flying/results_biphasic_het_10sec_1000CV_IMPEC_1294.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_1000 = data[5]
            n=1000
            x1000 = np.linspace(0+6.096/(2*n),6.096*(1-1/(2*n)),n)
            #x1000 = np.linspace(0,6.096*(1-2/(2*n)),n)
            #fSw_CMG = interp1d(x1000,Sw_10s_1000)
            #eL1_FOU_1000CV = (sum(abs(fSw_CMG(x100)-Sw_10s_100))*(L/100))

        datas = np.load('flying/results_biphasic_het_10sec_500CV_IMPEC_648.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_500 = data[5]
            x500 = np.linspace(0+6.096/1000,6.096*(1-1/1000),500)
            #x500[-1] = x_CMG[-1]
            #x500 = np.linspace(0,6.096*(1-2/1000),500)
            #import pdb; pdb.set_trace()
            Sw0_10sec_vals = np.array([0.48227077, 0.62165371, 0.28547398, 0.44326456, \
                                0.31831652, 0.63923074, 0.377705, 0.48091657, \
                                0.54459622, 0.23253678])
            eL1_FOU_500CV = (sum(abs(fSw_CMG(x500)-Sw_10s_500))*(L/500))

        datas = np.load('flying/results_biphasic_het_10sec_50CV_IMPEC_106.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_50 = data[5]
            n=50
            x50 = np.linspace(0+6.096/(2*n),6.096*(1-1/(2*n)),n)
            #x50 = np.linspace(0,6.096*(1-2/(2*n)),n)
            eL1_FOU_50CV = (sum(abs(fSw_CMG(x50)-Sw_10s_50))*(L/50))

        datas = np.load('flying/results_biphasic_het_10sec_100CV_IMPEC_133.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_100 = data[5]
            n=100
            x100 = np.linspace(0+6.096/(2*n),6.096*(1-1/(2*n)),n)
            #x100 = np.linspace(0,6.096*(1-2/(2*n)),n)
            eL1_FOU_100CV = (sum(abs(fSw_CMG(x100)-Sw_10s_100))*(L/n))

        datas = np.load('flying/results_biphasic_het_10sec_200CV_IMPEC_261.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_200 = data[5]
            n=200
            x200 = np.linspace(0+6.096/(2*n),6.096*(1-1/(2*n)),n)
            #x100 = np.linspace(0,6.096*(1-2/(2*n)),n)
            t_FOU_200 = data[3]
            eL1_FOU_200CV = (sum(abs(fSw_CMG(x200)-Sw_10s_200))*(L/n))

        datas = np.load('flying/results_biphasic_het_10sec_300CV_IMPEC_390.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_300 = data[5]
            n=300
            t_FOU_300 = data[2]
            x300 = np.linspace(0+6.096/(2*n),6.096*(1-1/(2*n)),n)
            #x100 = np.linspace(0,6.096*(1-2/(2*n)),n)
            #x300[-1] = x_CMG[300]
            eL1_FOU_300CV = (sum(abs(fSw_CMG(x300)-Sw_10s_300))*(L/n))

        n=150
        x150 = np.linspace(0+6.096/(2*n),6.096*(1-1/(2*n)),n)
        n=170
        x170 = np.linspace(0+6.096/(2*n),6.096*(1-1/(2*n)),n)

        """  MUSCL  """

        datas = np.load('flying/results_biphasic_het_10sec_50CV_IMPEC_MUSCL_LLF_VA_RK3_198.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_MUSCL_50 = data[5]
            t_MUSCL_50CV = data[2]
            eL1_MUSCL_50CV = (sum(abs(fSw_CMG(x50)-Sw_10s_MUSCL_50))*(L/50))

        datas = np.load('flying/results_biphasic_het_10sec_100CV_IMPEC_MUSCL_LLF_VA_RK3_393.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_MUSCL_100 = data[5]
            eL1_MUSCL_100CV = (sum(abs(fSw_CMG(x100)-Sw_10s_MUSCL_100))*(L/100))

        datas = np.load('flying/results_biphasic_het_10sec_170CV_IMPEC_MUSCL_LLF_VA_RK3_664.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_MUSCL_170 = data[5]
            t_MUSCL_170CV = data[2]
            eL1_MUSCL_170CV = (sum(abs(fSw_CMG(x170)-Sw_10s_MUSCL_170))*(L/170))

        datas = np.load('flying/results_biphasic_het_10sec_500CV_IMPEC_MUSCL_LLF_VL_RK3_1941.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_MUSCL_500 = data[5]
            eL1_MUSCL_500CV = (sum(abs(fSw_CMG(x500)-Sw_10s_MUSCL_500))*(L/500))

        """  FR P1  """
        datas = np.load('flying/results_biphasic_het_10sec_50CV_IMPEC_FRP1_LLF_RK3_200.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_FRP1_50 = data[5]
            eL1_FRP1_50CV = (sum(abs(fSw_CMG(x50)-Sw_10s_FRP1_50))*(L/50))

        datas = np.load('flying/results_biphasic_het_10sec_100CV_IMPEC_FRP1_LLF_RK3_396.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_FRP1_100 = data[5]
            eL1_FRP1_100CV = (sum(abs(fSw_CMG(x100)-Sw_10s_FRP1_100))*(L/100))

        datas = np.load('flying/results_biphasic_het_10sec_500CV_IMPEC_FRP1_LLF_RK3_1972.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_FRP1_500 = data[5]
            eL1_FRP1_500CV = (sum(abs(fSw_CMG(x500)-Sw_10s_FRP1_500))*(L/500))

        """  FR P2  """
        datas = np.load('flying/results_biphasic_het_10sec_50CV_IMPEC_FRP2_LLF_RK3_329.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_FRP2_50 = data[5]
            eL1_FRP2_50CV = (sum(abs(fSw_CMG(x50)-Sw_10s_FRP2_50))*(L/50))

        datas = np.load('flying/results_biphasic_het_10sec_100CV_IMPEC_FRP2_LLF_RK3_655.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_FRP2_100 = data[5]
            eL1_FRP2_100CV = (sum(abs(fSw_CMG(x100)-Sw_10s_FRP2_100))*(L/100))

        datas = np.load('flying/results_biphasic_het_10sec_500CV_IMPEC_FRP2_LLF_RK3_3268.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_FRP2_500 = data[5]
            eL1_FRP2_500CV = (sum(abs(fSw_CMG(x500)-Sw_10s_FRP2_500))*(L/500))

        """  FR P3  """
        datas = np.load('flying/results_biphasic_het_10sec_50CV_IMPEC_FRP3_LLF_RK3_466.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_FRP3_50 = data[5]
            eL1_FRP3_50CV = (sum(abs(fSw_CMG(x50)-Sw_10s_FRP3_50))*(L/50))

        datas = np.load('flying/results_biphasic_het_10sec_100CV_IMPEC_FRP3_LLF_RK3_930.npy', allow_pickle=True)
        datas = np.load('flying/results_biphasic_het_10sec_100CV_IMPEC_FRP3_LLF_RK3_930.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_FRP3_100 = data[5]
            t_FRP3_100CV = data[2]
            eL1_FRP3_100CV = (sum(abs(fSw_CMG(x100)-Sw_10s_FRP3_100))*(L/100))

        datas = np.load('flying/results_biphasic_het_10sec_500CV_IMPEC_FRP3_LLF_RK3_4634.npy', allow_pickle=True)
        for data in datas[i1:i2]:
            Sw_10s_FRP3_500 = data[5]
            eL1_FRP3_500CV = (sum(abs(fSw_CMG(x500)-Sw_10s_FRP3_500))*(L/500))


        Sw0_10sec = np.ones_like(Sw_10s_500)
        Sw0_10sec[0:50] = Sw0_10sec_vals[0]
        Sw0_10sec[50:100] = Sw0_10sec_vals[1]
        Sw0_10sec[100:150] = Sw0_10sec_vals[2]
        Sw0_10sec[150:200] = Sw0_10sec_vals[3]
        Sw0_10sec[200:250] = Sw0_10sec_vals[4]
        Sw0_10sec[250:300] = Sw0_10sec_vals[5]
        Sw0_10sec[300:350] = Sw0_10sec_vals[6]
        Sw0_10sec[350:400] = Sw0_10sec_vals[7]
        Sw0_10sec[400:450] = Sw0_10sec_vals[8]
        Sw0_10sec[450:500] = Sw0_10sec_vals[9]


        plt.figure(1)
        #plt.plot(x_CMG, Sw_10sec_CMG, 'k')
        plt.plot(x1000, Sw_10s_1000, 'k')
        plt.plot(x500, Sw_10s_500, '-b', mfc='none')
        plt.plot(x500, Sw_10s_MUSCL_500, '-m', mfc='none')
        plt.plot(x500, Sw_10s_FRP1_500, '-g', mfc='none')
        plt.plot(x500, Sw_10s_FRP2_500, '--y', mfc='none')
        plt.plot(x500, Sw_10s_FRP3_500, '-r', mfc='none')
        plt.grid()
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P3', 'FR-P3'))
        plt.ylabel('Sw')
        plt.xlabel('x')
        plt.savefig('results/compositional/Sw_biphasic_het_10sec.png')
        plt.close()

        plt.figure(3)
        #plt.plot(x_CMG, Sw_10sec_CMG, 'k')
        plt.plot(x1000, Sw_10s_1000, 'k')
        plt.plot(x100, Sw_10s_100, '-b', mfc='none')
        plt.plot(x100, Sw_10s_MUSCL_100, '-m', mfc='none')
        plt.plot(x100, Sw_10s_FRP1_100, '-g', mfc='none')
        plt.plot(x100, Sw_10s_FRP2_100, '--y', mfc='none')
        plt.plot(x100, Sw_10s_FRP3_100, '-r', mfc='none')
        plt.grid()
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P3', 'FR-P3'))
        plt.ylabel('Sw')
        plt.xlabel('x')
        plt.savefig('results/compositional/Sw_biphasic_het_10sec_100.png')
        plt.close()

        plt.figure(2)
        #plt.plot(x_CMG, Sw_10sec_CMG, '--k')
        plt.plot(x1000, Sw_10s_1000, 'k')
        plt.plot(x50, Sw_10s_50, '-b', mfc='none')
        plt.plot(x50, Sw_10s_MUSCL_50, '-m', mfc='none')
        plt.plot(x50, Sw_10s_FRP1_50, '-g', mfc='none')
        plt.plot(x50, Sw_10s_FRP2_50, '--y', mfc='none')
        plt.plot(x100, Sw_10s_FRP3_100, '-r', mfc='none')
        plt.grid()
        plt.legend(('Reference', 'FOU', 'MUSCL', 'FR-P1', 'FR-P3', 'FR-P3'))
        plt.ylabel('Sw')
        plt.xlabel('x')
        plt.savefig('results/compositional/Sw_biphasic_het_10sec_50CV.png')
        plt.close()

        plt.figure(2)
        plt.plot(x_CMG, Sw_10sec_CMG, '-k')
        #plt.plot(x1000, Sw_10s_1000, 'k')
        plt.plot(x300, Sw_10s_300, '--b', mfc='none')
        plt.plot(x170, Sw_10s_MUSCL_170, '-.m', mfc='none')
        plt.plot(x100, Sw_10s_FRP3_100, '.r', mfc='none', markersize=5)
        plt.grid()
        plt.legend(('Reference', 'FOU (300 CV)', 'MUSCL (170 CV)', 'FR-P3 (100 CV)'))
        plt.ylabel('Sw')
        plt.xlabel('Distance [m]')
        plt.savefig('results/compositional/Sw_biphasic_het_10sec_otim_comp.png')
        plt.close()

        plt.figure(2)
        plt.plot(x_CMG, Sw_10sec_CMG, '-k')
        #plt.plot(x1000, Sw_10s_1000, 'k')
        plt.plot(x100, Sw_10s_100, '--b', mfc='none')
        plt.plot(x100, Sw_10s_MUSCL_100, '-.m', mfc='none')
        plt.plot(x100, Sw_10s_FRP3_100, '.r', mfc='none', markersize=5)
        plt.grid()
        plt.legend(('Reference', 'FOU (100 CV)', 'MUSCL (100 CV)', 'FR-P3 (100 CV)'))
        plt.ylabel('Sw')
        plt.xlabel('Distance [m]')
        plt.savefig('results/compositional/Sw_biphasic_het_10sec_otim_comp_100CV.png')
        plt.close()

        plt.figure(2)
        plt.plot(x_CMG, Sw_10sec_CMG, '-k')
        #plt.plot(x1000, Sw_10s_1000, 'k')
        plt.plot(x50, Sw_10s_50, '--b', mfc='none')
        plt.plot(x50, Sw_10s_MUSCL_50, '-.m', mfc='none')
        plt.plot(x50, Sw_10s_FRP3_50, '.r', mfc='none', markersize=5)
        plt.grid()
        plt.legend(('Reference', 'FOU (50 CV)', 'MUSCL (50 CV)', 'FR-P3 (50 CV)'))
        plt.ylabel('Sw')
        plt.xlabel('Distance [m]')
        plt.savefig('results/compositional/Sw_biphasic_het_10sec_otim_comp_50CV.png')
        plt.close()

        plt.figure(2)
        plt.plot(x500, Sw0_10sec, 'k')
        plt.grid()
        plt.ylabel('$S_w^0$')
        plt.xlabel('Distance [m]')
        plt.savefig('results/compositional/Sw0_biphasic_het_10sec.png')
        plt.close()



        import pdb; pdb.set_trace()
