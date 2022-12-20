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
n_200 = 200
n_100 = 100
n_50 = 50
n_25 = 25

"""
fx = open('x_points_CMG_SchmallNEW.txt','r')
x_CMG = np.loadtxt('x_points_CMG_SchmallNEW.txt')

#fP = open('P1024.txt','r')
#P_CMG = [float(line.rstrip('\n\r')) for line in fP]

fSw = open('Sw_points_CMG_SchmallNEW.txt','r')
Sw_CMG = [float(line.rstrip('\n\r')) for line in fSw]

fSo = open('So_points_CMG_SchmallNEW.txt','r')
So_CMG = [float(line.rstrip('\n\r')) for line in fSo]

fSg = open('Sg_points_CMG_SchmallNEW.txt','r')
Sg_CMG = [float(line.rstrip('\n\r')) for line in fSg]
"""


x_CMG = np.loadtxt('x_points_CMG_SchmallNEW.txt')
x_CMG[-1] = 50.0
Sw_CMG = np.loadtxt('Sw_points_CMG_SchmallNEW.txt')
So_CMG = np.loadtxt('So_points_CMG_SchmallNEW.txt')
Sg_CMG = np.loadtxt('Sg_points_CMG_SchmallNEW.txt')

f_Sw = interp1d(x_CMG,Sw_CMG)
f_So = interp1d(x_CMG,So_CMG)
f_Sg = interp1d(x_CMG,Sg_CMG)


for arq in arquivos:
    if  arq.startswith(name):


        datas3 = np.load('flying/6k/25_0_0/results_6k_SchmallNEW_25_FI_20days_POCONEW_40.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI_25 = data[5]
            So_FI_25 = data[6]
            Sg_FI_25 = data[7]
            x1_25 = np.linspace(0, 50, n_25)

            e25_Sw_L1_FI = (sum(abs(f_Sw(x1_25)-Sw_FI_25))*(1/25))
            e25_Sw_L2_FI = np.sqrt(np.sum((f_Sw(x1_25)-Sw_FI_25)**2) * 1 / 25)

            e25_So_L1_FI = (sum(abs(f_So(x1_25)-So_FI_25))*(1/25))
            e25_So_L2_FI = np.sqrt(np.sum((f_So(x1_25)-So_FI_25)**2) * 1 / 25)

            e25_Sg_L1_FI = (sum(abs(f_Sg(x1_25)-Sg_FI_25))*(1/25))
            e25_Sg_L2_FI = np.sqrt(np.sum((f_Sg(x1_25)-Sg_FI_25)**2) * 1 / 25)


        datas3 = np.load('flying/6k/results_6k_SchmallNEW_50_FI_20days_2DTeste_203.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI_50 = data[5]
            So_FI_50 = data[6]
            Sg_FI_50 = data[7]
            x1_50 = np.linspace(0, 50, n_50)

            e50_Sw_L1_FI = (sum(abs(f_Sw(x1_50)-Sw_FI_50))*(1/50))
            R50_Sw_L1_FI = math.log(e25_Sw_L1_FI/e50_Sw_L1_FI,2)
            e50_Sw_L2_FI = np.sqrt(np.sum((f_Sw(x1_50)-Sw_FI_50)**2) * 1 / 50)
            R50_Sw_L2_FI = math.log(e25_Sw_L2_FI/e50_Sw_L2_FI,2)

            e50_So_L1_FI = (sum(abs(f_So(x1_50)-So_FI_50))*(1/50))
            R50_So_L1_FI = math.log(e25_So_L1_FI/e50_So_L1_FI,2)
            e50_So_L2_FI = np.sqrt(np.sum((f_So(x1_50)-So_FI_50)**2) * 1 / 50)
            R50_So_L2_FI = math.log(e25_So_L2_FI/e50_So_L2_FI,2)

            e50_Sg_L1_FI = (sum(abs(f_Sg(x1_50)-Sg_FI_50))*(1/50))
            R50_Sg_L1_FI = math.log(e25_Sg_L1_FI/e50_Sg_L1_FI,2)
            e50_Sg_L2_FI = np.sqrt(np.sum((f_Sg(x1_50)-Sg_FI_50)**2) * 1 / 50)
            R50_Sg_L2_FI = math.log(e25_Sg_L2_FI/e50_Sg_L2_FI,2)

        datas2 = np.load('flying/6k/100_0_0/results_6k_SchmallNEW_100_FI_20days_POCONEW_dtmenor_2003.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas2[1:]:
            Sw_FI_100 = data[5]
            So_FI_100 = data[6]
            Sg_FI_100 = data[7]
            x1_100 = np.linspace(0, 50, n_100)

            e100_Sw_L1_FI = (sum(abs(f_Sw(x1_100)-Sw_FI_100))*(1/100))
            R100_Sw_L1_FI = math.log(e50_Sw_L1_FI/e100_Sw_L1_FI,2)
            e100_Sw_L2_FI = np.sqrt(np.sum((f_Sw(x1_100)-Sw_FI_100)**2) * 1 / 100)
            R100_Sw_L2_FI = math.log(e50_Sw_L2_FI/e100_Sw_L2_FI,2)

            e100_So_L1_FI = (sum(abs(f_So(x1_100)-So_FI_100))*(1/100))
            R100_So_L1_FI = math.log(e50_So_L1_FI/e100_So_L1_FI,2)
            e100_So_L2_FI = np.sqrt(np.sum((f_So(x1_100)-So_FI_100)**2) * 1 / 100)
            R100_So_L2_FI = math.log(e50_So_L2_FI/e100_So_L2_FI,2)

            e100_Sg_L1_FI = (sum(abs(f_Sg(x1_100)-Sg_FI_100))*(1/100))
            R100_Sg_L1_FI = math.log(e50_Sg_L1_FI/e100_Sg_L1_FI,2)
            e100_Sg_L2_FI = np.sqrt(np.sum((f_Sg(x1_100)-Sg_FI_100)**2) * 1 / 100)
            R100_Sg_L2_FI = math.log(e50_Sg_L2_FI/e100_Sg_L2_FI,2)

        datas3 = np.load('flying/6k/200_0_0/results_6k_SchmallNEW_200_FI_20days_POCONEW_2003.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI_200 = data[5]
            So_FI_200 = data[6]
            Sg_FI_200 = data[7]
            x1_200 = np.linspace(0, 50, n_200)

            e200_Sw_L1_FI = (sum(abs(f_Sw(x1_200)-Sw_FI_200))*(1/200))
            R200_Sw_L1_FI = math.log(e100_Sw_L1_FI/e200_Sw_L1_FI,2)
            e200_Sw_L2_FI = np.sqrt(np.sum((f_Sw(x1_200)-Sw_FI_200)**2) * 1 / 200)
            R200_Sw_L2_FI = math.log(e100_Sw_L2_FI/e200_Sw_L2_FI,2)

            e200_So_L1_FI = (sum(abs(f_So(x1_200)-So_FI_200))*(1/200))
            R200_So_L1_FI = math.log(e100_So_L1_FI/e200_So_L1_FI,2)
            e200_So_L2_FI = np.sqrt(np.sum((f_So(x1_200)-So_FI_200)**2) * 1 / 200)
            R200_So_L2_FI = math.log(e100_So_L2_FI/e200_So_L2_FI,2)

            e200_Sg_L1_FI = (sum(abs(f_Sg(x1_200)-Sg_FI_200))*(1/200))
            R200_Sg_L1_FI = math.log(e100_Sg_L1_FI/e200_Sg_L1_FI,2)
            e200_Sg_L2_FI = np.sqrt(np.sum((f_Sg(x1_200)-Sg_FI_200)**2) * 1 / 200)
            R200_Sg_L2_FI = math.log(e100_Sg_L2_FI/e200_Sg_L2_FI,2)



        plt.figure(1)
        x = np.log2(np.array([25, 50, 100, 200]))
        Sw_FI = np.log2(np.array([e25_Sw_L1_FI, e50_Sw_L1_FI, e100_Sw_L1_FI, e200_Sw_L1_FI]))

        y_ref = -x-2.0

        plt.plot(x, Sw_FI, '-mP', mfc='none')
        plt.plot(x, y_ref, 'k')
        #plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{2}({E}_{L_1} Sw)$')
        plt.xlabel('$log_{2}(N)$')
        plt.legend(('FI', 'Convergência linear'))
        plt.grid()
        plt.savefig('results/6k_L1_convergence_FI_Sw.png')



        plt.figure(2)
        x = np.log2(np.array([25, 50, 100, 200]))
        So_FI = np.log2(np.array([e25_So_L1_FI, e50_So_L1_FI, e100_So_L1_FI, e200_So_L1_FI]))

        y_ref = -x-2.0

        plt.plot(x, So_FI, '-mP', mfc='none')
        plt.plot(x, y_ref, 'k')
        #plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{2}({E}_{L_1} So)$')
        plt.xlabel('$log_{2}(N)$')
        plt.legend(('FI', 'Convergência linear'))
        plt.grid()
        plt.savefig('results/6k_L1_convergence_FI_So.png')


        plt.figure(3)
        x = np.log2(np.array([25, 50, 100, 200]))
        Sg_FI = np.log2(np.array([e25_Sg_L1_FI, e50_Sg_L1_FI, e100_Sg_L1_FI, e200_Sg_L1_FI]))

        y_ref = -x-2.0

        plt.plot(x, Sg_FI, '-mP', mfc='none')
        plt.plot(x, y_ref, 'k')
        #plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{2}({E}_{L_1} Sg)$')
        plt.xlabel('$log_{2}(N)$')
        plt.legend(('FI', 'Convergência linear'))
        plt.grid()
        plt.savefig('results/6k_L1_convergence_FI_Sg.png')



        plt.figure(4)
        x = np.log2(np.array([25, 50, 100, 200]))
        Sw_FI_L2 = np.log2(np.array([e25_Sw_L2_FI, e50_Sw_L2_FI, e100_Sw_L2_FI, e200_Sw_L2_FI]))

        plt.plot(x, Sw_FI_L2,'-bo')

        #ref_line = x[0:]/2
        ref_line = x[0:]/1
        #plt.plot(x[0:],-ref_line-3.0,'-k')
        plt.plot(x[0:],-ref_line-0.0,'-k')
        plt.plot()
        #plt.legend(('Sw - FI','2$nd$ order'))
        plt.legend(('Sw FI','Primeira ordem'))
        plt.ylabel('$log_2({E}_{L_2})$')
        plt.xlabel('$log_2(N)$')
        plt.grid()
        plt.savefig('results/6k_L2_convergence_FI_Sw.png')


        plt.figure(5)
        x = np.log2(np.array([25, 50, 100, 200]))
        So_FI_L2 = np.log2(np.array([e25_So_L2_FI, e50_So_L2_FI, e100_So_L2_FI, e200_So_L2_FI]))

        plt.plot(x, So_FI_L2,'-ro')

        #ref_line = x[0:]/2
        ref_line = x[0:]/1
        #plt.plot(x[0:],-ref_line-3.0,'-k')
        plt.plot(x[0:],-ref_line-0.0,'-k')
        plt.plot()
        plt.legend(('So FI','Primeira ordem'))
        plt.ylabel('$log_2({E}_{L_2})$')
        plt.xlabel('$log_2(N)$')
        plt.grid()
        plt.savefig('results/6k_L2_convergence_FI_So.png')


        plt.figure(6)
        x = np.log2(np.array([25, 50, 100, 200]))
        Sg_FI_L2 = np.log2(np.array([e25_Sg_L2_FI, e50_Sg_L2_FI, e100_Sg_L2_FI, e200_Sg_L2_FI]))

        plt.plot(x, Sg_FI_L2,'-go')

        #ref_line = x[0:]/2
        ref_line = x[0:]/1
        #plt.plot(x[0:],-ref_line-6.0,'-k')
        plt.plot(x[0:],-ref_line-3.0,'-k')
        plt.plot()
        plt.legend(('Sg FI','Primeira ordem'))
        plt.ylabel('$log_2({E}_{L_2})$')
        plt.xlabel('$log_2(N)$')
        plt.grid()
        plt.savefig('results/6k_L2_convergence_FI_Sg.png')


        print('Done')
        import pdb; pdb.set_trace()
