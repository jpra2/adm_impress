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
n_400 = 400
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

Sw_CMG_400 = np.loadtxt('Sw_points_CMG_SchmallNEW_400.txt')
So_CMG_400 = np.loadtxt('So_points_CMG_SchmallNEW_400.txt')
Sg_CMG_400 = np.loadtxt('Sg_points_CMG_SchmallNEW_400.txt')

Sw_CMG_200 = np.loadtxt('Sw_points_CMG_SchmallNEW_200.txt')
So_CMG_200 = np.loadtxt('So_points_CMG_SchmallNEW_200.txt')
Sg_CMG_200 = np.loadtxt('Sg_points_CMG_SchmallNEW_200.txt')

Sw_CMG_100 = np.loadtxt('Sw_points_CMG_SchmallNEW_100.txt')
So_CMG_100 = np.loadtxt('So_points_CMG_SchmallNEW_100.txt')
Sg_CMG_100 = np.loadtxt('Sg_points_CMG_SchmallNEW_100.txt')

Sw_CMG_50 = np.loadtxt('Sw_points_CMG_SchmallNEW_50.txt')
So_CMG_50 = np.loadtxt('So_points_CMG_SchmallNEW_50.txt')
Sg_CMG_50 = np.loadtxt('Sg_points_CMG_SchmallNEW_50.txt')

Sw_CMG_25 = np.loadtxt('Sw_points_CMG_SchmallNEW_25.txt')
So_CMG_25 = np.loadtxt('So_points_CMG_SchmallNEW_25.txt')
Sg_CMG_25 = np.loadtxt('Sg_points_CMG_SchmallNEW_25.txt')

fzC1 = open('ZC1_6k.txt','r')
zC1_CMG = [float(line.rstrip('\n\r')) for line in fzC1]

fzC20 = open('ZC20_6k.txt','r')
zC20_CMG = [float(line.rstrip('\n\r')) for line in fzC20]

for arq in arquivos:
    if  arq.startswith(name):


        datas3 = np.load('flying/6k/25_0_0/results_6k_SchmallNEW_25_FI_20days_POCONEW_40.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI_25 = data[5]
            So_FI_25 = data[6]
            Sg_FI_25 = data[7]
            zC1_FI_25 = data[10][0]
            zC20_FI_25 = data[10][5]
            x1_25 = np.linspace(0, 50, n_25)

            e25_Sw_L1_FI = (sum(abs(f_Sw(x1_25)-Sw_FI_25))*(1/25))
            e25_Sw_L2_FI = np.sqrt(np.sum((f_Sw(x1_25)-Sw_FI_25)**2) * 1 / 25)

            e25_So_L1_FI = (sum(abs(f_So(x1_25)-So_FI_25))*(1/25))
            e25_So_L2_FI = np.sqrt(np.sum((f_So(x1_25)-So_FI_25)**2) * 1 / 25)

            e25_Sg_L1_FI = (sum(abs(f_Sg(x1_25)-Sg_FI_25))*(1/25))
            e25_Sg_L2_FI = np.sqrt(np.sum((f_Sg(x1_25)-Sg_FI_25)**2) * 1 / 25)


            # New
            e25_Sw_L1_FI_new = (sum(abs(Sw_CMG_25-Sw_FI_25))*(1/25))
            e25_Sw_L2_FI_new = np.sqrt(np.sum((Sw_CMG_25-Sw_FI_25)**2) * 1 / 25)

            e25_So_L1_FI_new = (sum(abs(So_CMG_25-So_FI_25))*(1/25))
            e25_So_L2_FI_new = np.sqrt(np.sum((So_CMG_25-So_FI_25)**2) * 1 / 25)

            e25_Sg_L1_FI_new = (sum(abs(Sg_CMG_25-Sg_FI_25))*(1/25))
            e25_Sg_L2_FI_new = np.sqrt(np.sum((Sg_CMG_25-Sg_FI_25)**2) * 1 / 25)




        datas3 = np.load('flying/6k/results_6k_SchmallNEW_50_FI_20days_2DTeste_203.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI_50 = data[5]
            So_FI_50 = data[6]
            Sg_FI_50 = data[7]
            zC1_FI_50 = data[10][0]
            zC20_FI_50 = data[10][5]
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



            # New
            e50_Sw_L1_FI_new = (sum(abs(Sw_CMG_50-Sw_FI_50))*(1/50))
            R50_Sw_L1_FI_new = math.log(e25_Sw_L1_FI_new/e50_Sw_L1_FI_new,2)
            e50_Sw_L2_FI_new = np.sqrt(np.sum((Sw_CMG_50-Sw_FI_50)**2) * 1 / 50)
            R50_Sw_L2_FI_new = math.log(e25_Sw_L2_FI_new/e50_Sw_L2_FI_new,2)

            e50_So_L1_FI_new = (sum(abs(So_CMG_50-So_FI_50))*(1/50))
            R50_So_L1_FI_new = math.log(e25_So_L1_FI_new/e50_So_L1_FI_new,2)
            e50_So_L2_FI_new = np.sqrt(np.sum((So_CMG_50-So_FI_50)**2) * 1 / 50)
            R50_So_L2_FI_new = math.log(e25_So_L2_FI_new/e50_So_L2_FI_new,2)

            e50_Sg_L1_FI_new = (sum(abs(Sg_CMG_50-Sg_FI_50))*(1/50))
            R50_Sg_L1_FI_new = math.log(e25_Sg_L1_FI_new/e50_Sg_L1_FI_new,2)
            e50_Sg_L2_FI_new = np.sqrt(np.sum((Sg_CMG_50-Sg_FI_50)**2) * 1 / 50)
            R50_Sg_L2_FI_new = math.log(e25_Sg_L2_FI_new/e50_Sg_L2_FI_new,2)

        datas2 = np.load('flying/6k/100_0_0/results_6k_SchmallNEW_100_FI_20days_POCONEW_dtmenor_2003.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas2[1:]:
            Sw_FI_100 = data[5]
            So_FI_100 = data[6]
            Sg_FI_100 = data[7]
            zC1_FI_100 = data[10][0]
            zC20_FI_100 = data[10][5]
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



            # New
            e100_Sw_L1_FI_new = (sum(abs(Sw_CMG_100-Sw_FI_100))*(1/100))
            R100_Sw_L1_FI_new = math.log(e50_Sw_L1_FI_new/e100_Sw_L1_FI_new,2)
            e100_Sw_L2_FI_new = np.sqrt(np.sum((Sw_CMG_100-Sw_FI_100)**2) * 1 / 100)
            R100_Sw_L2_FI_new = math.log(e50_Sw_L2_FI_new/e100_Sw_L2_FI_new,2)

            e100_So_L1_FI_new = (sum(abs(So_CMG_100-So_FI_100))*(1/100))
            R100_So_L1_FI_new = math.log(e50_So_L1_FI_new/e100_So_L1_FI_new,2)
            e100_So_L2_FI_new = np.sqrt(np.sum((So_CMG_100-So_FI_100)**2) * 1 / 100)
            R100_So_L2_FI_new = math.log(e50_So_L2_FI_new/e100_So_L2_FI_new,2)

            e100_Sg_L1_FI_new = (sum(abs(Sg_CMG_100-Sg_FI_100))*(1/100))
            R100_Sg_L1_FI_new = math.log(e50_Sg_L1_FI_new/e100_Sg_L1_FI_new,2)
            e100_Sg_L2_FI_new = np.sqrt(np.sum((Sg_CMG_100-Sg_FI_100)**2) * 1 / 100)
            R100_Sg_L2_FI_new = math.log(e50_Sg_L2_FI_new/e100_Sg_L2_FI_new,2)

        datas3 = np.load('flying/6k/200_0_0/results_6k_SchmallNEW_200_FI_20days_POCONEW_2003.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas3[1:]:
            Sw_FI_200 = data[5]
            So_FI_200 = data[6]
            Sg_FI_200 = data[7]
            zC1_FI_200 = data[10][0]
            zC20_FI_200 = data[10][5]
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


            # New
            e200_Sw_L1_FI_new = (sum(abs(Sw_CMG_200-Sw_FI_200))*(1/200))
            R200_Sw_L1_FI_new = math.log(e100_Sw_L1_FI_new/e200_Sw_L1_FI_new,2)
            e200_Sw_L2_FI_new = np.sqrt(np.sum((Sw_CMG_200-Sw_FI_200)**2) * 1 / 200)
            R200_Sw_L2_FI_new = math.log(e100_Sw_L2_FI_new/e200_Sw_L2_FI_new,2)

            e200_So_L1_FI_new = (sum(abs(So_CMG_200-So_FI_200))*(1/200))
            R200_So_L1_FI_new = math.log(e100_So_L1_FI_new/e200_So_L1_FI_new,2)
            e200_So_L2_FI_new = np.sqrt(np.sum((So_CMG_200-So_FI_200)**2) * 1 / 200)
            R200_So_L2_FI_new = math.log(e100_So_L2_FI_new/e200_So_L2_FI_new,2)

            e200_Sg_L1_FI_new = (sum(abs(Sg_CMG_200-Sg_FI_200))*(1/200))
            R200_Sg_L1_FI_new = math.log(e100_Sg_L1_FI_new/e200_Sg_L1_FI_new,2)
            e200_Sg_L2_FI_new = np.sqrt(np.sum((Sg_CMG_200-Sg_FI_200)**2) * 1 / 200)
            R200_Sg_L2_FI_new = math.log(e100_Sg_L2_FI_new/e200_Sg_L2_FI_new,2)



        datas4 = np.load('flying/results_6k_SchmallNEW_400_FI_20days_CFL3_0_558.npy', allow_pickle=True)
        #import pdb; pdb.set_trace()
        for data in datas4[1:]:
            Sw_FI_400 = data[5]
            So_FI_400 = data[6]
            Sg_FI_400 = data[7]
            zC1_FI_400 = data[10][0]
            zC20_FI_400 = data[10][5]
            x1_400 = np.linspace(0, 50, n_400)

            e400_Sw_L1_FI = (sum(abs(f_Sw(x1_400)-Sw_FI_400))*(1/400))
            R400_Sw_L1_FI = math.log(e200_Sw_L1_FI/e400_Sw_L1_FI,2)
            e400_Sw_L2_FI = np.sqrt(np.sum((f_Sw(x1_400)-Sw_FI_400)**2) * 1 / 400)
            R400_Sw_L2_FI = math.log(e200_Sw_L2_FI/e400_Sw_L2_FI,2)

            e400_So_L1_FI = (sum(abs(f_So(x1_400)-So_FI_400))*(1/400))
            R400_So_L1_FI = math.log(e200_So_L1_FI/e400_So_L1_FI,2)
            e400_So_L2_FI = np.sqrt(np.sum((f_So(x1_400)-So_FI_400)**2) * 1 / 400)
            R400_So_L2_FI = math.log(e200_So_L2_FI/e400_So_L2_FI,2)

            e400_Sg_L1_FI = (sum(abs(f_Sg(x1_400)-Sg_FI_400))*(1/400))
            R400_Sg_L1_FI = math.log(e200_Sg_L1_FI/e400_Sg_L1_FI,2)
            e400_Sg_L2_FI = np.sqrt(np.sum((f_Sg(x1_400)-Sg_FI_400)**2) * 1 / 400)
            R400_Sg_L2_FI = math.log(e200_Sg_L2_FI/e400_Sg_L2_FI,2)

            # New
            e400_Sw_L1_FI_new = (sum(abs(Sw_CMG_400-Sw_FI_400))*(1/400))
            R400_Sw_L1_FI_new = math.log(e200_Sw_L1_FI_new/e400_Sw_L1_FI_new,2)
            e400_Sw_L2_FI_new = np.sqrt(np.sum((Sw_CMG_400-Sw_FI_400)**2) * 1 / 400)
            R400_Sw_L2_FI_new = math.log(e200_Sw_L2_FI_new/e400_Sw_L2_FI_new,2)

            e400_So_L1_FI_new = (sum(abs(So_CMG_400-So_FI_400))*(1/400))
            R400_So_L1_FI_new = math.log(e200_So_L1_FI_new/e400_So_L1_FI_new,2)
            e400_So_L2_FI_new = np.sqrt(np.sum((So_CMG_400-So_FI_400)**2) * 1 / 400)
            R400_So_L2_FI_new = math.log(e200_So_L2_FI_new/e400_So_L2_FI_new,2)

            e400_Sg_L1_FI_new = (sum(abs(Sg_CMG_400-Sg_FI_400))*(1/400))
            R400_Sg_L1_FI_new = math.log(e200_Sg_L1_FI_new/e400_Sg_L1_FI_new,2)
            e400_Sg_L2_FI_new = np.sqrt(np.sum((Sg_CMG_400-Sg_FI_400)**2) * 1 / 400)
            R400_Sg_L2_FI_new = math.log(e200_Sg_L2_FI_new/e400_Sg_L2_FI_new,2)



        plt.figure(1)
        x = np.log2(np.array([25, 50, 100, 200]))
        Sw_FI = np.log2(np.array([e25_Sw_L1_FI, e50_Sw_L1_FI, e100_Sw_L1_FI, e200_Sw_L1_FI]))
        #Sw_FI = np.log2(np.array([e25_Sw_L1_FI, e50_Sw_L1_FI, e100_Sw_L1_FI, e200_Sw_L1_FI, e400_Sw_L1_FI]))
        #Sw_FI = np.log2(np.array([e25_Sw_L1_FI_new, e50_Sw_L1_FI_new, e100_Sw_L1_FI_new, e200_Sw_L1_FI_new, e400_Sw_L1_FI_new]))

        y_ref = -x-2.0

        plt.plot(x, Sw_FI, '-bo')
        plt.plot(x, y_ref, 'k')
        #plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{2}({E}_{L_1})$')
        plt.xlabel('$log_{2}(N)$')
        plt.legend(('Sw FI', 'Convergência linear'))
        plt.grid()
        plt.savefig('results/6k/6k_400/6k_L1_convergence_FI_Sw_new_.png')



        plt.figure(2)
        x = np.log2(np.array([25, 50, 100, 200]))
        So_FI = np.log2(np.array([e25_So_L1_FI, e50_So_L1_FI, e100_So_L1_FI, e200_So_L1_FI]))
        #So_FI = np.log2(np.array([e25_So_L1_FI, e50_So_L1_FI, e100_So_L1_FI, e200_So_L1_FI, e400_So_L1_FI]))
        #So_FI = np.log2(np.array([e25_So_L1_FI_new, e50_So_L1_FI_new, e100_So_L1_FI_new, e200_So_L1_FI_new, e400_So_L1_FI_new]))

        y_ref = -x-2.0

        plt.plot(x, So_FI, '-ro')
        plt.plot(x, y_ref, 'k')
        #plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{2}({E}_{L_1})$')
        plt.xlabel('$log_{2}(N)$')
        plt.legend(('So FI', 'Convergência linear'))
        plt.grid()
        plt.savefig('results/6k/6k_400/6k_L1_convergence_FI_So_new.png')


        plt.figure(3)
        x = np.log2(np.array([25, 50, 100, 200]))
        Sg_FI = np.log2(np.array([e25_Sg_L1_FI, e50_Sg_L1_FI, e100_Sg_L1_FI, e200_Sg_L1_FI]))
        #Sg_FI = np.log2(np.array([e25_Sg_L1_FI, e50_Sg_L1_FI, e100_Sg_L1_FI, e200_Sg_L1_FI, e400_Sg_L1_FI]))
        #Sg_FI = np.log2(np.array([e25_Sg_L1_FI_new, e50_Sg_L1_FI_new, e100_Sg_L1_FI_new, e200_Sg_L1_FI_new, e400_Sg_L1_FI_new]))

        y_ref = -x-4.0

        plt.plot(x, Sg_FI, '-go')
        plt.plot(x, y_ref, 'k')
        #plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{2}({E}_{L_1})$')
        plt.xlabel('$log_{2}(N)$')
        plt.legend(('Sg FI', 'Convergência linear'))
        plt.grid()
        plt.savefig('results/6k/6k_400/6k_L1_convergence_FI_Sg_new.png')


        """
        plt.figure(4)
        x = np.log2(np.array([25, 50, 100, 200, 400]))
        #Sw_FI_L2 = np.log2(np.array([e25_Sw_L2_FI, e50_Sw_L2_FI, e100_Sw_L2_FI, e200_Sw_L2_FI, e400_Sw_L2_FI]))
        Sw_FI_L2 = np.log2(np.array([e25_Sw_L2_FI_new, e50_Sw_L2_FI_new, e100_Sw_L2_FI_new, e200_Sw_L2_FI_new, e400_Sw_L2_FI_new]))

        plt.plot(x, Sw_FI_L2,'-bo')

        #ref_line = x[0:]/2
        ref_line = x[0:]/1
        #plt.plot(x[0:],-ref_line-3.0,'-k')
        plt.plot(x[0:],-ref_line-0.5,'-k')
        plt.plot()
        #plt.legend(('Sw - FI','2$nd$ order'))
        plt.legend(('Sw FI','Primeira ordem'))
        plt.ylabel('$log_2({E}_{L_2})$')
        plt.xlabel('$log_2(N)$')
        plt.grid()
        plt.savefig('results/6k/6k_400/6k_L2_convergence_FI_Sw_new.png')


        plt.figure(5)
        x = np.log2(np.array([25, 50, 100, 200, 400]))
        #So_FI_L2 = np.log2(np.array([e25_So_L2_FI, e50_So_L2_FI, e100_So_L2_FI, e200_So_L2_FI, e400_So_L2_FI]))
        So_FI_L2 = np.log2(np.array([e25_So_L2_FI_new, e50_So_L2_FI_new, e100_So_L2_FI_new, e200_So_L2_FI_new, e400_So_L2_FI_new]))

        plt.plot(x, So_FI_L2,'-ro')

        #ref_line = x[0:]/2
        ref_line = x[0:]/1
        #plt.plot(x[0:],-ref_line-3.0,'-k')
        plt.plot(x[0:],-ref_line-0.5,'-k')
        plt.plot()
        plt.legend(('So FI','Primeira ordem'))
        plt.ylabel('$log_2({E}_{L_2})$')
        plt.xlabel('$log_2(N)$')
        plt.grid()
        plt.savefig('results/6k/6k_400/6k_L2_convergence_FI_So_new.png')


        plt.figure(6)
        x = np.log2(np.array([25, 50, 100, 200, 400]))
        #Sg_FI_L2 = np.log2(np.array([e25_Sg_L2_FI, e50_Sg_L2_FI, e100_Sg_L2_FI, e200_Sg_L2_FI, e400_Sg_L2_FI]))
        Sg_FI_L2 = np.log2(np.array([e25_Sg_L2_FI_new, e50_Sg_L2_FI_new, e100_Sg_L2_FI_new, e200_Sg_L2_FI_new, e400_Sg_L2_FI_new]))

        plt.plot(x, Sg_FI_L2,'-go')

        #ref_line = x[0:]/2
        ref_line = x[0:]/1
        #plt.plot(x[0:],-ref_line-6.0,'-k')
        plt.plot(x[0:],-ref_line-3.5,'-k')
        plt.plot()
        plt.legend(('Sg FI','Primeira ordem'))
        plt.ylabel('$log_2({E}_{L_2})$')
        plt.xlabel('$log_2(N)$')
        plt.grid()
        plt.savefig('results/6k/6k_400/6k_L2_convergence_FI_Sg_new.png')
        """

        plt.figure(7)
        #plt.title('t = 10 days - 50x1x1 mesh')
        #plt.plot(x1, zC1_FOU, 'k')
        plt.plot(x_CMG, zC1_CMG, 'k')
        #plt.plot(x1_50, zC1_FI_50, '--y')
        plt.plot(x1_100, zC1_FI_100, ':b')
        plt.plot(x1_200, zC1_FI_200, '-.g')
        plt.plot(x1_400, zC1_FI_400, '--r')
        #plt.legend(('CMG - 5000 CVs', 'FI 50 CVs', 'FI 100 CVs', 'FI 200 CVs', 'FI 400 CVs'), loc=1)
        plt.legend(('CMG - 5000 CVs', 'FI 50 CVs', 'FI 100 CVs', 'FI 200 CVs'), loc=1)
        plt.ylabel('Fração molar global do C1')
        plt.xlabel('Distância (m)')
        plt.grid()
        plt.ylim((0.38, 0.45))
        plt.savefig('results/6k/6k_400/zC1_6k_FI_20days_lim' + '.png')


        plt.figure(8)
        #plt.title('t = 10 days - 50x1x1 mesh')
        #plt.plot(x1, zC1_FOU, 'k')
        plt.plot(x_CMG, zC20_CMG, 'k')
        #plt.plot(x1_50, zC20_FI_50, '--y')
        plt.plot(x1_100, zC20_FI_100, ':b')
        plt.plot(x1_200, zC20_FI_200, '-.g')
        plt.plot(x1_400, zC20_FI_400, '--r')
        #plt.legend(('CMG - 5000 CVs', 'FI 50 CVs', 'FI 100 CVs', 'FI 200 CVs', 'FI 400 CVs'), loc=1)
        plt.legend(('CMG - 5000 CVs', 'FI 50 CVs', 'FI 100 CVs', 'FI 200 CVs'), loc=4)
        plt.ylabel('Fração molar global do C20')
        plt.xlabel('Distância (m)')
        plt.grid()
        plt.ylim((0.055, 0.0625))
        plt.savefig('results/6k/6k_400/zC20_6k_FI_20days_lim' + '.png')


        print('Done')
        import pdb; pdb.set_trace()
