import numpy as np
import matplotlib.pyplot as plt
import os
import math

flying = 'flying'
name = 'results_'
arquivos = os.listdir(flying)
i=1

Q = 0.04 * 0.3048**3/86400#0.178107607
g = 9.80665
L = 2*0.3048
rho = 44.7*0.45359237/0.3048**3
gamma = 44.7/62.42
A = 0.01*0.3048**2
k = 500/1000*9.869233*10**(-13)
por = 0.2
vis = 0.249*10**(-3)
n = 100
l = L/100
P = np.zeros(n)

P[99] = 2000.*6894.757

x_ans = np.linspace(0,1,100)

for i in range(1,n):
    #P[n-1-i] = P[n-i] + (Q*vis/(1.1271*k*A)-0.43333333333*gamma)*l
    P[n-1-i] = P[n-i] + (Q*vis/(k*A)- g*rho)*l

for arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_incompressive_vertical_case_22.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure = (data[4] - data[4][99] * np.ones(100))/6894.757

            P = (P - P[n-1])/6894.75729
            time = data[3]
            loop = data[0]
            #flux = data[5]
            #flux_vector = data[6]
            i=loop
            #flux_vols = data[5]
            x = np.linspace(0,1,100)
        #    p_resp = np.linspace(0.623843,0,100)
            plt.figure(1)
            #plt.title('t = 0.01 days')
            plt.plot(x, pressure, 'r', x_ans, P, 'g')
            plt.grid()
            plt.legend(('PADMEC', 'Analytical Solution'))
            plt.ylabel('Pressure Drop (psi)')
            plt.xlabel('Dimensionless distance')
            plt.title('Water Vertical Displacement')
            plt.savefig('results/compositional/pressure_vert_analytical_w' + str(loop) + '.png')
            import pdb; pdb.set_trace()
