import numpy as np
import matplotlib.pyplot as plt
import os
import math

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
i=1

Q = 0.04 * 0.3048**3/86400
L = 2*0.3048
rho = 44.7*0.45359237/0.3048**3
gamma = 44.7/62.42
A = 0.01*0.3048**2
k = 300/1000*9.869233*10**(-13)
por = 0.2
vis = 0.249*10**(-3)
n = 100
l = L/100
P = np.zeros(n)
H = np.zeros(n)
P[99] = 2000.*6894.757
H = np.zeros(n)
#H[0] = 0.01

for i in range(1,n):
    P[n-1-i] = P[n-i] + (Q*vis*l/k/A)

x_ans = np.linspace(0,1,100)

Q1 = 0.04/86400#*0.178107607
vis = 0.00249*0.000014503774389728
l = 0.02
A = 0.01
k = 0.5*0.9869233*10**(-12)/0.3048**2
P1 = np.zeros(n)
P1[99] = 2000
for i in range(1,n):
    P1[n-1-i] = P1[n-i] + (Q1*vis*l/k/A)
#import pdb; pdb.set_trace()

for arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_monophasic_compressiblew_6.npy', allow_pickle=True)

        pressure = np.zeros([6,100])
        time = np.zeros(6)
        b=0
        for data in datas[1:]:
            pressure[b,:] = (data[4] - data[4][99] * np.ones(100))/6894.75729
            P = (P - P[n-1])/6894.75729
            time[b] = data[3]
            loop = data[0]
            #flux = data[5]
            #flux_vector = data[6]
            i=loop
            #flux_vols = data[5]
            x = np.linspace(0,1,100)
            b = b+1

        #    p_resp = np.linspace(0.623843,0,100)
        plt.figure(1)
        #plt.title('t = 0.01 days')

        for i in range(0,5):
            plt.plot(x, pressure[i,:], label = '{}s'.format(time[i]))#, x_ans, P, 'g')
        plt.plot(x, pressure[5,:], label = 'solution')
        plt.grid()
        plt.legend(bbox_to_anchor=(0.95, 0.8), loc=7, borderaxespad=0., ncol = 2, handletextpad = 0.2)
        #plt.legend(('PADMEC', 'Analytical Solution'))
        plt.ylabel('Pressure Drop (psi)')
        plt.xlabel('Dimensionless distance')
        plt.title('Water horizontal displacement')
        plt.savefig('results/compositional/pressure_hor_analytical_comp_w' + str(loop) + '.png')
        import pdb; pdb.set_trace()
