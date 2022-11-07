import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy import interpolate

flying = 'flying'
name = 'results_'
arquivos = os.listdir(flying)

""" Analytical Solution """
a = 2000
b = 2000
l = 1000
q = 1000
h = 1

Pi = 2000
Bo = 1.1178
miu = 0.249
por = 0.2
k = 1.5
cf = 1.04e-5
cr = 5e-4
c = cf+cr
N = 100

Qj = 8.3
Q = -(1/Bo)*Qj/5.614/h

alpha = 157.952*por*c*miu/k
beta = 886.905*Bo*miu/k

t = 365
x = np.linspace(0,2000,N)
y = 840
P = np.zeros(N)

for i in range(0,N):
    d = 0
    for m in range(1,100+1):
        s = 1/(math.pi**2 * (m**2/a**2)) * (1 - np.exp(-math.pi**2/alpha*(m**2/a**2)*t))*np.cos(m*math.pi*l/a)*np.cos(m*math.pi*x[i]/a)
        d = d+s

    f=0
    for n in range(1,100+1):
        j = 1/(math.pi**2 * (n**2/b**2)) * (1 - np.exp(-math.pi**2/alpha*(n**2/b**2)*t))*np.cos(n*math.pi*q/b)*np.cos(n*math.pi*y/b)
        f = f+j

    g=0
    for n in range(1,100+1):
        for m in range(1,100+1):
            z = 1/(math.pi**2*(m**2/a**2+n**2/b**2)) * (1-np.exp(-math.pi**2/alpha*(m**2/a**2+n**2/b**2)*t)) * np.cos(m*math.pi*l/a) * \
            np.cos(n*math.pi*q/b) * np.cos(m*math.pi*x[i]/a) * np.cos(n*math.pi*y/b)
            g = g+z
    P[i] = Pi - beta*Q/(a*b)*(t/alpha+2*d+2*f+4*g)

P *= 6894.757/1e6
#    p_resp = np.linspace(0.623843,0,100)

for  arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/Li_Case4_2D/results_2d_Li_Case4_5x5x1_FI_CFL_0_5_456.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure_5x5x1 = data[4] / 1e6#6894.757
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[11]
            p0 = [0,195.072,-0.3048]
            p1 = [609.6,316.992,0.]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)

            pressure_5x5x1_1 = pressure_5x5x1[ind_ans]
            pressure_5x5x1_1 = pressure_5x5x1_1[ind_ans_sort]

            p0 = [0,121.92,-0.3048]
            p1 = [609.6,243.84,0.]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure_5x5x1_2 = pressure_5x5x1[ind_ans]
            pressure_5x5x1_2 = pressure_5x5x1_2[ind_ans_sort]

            pressure_5x5x1 = 0.6*(pressure_5x5x1_2 - pressure_5x5x1_1) + pressure_5x5x1_1
            x1 = np.linspace(0,2000,5)
            tck = interpolate.splrep(x1,pressure_5x5x1,s=0)
            p1 = interpolate.splev(x,tck,der=0)
            e1 = (sum((P-p1)**2)/(5*5))**(1/2)

        datas = np.load('flying/Li_Case4_2D/results_2d_Li_Case4_15x15x1_FI_CFL_0_5_456.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure_15x15x1 = data[4] / 1e6#6894.757
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[11]
            p0 = [0,235.712,-0.3048]
            p1 = [609.6,276.352,0.]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure_15x15x1_1 = pressure_15x15x1[ind_ans]
            pressure_15x15x1_1 = pressure_15x15x1_1[ind_ans_sort]

            p0 = [0,203.2,-0.3048]
            p1 = [609.6,243.84,0.]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure_15x15x1_2 = pressure_15x15x1[ind_ans]
            pressure_15x15x1_2 = pressure_15x15x1_2[ind_ans_sort]

            pressure_15x15x1 = 0.8*(pressure_15x15x1_2 - pressure_15x15x1_1) + pressure_15x15x1_1
            x2 = np.linspace(0,2000,15)

            tck = interpolate.splrep(x2,pressure_15x15x1,s=0)
            p2 = interpolate.splev(x,tck,der=0)

            e2 = (sum((P-p2)**2)/(15*15))**(1/2)
            R2 = np.log(e2/e1)/np.log(5/15)

        datas = np.load('flying/Li_Case4_2D/results_2d_Li_Case4_45x45x1_FI_CFL_0_5_456.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure_45x45x1 = data[4]  / 1e6#6894.757
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[11]
            p0 = [0,249.2586667,-0.3048]
            p1 = [609.6,262.8053333,0.0]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure_45x45x1_1 = pressure_45x45x1[ind_ans]
            pressure_45x45x1_1 = pressure_45x45x1_1[ind_ans_sort]

            p0 = [0,257.3866667,-0.3048]
            p1 = [609.6,270.9333333,0.0]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure_45x45x1_2 = pressure_45x45x1[ind_ans]
            pressure_45x45x1_2 = pressure_45x45x1_2[ind_ans_sort]

            pressure_45x45x1 = 0.40059*(pressure_45x45x1_2 - pressure_45x45x1_1) + pressure_45x45x1_1
            x3 = np.linspace(0,2000,45)
            tck = interpolate.splrep(x3,pressure_45x45x1,s=0)
            p3 = interpolate.splev(x,tck,der=0)

            e3 = (sum((P-p3)**2)/(45*45))**(1/2)
            R3 = np.log(e3/e2)/np.log(15/45)


        #    p_resp = np.linspace(0.623843,0,100)
        sizeletter = 12
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300

        plt.rcParams.update({'font.size': sizeletter})
        plt.figure(1)
        #plt.title('t = 365 dias')
        plt.plot(x1*0.3048, pressure_5x5x1, '-ro', x2*0.3048, pressure_15x15x1, '-bs', \
            x3*0.3048, pressure_45x45x1, '-gP', x*0.3048, P, 'k', mfc='none', markersize=5)
        plt.legend(('FI-5x5x1', 'FI-15x15x1', 'FI-45x45x1',
            'Solução Analítica'), prop={'size': sizeletter-1})
        #plt.figure(2)
        #plt.plot( x4, pressure4, 'g', x, P, 'y')
        plt.grid()
        plt.ylabel('Pressão [MPa]')
        plt.xlabel('Distância na direção x [m]')
        plt.savefig('results/pressure_2d_FI.png')

        import pdb; pdb.set_trace()
