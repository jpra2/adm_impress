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



        datas = np.load('flying/results_2d_injection_case_4_Li_IMPEC_FOU_452.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure_FOU = data[4]/1e6#6894.757
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[13]

            p0 = [0,243.84,-0.3048] #[0,243.84,-0.3048]
            p1 = [609.6,268.224,0.0] #[609.6,268.224,0.0]
            import pdb; pdb.set_trace()
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure_FOU_1 = pressure_FOU[ind_ans]
            pressure_FOU = pressure_FOU_1[ind_ans_sort]

            #pressure4 = 0.2*(pressure4_2-pressure4_1) + pressure4_1
            x4 = np.linspace(0,2000,25)
            #tck = interpolate.splrep(x4,pressure_FOU,s=0)
            #p4 = interpolate.splev(x,tck,der=0)
            #e4 = (sum((P-p4)**2)/(25*25))**(1/2)

        datas = np.load('flying/results_2d_injection_case_4_Li_IMPEC_MUSCL_452.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure_MUSCL = data[4]/1e6#6894.757
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[13]

            p0 = [0,292.608,-0.3048] #[0,243.84,-0.3048]
            p1 = [609.6,316.992,0.0] #[609.6,268.224,0.0]
            import pdb; pdb.set_trace()
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure_MUSCL_1 = pressure_MUSCL[ind_ans]
            pressure_MUSCL = pressure_MUSCL_1[ind_ans_sort]

            #pressure4 = 0.2*(pressure4_2-pressure4_1) + pressure4_1
            x4 = np.linspace(0,2000,25)
            tck = interpolate.splrep(x4,pressure_MUSCL,s=0)
            p4 = interpolate.splev(x,tck,der=0)
            e4 = (sum((P-p4)**2)/(25*25))**(1/2)

        datas = np.load('flying/results_2d_injection_case_4_Li_IMPEC_MUSCLu_452.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure_MUSCLu = data[4]/1e6#6894.757
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[13]

            p0 = [0,292.608,-0.3048] #[0,243.84,-0.3048]
            p1 = [609.6,316.992,0.0] #[609.6,268.224,0.0]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure_MUSCLu_1 = pressure_MUSCLu[ind_ans]
            pressure_MUSCLu = pressure_MUSCLu_1[ind_ans_sort]

            #pressure4 = 0.2*(pressure4_2-pressure4_1) + pressure4_1
            x42u = np.linspace(0,2000,len(pressure_MUSCLu))
            #x4 = cent_ind
            #tck = interpolate.splrep(x4,pressure_MUSCL,s=0)
            #p4 = interpolate.splev(x,tck,der=0)
            #e4 = (sum((P-p4)**2)/(25*25))**(1/2)

        #    p_resp = np.linspace(0.623843,0,100)
        sizeletter = 12
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300
        import pdb; pdb.set_trace()

        plt.rcParams.update({'font.size': sizeletter})
        plt.figure(1)
        #plt.title('t = 365 dias')
        #plt.plot(x4*0.3048, pressure_FOU, 'yo', x4*0.3048, pressure_MUSCL, 'g',x4*0.3048, pressure_MUSCLu, 'r', x*0.3048, P, 'k', mfc='none', markersize=5)
        plt.plot(x4*0.3048, pressure_FOU, 'yo', x4*0.3048, pressure_MUSCL, 'gs', x*0.3048, P, 'k', mfc='none', markersize=5)
        #plt.plot(x4*0.3048, pressure_FOU, 'yo', x*0.3048, P, 'k', mfc='none', markersize=5)
        plt.legend(('IMPEC FOU-2025', 'IMPEC MUSCL-2025', 'IMPEC MUSCLu-2025', 'Solução Analítica'), prop={'size': sizeletter-1})
        #plt.figure(2)
        #plt.plot( x4, pressure4, 'g', x, P, 'y')
        plt.grid()
        plt.ylabel('Pressão [MPa]')
        plt.xlabel('Distância na direção x [m]')
        plt.savefig('results/compositional/pressure_2d_case4_Li_FOU.png')

        import pdb; pdb.set_trace()
