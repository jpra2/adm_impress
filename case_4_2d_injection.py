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


#    p_resp = np.linspace(0.623843,0,100)
for  arq in arquivos:
    if  arq.startswith(name):

        datas = np.load('flying/results_2d_injection_5_case_452.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure1 = data[4] / 6894.757
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[6]
            p0 = [0,195.072,-0.3048]
            p1 = [609.6,316.992,0.]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure1_1 = pressure1[ind_ans]
            pressure1_1 = pressure1_1[ind_ans_sort]

            p0 = [0,121.9,-0.3048]
            p1 = [609.6,243.8,0.]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure1_2 = pressure1[ind_ans]
            pressure1_2 = pressure1_2[ind_ans_sort]

            pressure1 = 0.6*(pressure1_2 - pressure1_1) + pressure1_1
            x1 = np.linspace(0,2000,5)
            tck = interpolate.splrep(x1,pressure1,s=0)
            p1 = interpolate.splev(x,tck,der=0)
            e1 = (sum((P-p1)**2)/(5*5))**(1/2)

        datas = np.load('flying/results_2d_injection_15_case_452.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure2 = data[4] / 6894.757
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[6]
            p0 = [0,235.712,-0.3048]
            p1 = [609.6,276.712,0.]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure2_1 = pressure2[ind_ans]
            pressure2_1 = pressure2_1[ind_ans_sort]

            p0 = [0,203.2,-0.3048]
            p1 = [609.6,243.84,0.]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure2_2 = pressure2[ind_ans]
            pressure2_2 = pressure2_2[ind_ans_sort]

            pressure2 = 0.8*(pressure2_2 - pressure2_1) + pressure2_1
            x2 = np.linspace(0,2000,15)

            tck = interpolate.splrep(x2,pressure2,s=0)
            p2 = interpolate.splev(x,tck,der=0)

            e2 = (sum((P-p2)**2)/(15*15))**(1/2)
            R2 = np.log(e2/e1)/np.log(5/15)

        datas = np.load('flying/results_2d_injection_45_case_452.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure3 = data[4] / 6894.757
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[6]
            p0 = [0,249.2586667,-0.3048]
            p1 = [609.6,262.8053333,0.0]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure3_1 = pressure3[ind_ans]
            pressure3_1 = pressure3_1[ind_ans_sort]

            p0 = [0,257.3866667,-0.3048]
            p1 = [609.6,270.9333333,0.0]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure3_2 = pressure3[ind_ans]
            pressure3_2 = pressure3_2[ind_ans_sort]

            pressure3 = 0.40059*(pressure3_2 - pressure3_1) + pressure3_1
            x3 = np.linspace(0,2000,45)
            tck = interpolate.splrep(x3,pressure3,s=0)
            p3 = interpolate.splev(x,tck,der=0)

            e3 = (sum((P-p3)**2)/(45*45))**(1/2)
            R3 = np.log(e3/e2)/np.log(15/45)

        datas = np.load('flying/results_2d_injection_25_case_452.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure4 = data[4] / 6894.75729
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[6]
            p0 = [0,243.84,-0.3048]
            p1 = [609.6,268.224,0.0]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure4_1 = pressure4[ind_ans]
            pressure4 = pressure4_1[ind_ans_sort]
            #pressure4 = 0.2*(pressure4_2-pressure4_1) + pressure4_1
            x4 = np.linspace(0,2000,25)
            tck = interpolate.splrep(x4,pressure4,s=0)
            p4 = interpolate.splev(x,tck,der=0)
            e4 = (sum((P-p4)**2)/(25*25))**(1/2)

        '''datas = np.load('flying/results_2d_injection_case_35__74.npy', allow_pickle=True)

        for data in datas[1:]:
            pressure5 = data[4] / 6894.757
            """ Just for the 2D case """
            from packs.utils.utils_old import get_box
            centroids = data[6]
            p0 = [0,247.3234286,0]
            p1 = [609.6,264.74057,0.3048]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure5_1 = pressure5[ind_ans]
            pressure5_1 = pressure5_1[ind_ans_sort]

            p0 = [0,261.2571429,0]
            p1 = [609.6,278.6742857,0.3048]
            ind_ans = get_box(centroids,np.array([p0,p1]))
            cent_ind = centroids[ind_ans]
            cent_mix = cent_ind[:,0]
            ind_ans_sort = np.argsort(cent_mix)
            pressure5_2 = pressure5[ind_ans]
            pressure5_2 = pressure5_2[ind_ans_sort]

            pressure5 = 0.2*(pressure5_2-pressure5_1) + pressure5_1
            x5 = np.linspace(0,2000,35)
            tck = interpolate.splrep(x5,pressure5,s=0)
            p5 = interpolate.splev(x,tck,der=0)
            e5 = (sum((P-p5)**2)/(25*25))**(1/2)'''


        #    p_resp = np.linspace(0.623843,0,100)
        plt.figure(1)
        plt.title('t = 365 days')
        plt.plot(x1, pressure1, 'r', x2, pressure2, 'b', x3, pressure3, 'g', x, P, 'y')

        plt.legend(('25 blocks', '225 blocks', '2025 blocks', 'Analytical Solution'))
        #plt.figure(2)
        #plt.plot( x4, pressure4, 'g', x, P, 'y')
        plt.grid()
        plt.ylabel('Pressure (psi)')
        plt.xlabel('Distance in X - Direction (ft)')
        plt.savefig('results/compositional/pressure_2d1_'  + '.png')
        import pdb; pdb.set_trace()
