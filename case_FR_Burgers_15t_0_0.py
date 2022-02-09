import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
import quadpy
from scipy.optimize import root_scalar


flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
t = 1.5

def f(x, pos, t):
    return x - (1/(2*math.pi))*np.sin(2*math.pi*(pos-t*x))

'''datas = np.load('flying/results_Burger_10000_15t_upw_19894.npy', allow_pickle=True)

for data in datas[1:]:
    Nk_FOU = data[12][0]'''

x = np.linspace(0+1/4000,1-1/4000, 2000)
Nk_ans = np.empty_like(x)
for i in range(2000):
    Nk_ans[i] = root_scalar(f, args=(x[i], t), method='toms748', bracket=[-1, 1]).root


plt.figure(1)
plt.plot(x,Nk_ans)
plt.grid()
plt.savefig('results/compositional/FR_paper/Burger/Burguer_ref_solution.png')
with open('Sw_Burguer_ref_sol.txt', 'w') as f2:
    for item in Nk_ans:
        f2.write("%s\n" % item)

with open('x_Burguer_re_sol.txt', 'w') as f2:
    for item in x:
        f2.write("%s\n" % item)


for arq in arquivos:
    if  arq.startswith(name):

        '---------------------------Flux Reconstruction------------------------'
        datas = np.load('flying/results_Burger_8_15t_FR2_1501.npy', allow_pickle=True)


        for data in datas[1:]:
            Nk8_FR = data[12][0].flatten()
            n = 8
            x8_2 = np.empty((n,2))
            GL = quadpy.c1.gauss_lobatto(2)
            points = GL.points
            for i in range(len(points)):
                x8_2[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x8_2 = x8_2.flatten()
            Nk8_ans_2 = np.empty_like(Nk8_FR)
            for i in range(n*2):
                Nk8_ans_2[i] = root_scalar(f, args=(x8_2[i], t), method='toms748', bracket=[-1, 1]).root

            e8_L1 = np.sum(abs(Nk8_ans_2 - Nk8_FR)) * 1 / n
            e8_L2 = np.sqrt(np.sum((Nk8_ans_2 - Nk8_FR)**2) * 1 / n)

        datas = np.load('flying/results_Burger_16_15t_FR2_3001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk16_FRx = data[12][0]
            Nk16_FR2med = 1/2*np.sum(GL.weights*Nk16_FRx, axis=-1)
            Nk16_FR = Nk16_FRx.flatten()
            n = 16
            x16_2 = np.empty((n,2))
            for i in range(len(points)):
                x16_2[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x16_2 = x16_2.flatten()
            Nk16_ans_2 = np.empty_like(Nk16_FR)
            for i in range(n*2):
                Nk16_ans_2[i] = root_scalar(f, args=(x16_2[i], t), method='toms748', bracket=[-1, 1]).root

            e16_L1 = np.sum(abs(Nk16_ans_2 - Nk16_FR)) * 1 / n
            e16_L2 = np.sqrt(np.sum((Nk16_ans_2 - Nk16_FR)**2) * 1 / n)
            R16_L1 = math.log(e8_L1/e16_L1,2)
            R16_L2 = math.log(e8_L2/e16_L2,2)


        datas = np.load('flying/results_Burger_32_15t_FR2_5001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk32_FRx = data[12][0]
            Nk32_FR2med = 1/2*np.sum(GL.weights*Nk32_FRx, axis=-1)
            Nk32_FR = data[12][0].flatten()
            n = 32
            x32_2 = np.empty((n,2))
            for i in range(len(points)):
                x32_2[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x32_2 = x32_2.flatten()
            Nk32_ans_2 = np.empty_like(Nk32_FR)
            for i in range(n*2):
                Nk32_ans_2[i] = root_scalar(f, args=(x32_2[i], t), method='toms748', bracket=[-1, 1]).root

            e32_L1 = np.sum(abs(Nk32_ans_2 - Nk32_FR)) * 1 / n
            e32_L2 = np.sqrt(np.sum((Nk32_ans_2 - Nk32_FR)**2) * 1 / n)
            R32_L1 = math.log(e16_L1/e32_L1,2)
            R32_L2 =  math.log(e16_L2/e32_L2,2)

        datas = np.load('flying/results_Burger_64_15t_FR2_3001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk64_FRx = data[12][0]
            Nk64_FR2med = 1/2*np.sum(GL.weights*Nk64_FRx, axis=-1)
            Nk64_FR = Nk64_FRx.flatten()
            n = 64
            x64_2 = np.empty((n,2))
            for i in range(len(points)):
                x64_2[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x64_2 = x64_2.flatten()
            Nk64_ans_2 = np.empty_like(Nk64_FR)
            for i in range(n*2):
                Nk64_ans_2[i] = root_scalar(f, args=(x64_2[i], t), method='toms748', bracket=[-1, 1]).root

            e64_L1 = np.sum(abs(Nk64_ans_2 - Nk64_FR)) * 1 / n
            e64_L2 = np.sqrt(np.sum((Nk64_ans_2 - Nk64_FR)**2) * 1 / n)
            R64_L1 = math.log(e32_L1/e64_L1,2)
            R64_L2 = math.log(e32_L2/e64_L2,2)


        datas = np.load('flying/results_Burger_128_15t_FR2_3751.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk128_FR = data[12][0].flatten()
            n = 128
            x128_2 = np.empty((n,2))
            for i in range(len(points)):
                x128_2[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x128_2 = x128_2.flatten()
            Nk128_ans_2 = np.empty_like(Nk128_FR)
            for i in range(n*2):
                Nk128_ans_2[i] = root_scalar(f, args=(x128_2[i], t), method='toms748', bracket=[-1, 1]).root

            e128_L1 = np.sum(abs(Nk128_ans_2 - Nk128_FR)) * 1 / n
            e128_L2 = np.sqrt(np.sum((Nk128_ans_2 - Nk128_FR)**2) * 1 / n)
            R128_L1 = math.log(e64_L1/e128_L1,2)
            R128_L2 = math.log(e64_L2/e128_L2,2)


        datas = np.load('flying/results_Burger_256_15t_FR2_15001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk256_FR = data[12][0].flatten()
            n = 256
            x256_2 = np.empty((n,2))
            for i in range(len(points)):
                x256_2[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x256_2 = x256_2.flatten()
            Nk256_ans_2 = np.empty_like(Nk256_FR)
            for i in range(n*2):
                Nk256_ans_2[i] = root_scalar(f, args=(x256_2[i], t), method='toms748', bracket=[-1, 1]).root

            e256_L1 = np.sum(abs(Nk256_ans_2 - Nk256_FR)) * 1 / n
            e256_L2 = np.sqrt(np.sum((Nk256_ans_2 - Nk256_FR)**2) * 1 / n)
            R256_L1 = math.log(e128_L1/e256_L1,2)
            R256_L2 = math.log(e128_L2/e256_L2,2)

        datas = np.load('flying/results_Burger_512_15t_FR2_7501.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk512_FR = data[12][0].flatten()
            n = 512
            x512_2 = np.empty((n,2))
            for i in range(len(points)):
                x512_2[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x512_2 = x512_2.flatten()
            Nk512_ans_2 = np.empty_like(Nk512_FR)
            for i in range(n*2):
                Nk512_ans_2[i] = root_scalar(f, args=(x512_2[i], t), method='toms748', bracket=[-1, 1]).root

            e512_L1 = np.sum(abs(Nk512_ans_2 - Nk512_FR)) * 1 / n
            e512_L2 = np.sqrt(np.sum((Nk512_ans_2 - Nk512_FR)**2) * 1 / n)
            R512_L1 = (np.log(e256_L1) - np.log(e512_L1))/(np.log(2/256) - np.log(2/512))
            R512_L2 = math.log(e256_L2/e512_L2,2)

        datas = np.load('flying/results_Burger_80_FR2_1501.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk80_FR2med = 1/2*np.sum(GL.weights*data[12][0], axis=-1)


        '---------------------------Flux Reconstruction------------------------'
        datas = np.load('flying/results_Burger_8_15t_FR3_1501.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk8_FR3 = data[12][0].flatten()
            n = 8
            x8_3 = np.empty((n,3))
            GL = quadpy.c1.gauss_lobatto(3)
            points = GL.points
            for i in range(len(points)):
                x8_3[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x8_3 = x8_3.flatten()
            Nk8_ans_3 = np.empty_like(Nk8_FR3)
            for i in range(n*3):
                Nk8_ans_3[i] = root_scalar(f, args=(x8_3[i], t), method='toms748', bracket=[-1, 1]).root

            e8_L1_3 = np.sum(abs(Nk8_ans_3 - Nk8_FR3)) * 1 / n

        datas = np.load('flying/results_Burger_16_15t_FR3_3001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk16_FR3 = data[12][0].flatten()
            n = 16
            x16_3 = np.empty((n,3))
            for i in range(len(points)):
                x16_3[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x16_3 = x16_3.flatten()
            Nk16_ans_3 = np.empty_like(Nk16_FR3)
            for i in range(n*3):
                Nk16_ans_3[i] = root_scalar(f, args=(x16_3[i], t), method='toms748', bracket=[-1, 1]).root

            e16_L1_3 = np.sum(abs(Nk16_ans_3 - Nk16_FR3)) * 1 / n
            R16_3 = math.log(e8_L1_3/e16_L1_3,2)

        datas = np.load('flying/results_Burger_32_15t_FR3_5001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk32_FRx = data[12][0]
            Nk32_FR3med = 1/2*np.sum(GL.weights*Nk32_FRx, axis=-1)
            Nk32_FR3 = data[12][0].flatten()
            n = 32
            x32_3 = np.empty((n,3))
            for i in range(len(points)):
                x32_3[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x32_3 = x32_3.flatten()
            Nk32_ans_3 = np.empty_like(Nk32_FR3)
            for i in range(n*3):
                Nk32_ans_3[i] = root_scalar(f, args=(x32_3[i], t), method='toms748', bracket=[-1, 1]).root

            e32_L1_3 = np.sum(abs(Nk32_ans_3 - Nk32_FR3)) * 1 / n
            R32_3 = math.log(e16_L1_3/e32_L1_3,2)

        datas = np.load('flying/results_Burger_64_15t_FR3_3001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk64_FRx = data[12][0]
            Nk64_FR3med = 1/2*np.sum(GL.weights*Nk64_FRx, axis=-1)
            Nk64_FR3 = data[12][0].flatten()
            n = 64
            x64_3 = np.empty((n,3))
            for i in range(len(points)):
                x64_3[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x64_3 = x64_3.flatten()
            Nk64_ans_3 = np.empty_like(Nk64_FR3)
            for i in range(n*3):
                Nk64_ans_3[i] = root_scalar(f, args=(x64_3[i], t), method='toms748', bracket=[-1, 1]).root

            e64_L1_3 = np.sum(abs(Nk64_ans_3 - Nk64_FR3)) * 1 / n
            R64_3 = math.log(e32_L1_3/e64_L1_3,2)

        datas = np.load('flying/results_Burger_128_15t_FR3_3751.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk128_FR3 = data[12][0].flatten()
            n = 128
            x128_3 = np.empty((n,3))
            for i in range(len(points)):
                x128_3[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x128_3 = x128_3.flatten()
            Nk128_ans_3 = np.empty_like(Nk128_FR3)
            for i in range(n*3):
                Nk128_ans_3[i] = root_scalar(f, args=(x128_3[i], t), method='toms748', bracket=[-1, 1]).root

            e128_L1_3 = np.sum(abs(Nk128_ans_3 - Nk128_FR3)) * 1 / n
            R128_3 = math.log(e64_L1_3/e128_L1_3,2)

        datas = np.load('flying/results_Burger_256_15t_FR3_15001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk256_FR3 = data[12][0].flatten()
            n = 256
            x256_3 = np.empty((n,3))
            for i in range(len(points)):
                x256_3[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x256_3 = x256_3.flatten()
            Nk256_ans_3 = np.empty_like(Nk256_FR3)
            for i in range(n*3):
                Nk256_ans_3[i] = root_scalar(f, args=(x256_3[i], t), method='toms748', bracket=[-1, 1]).root

            e256_L1_3 = np.sum(abs(Nk256_ans_3 - Nk256_FR3)) * 1 / n
            R256_3 = math.log(e128_L1_3/e256_L1_3,2)

        datas = np.load('flying/results_Burger_512_15t_FR3_7501.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk512_FR3 = data[12][0].flatten()
            n = 512
            x512_3 = np.empty((n,3))
            for i in range(len(points)):
                x512_3[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x512_3 = x512_3.flatten()
            Nk512_ans_3 = np.empty_like(Nk512_FR3)
            for i in range(n*3):
                Nk512_ans_3[i] = root_scalar(f, args=(x512_3[i], t), method='toms748', bracket=[-1, 1]).root

            e512_L1_3 = np.sum(abs(Nk512_ans_3 - Nk512_FR3)) * 1 / n
            R512_3 = (np.log(e256_L1_3) - np.log(e512_L1_3))/(np.log(2/256) - np.log(2/512))

        datas = np.load('flying/results_Burger_80_FR3_3001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk80_FR3med = 1/2*np.sum(GL.weights*data[12][0], axis=-1)


        '''----------------------------FR4-------------------------'''

        datas = np.load('flying/results_Burger_8_15t_FR4_1501.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk8_FR4 = data[12][0].flatten()
            n = 8
            x8_4 = np.empty((n,4))
            GL = quadpy.c1.gauss_lobatto(4)
            points = GL.points
            for i in range(len(points)):
                x8_4[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x8_4 = x8_4.flatten()
            Nk8_ans_4 = np.empty_like(Nk8_FR4)
            for i in range(n*4):
                Nk8_ans_4[i] = root_scalar(f, args=(x8_4[i], t), method='toms748', bracket=[-1, 1]).root

            e8_L1_4 = np.sum(abs(Nk8_ans_4 - Nk8_FR4)) * 1 / n

        datas = np.load('flying/results_Burger_16_15t_FR4_3001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk16_FR4 = data[12][0].flatten()
            n = 16
            x16_4 = np.empty((n,4))
            for i in range(len(points)):
                x16_4[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x16_4 = x16_4.flatten()
            Nk16_ans_4 = np.empty_like(Nk16_FR4)
            for i in range(n*4):
                Nk16_ans_4[i] = root_scalar(f, args=(x16_4[i], t), method='toms748', bracket=[-1, 1]).root

            e16_L1_4 = np.sum(abs(Nk16_ans_4 - Nk16_FR4)) * 1 / n
            R16_4 = math.log(e8_L1_4/e16_L1_4, 2)

        datas = np.load('flying/results_Burger_32_15t_FR4_5001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk32_FRx = data[12][0]
            Nk32_FR4med = 1/2*np.sum(GL.weights*Nk32_FRx, axis=-1)

            Nk32_FR4 = data[12][0].flatten()
            n = 32
            x32_4 = np.empty((n,4))
            for i in range(len(points)):
                x32_4[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x32_4 = x32_4.flatten()
            Nk32_ans_4 = np.empty_like(Nk32_FR4)
            for i in range(n*4):
                Nk32_ans_4[i] = root_scalar(f, args=(x32_4[i], t), method='toms748', bracket=[-1, 1]).root

            e32_L1_4 = np.sum(abs(Nk32_ans_4 - Nk32_FR4)) * 1 / n
            R32_4 = math.log(e16_L1_4/e32_L1_4,2)

        datas = np.load('flying/results_Burger_64_15t_FR4_3001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk64_FRx = data[12][0]
            Nk64_FR4med = 1/2*np.sum(GL.weights*Nk64_FRx, axis=-1)
            Nk64_FR4 = data[12][0].flatten()
            n = 64
            x64_4 = np.empty((n,4))
            for i in range(len(points)):
                x64_4[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x64_4 = x64_4.flatten()
            Nk64_ans_4 = np.empty_like(Nk64_FR4)
            for i in range(n*4):
                Nk64_ans_4[i] = root_scalar(f, args=(x64_4[i], t), method='toms748', bracket=[-1, 1]).root

            e64_L1_4 = np.sum(abs(Nk64_ans_4 - Nk64_FR4)) * 1 / n
            R64_4 = math.log(e32_L1_4/e64_L1_4,2)

        datas = np.load('flying/results_Burger_128_15t_FR4_3751.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk128_FR4 = data[12][0].flatten()
            n = 128
            x128_4 = np.empty((n,4))
            for i in range(len(points)):
                x128_4[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x128_4 = x128_4.flatten()
            Nk128_ans_4 = np.empty_like(Nk128_FR4)
            for i in range(n*4):
                Nk128_ans_4[i] = root_scalar(f, args=(x128_4[i], t), method='toms748', bracket=[-1, 1]).root

            e128_L1_4 = np.sum(abs(Nk128_ans_4 - Nk128_FR4)) * 1 / n
            R128_4 = math.log(e64_L1_4/e128_L1_4,2)

        datas = np.load('flying/results_Burger_256_15t_FR4_15001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk256_FR4 = data[12][0].flatten()
            n = 256
            x256_4 = np.empty((n,4))
            for i in range(len(points)):
                x256_4[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x256_4 = x256_4.flatten()
            Nk256_ans_4 = np.empty_like(Nk256_FR4)
            for i in range(n*4):
                Nk256_ans_4[i] = root_scalar(f, args=(x256_4[i], t), method='toms748', bracket=[-1, 1]).root

            e256_L1_4 = np.sum(abs(Nk256_ans_4 - Nk256_FR4)) * 1 / n
            R256_4 = math.log(e128_L1_4/e256_L1_4,2)

        datas = np.load('flying/results_Burger_512_15t_FR4_7501.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk512_FR4 = data[12][0].flatten()
            n = 512
            x512_4 = np.empty((n,4))
            for i in range(len(points)):
                x512_4[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x512_4 = x512_4.flatten()
            Nk512_ans_4 = np.empty_like(Nk512_FR4)
            for i in range(n*4):
                Nk512_ans_4[i] = root_scalar(f, args=(x512_4[i], t), method='toms748', bracket=[-1, 1]).root

            e512_L1_4 = np.sum(abs(Nk512_ans_4 - Nk512_FR4)) * 1 / n
            R512_4 = math.log(e256_L1_4/e512_L1_4,2)

        datas = np.load('flying/results_Burger_80_FR4_5001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk80_FR4med = 1/2*np.sum(GL.weights*data[12][0], axis=-1)


        plt.figure(1)
        plt.plot(x8_2, Nk8_FR, '-b^', x16_2, Nk16_FR, '-yo', x32_2, Nk32_FR, '-rv',
                x64_2, Nk64_FR, '-gs', x128_2, Nk128_FR, '-m<', x256_2, Nk256_FR, '-cP')
        plt.plot(x512_2, Nk512_FR, 'tab:pink',marker='*')
        plt.plot(x, Nk_ans, 'k-')
        plt.grid()
        plt.xlabel('x')
        plt.ylabel('Nk [mole]')
        plt.legend(('FR - 8 CV', 'FR - 16 CV', 'FR - 32 CV', 'FR - 64 CV',
                    'FR - 128 CV', 'FR - 256 CV', 'FR - 512 CV', 'Reference solution'))
        plt.title('Results for 2nd order approximation')
        plt.savefig('results/compositional/FR_paper/Burger/Nk_Burgers_comparison_15t_order2.eps', format='eps')

        plt.figure(2)
        plt.plot(x8_3, Nk8_FR3, '-b^', x16_3, Nk16_FR3, '-yo', x32_3, Nk32_FR3, '-rv',
                x64_3, Nk64_FR3, '-gs', x128_3, Nk128_FR3, '-m<', x256_3, Nk256_FR3, '-cP')
        plt.plot(x512_3, Nk512_FR3, 'tab:pink',marker='*')
        plt.plot(x, Nk_ans, 'k-')
        plt.grid()
        plt.xlabel('x')
        plt.ylabel('Nk [mole]')
        plt.legend(('FR - 8 CV', 'FR - 16 CV', 'FR - 32 CV', 'FR - 64 CV',
                    'FR - 128 CV', 'FR - 256 CV', 'FR - 512 CV', 'Reference solution'))
        plt.title('Results for 3rd order approximation')
        plt.savefig('results/compositional/FR_paper/Burger/Nk_Burgers_comparison_15t_order3.eps', format='eps')

        plt.figure(3)
        plt.plot(x8_4, Nk8_FR4, '-b^', x16_4, Nk16_FR4, '-yo', x32_4, Nk32_FR4, '-rv',
                x64_4, Nk64_FR4, '-gs', x128_4, Nk128_FR4, '-m<', x256_4, Nk256_FR4, '-cP')
        plt.plot(x512_4, Nk512_FR4, 'tab:pink',marker='*')
        plt.plot(x, Nk_ans, 'k-')
        plt.grid()
        plt.xlabel('x')
        plt.ylabel('Nk [mole]')
        plt.legend(('FR - 8 CV', 'FR - 16 CV', 'FR - 32 CV', 'FR - 64 CV',
                    'FR - 128 CV', 'FR - 256 CV', 'FR - 512 CV', 'Reference solution'))
        plt.title('Results for 4th order approximation')
        plt.savefig('results/compositional/FR_paper/Burger/Nk_Burgers_comparison_15t_order4.eps', format='eps')

        plt.figure(4)
        x80 = np.linspace(0+1/160, 1-1/160, 80)
        plt.plot(x80, Nk80_FR2med, '-bo')
        plt.plot(x, Nk_ans, 'k-')
        plt.grid()
        plt.xlabel('x')
        plt.ylabel('Nk [mole]')
        plt.legend(('2nd order', 'Reference solution'))
        plt.title('Results for 80 control volumes mesh')
        plt.savefig('results/compositional/FR_paper/Burger/Nk_Burgers_15t_80CV_2.eps', format='eps')

        plt.figure(5)
        plt.plot(x80, Nk80_FR3med, '-bo')
        plt.plot(x, Nk_ans, 'k-')
        plt.grid()
        plt.xlabel('x')
        plt.ylabel('Nk [mole]')
        plt.legend(('3rd order', 'Reference solution'))
        plt.title('Results for 80 control volumes mesh')
        plt.savefig('results/compositional/FR_paper/Burger/Nk_Burgers_15t_80CV_3.eps', format='eps')

        plt.figure(7)
        plt.plot(x80, Nk80_FR4med, 'bo', mfc='none')
        plt.plot(x, Nk_ans, 'k-')
        plt.grid()
        plt.xlabel('Distância')
        plt.ylabel('Nk [mol]')
        plt.legend(('CPR-4$^a$ ordem', 'Semi-Analítica'))
        plt.title('Resultados para malha 80x1x1')
        plt.savefig('results/compositional/FR_paper/Burger/Nk_Burgers_15t_80CV_4.png')

        plt.figure(8)
        x64 = np.linspace(0+1/128, 1-1/128, 64)
        plt.plot(x64[16:32], Nk64_FR4med[16:32], '-bo', x64[16:32], Nk64_FR3med[16:32], '-rv', x64[16:32], Nk64_FR2med[16:32], '-g<')
        plt.plot(x[500:1000], Nk_ans[500:1000], 'k-')
        plt.grid()
        plt.xlabel('Distância')
        plt.ylabel('Nk [mol]')
        plt.legend(('4th order', '3rd order', '2nd order', 'Reference solution'))
        plt.title('Results for 64 control volumes mesh')
        plt.savefig('results/compositional/FR_paper/Burger/Nk_Burgers_comparison_15t_64CV.eps', format='eps')

        plt.figure(9)
        x32 = np.linspace(0+1/64, 1-1/64, 32)
        plt.plot(x32, Nk32_FR2med, '-bs')
        plt.plot(x, Nk_ans, 'k-')
        plt.grid()
        plt.xlabel('x')
        plt.ylabel('Nk [mole]')
        plt.legend(('2nd order', 'Reference solution'))
        plt.title('Results for 32 control volumes mesh')
        plt.savefig('results/compositional/FR/Nk_Burgers_FR2_15t_32CV.png')

        plt.figure(11)
        x16 = np.linspace(0+1/32, 1-1/32, 16)
        plt.plot(x16, Nk16_FR2med, '-bs')
        plt.plot(x, Nk_ans, 'k-')
        plt.grid()
        plt.xlabel('x')
        plt.ylabel('Nk [mole]')
        plt.legend(('2nd order', 'Reference solution'))
        plt.title('Results for 16 control volumes mesh')
        plt.savefig('results/compositional/FR/Nk_Burgers_FR2_15t_16CV.png')

        plt.figure(10)
        x80 = np.linspace(0+1/160, 1-1/160, 80)
        plt.plot(x80[32:39], Nk80_FR4med[32:39], 'bo', x80[32:39], Nk80_FR3med[32:39],
            'rv', x80[32:39], Nk80_FR2med[32:39], 'g<', mfc='none')
        plt.plot(x[800:999], Nk_ans[800:999], 'k-')
        plt.grid()
        plt.xlabel('Distância')
        plt.ylabel('Nk [mol]')
        plt.legend(('CPR-4$^a$ ordem', 'CPR-3$^a$ ordem', 'CPR-2$^a$ ordem', 'Semi-Analítica'))
        plt.title('Resultados para malha 80x1x1')
        plt.savefig('results/compositional/FR_paper/Burger/Nk_Burgers_comparison_15t_80CV.png')

        import pdb; pdb.set_trace()
        plt.figure(6)
        x = np.log10(np.array([8,16,32,64,128,256,512]))
        y = np.log10(np.array([e8_L1, e16_L1, e32_L1, e64_L1, e128_L1, e256_L1, e512_L1]))
        y_FR3 = np.log10(np.array([e8_L1_3, e16_L1_3, e32_L1_3, e64_L1_3, e128_L1_3, e256_L1_3, e512_L1_3]))
        y_FR4 = np.log10(np.array([e8_L1_4, e16_L1_4, e32_L1_4, e64_L1_4, e128_L1_4, e256_L1_4, e512_L1_4]))
        #y_FR5 = np.log10(np.array([e8_L1_5, e16_L1_5, e32_L1_5, e64_L1_5, e128_L1_5, e256_L1_5]))
        #y_MUSCL = np.log10(np.array([e32_L1_MUSCL, e64_L1_MUSCL, e128_L1_MUSCL, e256_L1_MUSCL, e512_L1_MUSCL]))
        y_ref = -2*x +2
        plt.plot(x, y, 'r', x, y_FR3, 'g', x, y_FR4, 'y')
        plt.grid()
        plt.ylabel('$log_{10}({E}_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.legend(('FR-2nd order', 'FR-3rd order', 'FR-4th order', 'FR-5th order'))
        plt.savefig('results/compositional/FR_paper/Burger/Nk_Burgers_convergence_15t.eps', format='eps')
        import pdb; pdb.set_trace()
