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
t = 0.3

def f(x, pos, t):
    return x - (1/(2*math.pi))*np.sin(2*math.pi*(pos-t*x))

datas = np.load('flying/results_Burger_10000_upw_2060.npy', allow_pickle=True)

for data in datas[1:]:
    Nk_FOU = data[12][0]

for arq in arquivos:
    if  arq.startswith(name):


        '---------------------------Flux Reconstruction------------------------'
        datas = np.load('flying/results_Burger_8_FR2_600.npy', allow_pickle=True)

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

        datas = np.load('flying/results_Burger_16_FR2_600.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk16_FR = data[12][0].flatten()
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


        datas = np.load('flying/results_Burger_32_FR2_1000.npy', allow_pickle=True)

        for data in datas[1:]:
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
            Nk32_FR_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        datas = np.load('flying/results_Burger_64_FR2_600.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk64_FR = data[12][0].flatten()
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
            Nk64_FR_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2


        datas = np.load('flying/results_Burger_128_FR2_750.npy', allow_pickle=True)

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


        datas = np.load('flying/results_Burger_256_FR2_3001.npy', allow_pickle=True)

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

        datas = np.load('flying/results_Burger_512_FR2_1501.npy', allow_pickle=True)

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
            import pdb; pdb.set_trace()

        '---------------------------Flux Reconstruction------------------------'
        #datas = np.load('flying/results_Burger_8_FR3_600.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_8_FR3t_300.npy', allow_pickle=True)
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

        #datas = np.load('flying/results_Burger_16_FR3_1000.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_16_FR3t_600.npy', allow_pickle=True)

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
            import pdb; pdb.set_trace()

        datas = np.load('flying/results_Burger_32_FR3_3001.npy', allow_pickle=True)

        for data in datas[1:]:
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
            Nk32_FR3_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        datas = np.load('flying/results_Burger_64_FR3_600.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk64_FR3 = data[12][0].flatten()
            n = 64
            x64_3 = np.empty((n,3))
            for i in range(len(points)):
                x64_3[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x64_3 = x64_3.flatten()
            Nk64_ans_3 = np.empty_like(Nk64_FR3)
            for i in range(n*3):
                Nk64_ans_3[i] = root_scalar(f, args=(x64_3[i], t), method='toms748', bracket=[-1, 1]).root

            e64_L1_3 = np.sum(abs(Nk64_ans_3 - Nk64_FR3)) * 2 / n
            R64_3 = math.log(e32_L1_3/e64_L1_3,2)
            Nk64_FR3_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        datas = np.load('flying/results_Burger_128_FR3_750.npy', allow_pickle=True)

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

        datas = np.load('flying/results_Burger_256_FR3_3001.npy', allow_pickle=True)

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

        datas = np.load('flying/results_Burger_512_FR3_1501.npy', allow_pickle=True)

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

        '''----------------------------FR4-------------------------'''

        datas = np.load('flying/results_Burger_8_FR4_1000.npy', allow_pickle=True)

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

            e8_L1_4 = np.sum(abs(Nk8_ans_4 - Nk8_FR4)) * 2 / n

        datas = np.load('flying/results_Burger_16_FR4_3001.npy', allow_pickle=True)

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

            e16_L1_4 = np.sum(abs(Nk16_ans_4 - Nk16_FR4)) * 2 / n
            R16_4 = math.log(e8_L1_4/e16_L1_4, 2)

        datas = np.load('flying/results_Burger_32_FR4_3001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk32_FR4 = data[12][0].flatten()
            n = 32
            x32_4 = np.empty((n,4))
            for i in range(len(points)):
                x32_4[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x32_4 = x32_4.flatten()
            Nk32_ans_4 = np.empty_like(Nk32_FR4)
            for i in range(n*4):
                Nk32_ans_4[i] = root_scalar(f, args=(x32_4[i], t), method='toms748', bracket=[-1, 1]).root

            e32_L1_4 = np.sum(abs(Nk32_ans_4 - Nk32_FR4)) * 2 / n
            R32_4 = math.log(e16_L1_4/e32_L1_4,2)
            Nk32_FR4_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        datas = np.load('flying/results_Burger_64_FR4_600.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk64_FR4 = data[12][0].flatten()
            n = 64
            x64_4 = np.empty((n,4))
            for i in range(len(points)):
                x64_4[:,i] = np.linspace(0+(points[i] + 1)*1/(2*n),1-(1-points[i])*1/(2*n),n)
            x64_4 = x64_4.flatten()
            Nk64_ans_4 = np.empty_like(Nk64_FR4)
            for i in range(n*4):
                Nk64_ans_4[i] = root_scalar(f, args=(x64_4[i], t), method='toms748', bracket=[-1, 1]).root

            e64_L1_4 = np.sum(abs(Nk64_ans_4 - Nk64_FR4)) * 2 / n
            R64_4 = math.log(e32_L1_4/e64_L1_4,2)
            Nk64_FR4_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        datas = np.load('flying/results_Burger_128_FR4_750.npy', allow_pickle=True)

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

            e128_L1_4 = np.sum(abs(Nk128_ans_4 - Nk128_FR4)) * 2 / n
            R128_4 = math.log(e64_L1_4/e128_L1_4,2)

        datas = np.load('flying/results_Burger_256_FR4_3001.npy', allow_pickle=True)

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

            e256_L1_4 = np.sum(abs(Nk256_ans_4 - Nk256_FR4)) * 2 / n
            R256_4 = math.log(e128_L1_4/e256_L1_4,2)

        datas = np.load('flying/results_Burger_512_FR4_1501.npy', allow_pickle=True)

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

            e512_L1_4 = np.sum(abs(Nk512_ans_4 - Nk512_FR4)) * 2 / n
            R512_4 = math.log(e256_L1_4/e512_L1_4,2)


        plt.figure(1)
        x = np.linspace(0+1/20000, 1-1/20000, 10000)
        plt.plot(x, Nk_FOU, 'k', x512_4, Nk512_ans_4, 'b')
        plt.legend(('FOU', 'Analytical'))
        plt.savefig('results/compositional/FR/Nk_Burgers_sol_03t')

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
        plt.savefig('results/compositional/FR/Nk_Burgers_convergence.png')

        plt.figure(7)
        x32 = np.linspace(0+1/64,1-1/64,32)
        plt.plot(x32, Nk32_FR_avg, 'r^', x32, Nk32_FR3_avg, 'go', x32, Nk32_FR4_avg, 'ys', mfc='none')
        plt.plot(x512_4, Nk512_ans_4, 'k')
        plt.grid()
        plt.ylabel('$N_k$')
        plt.xlabel('Distância')
        plt.title('Resultados para t=0.3 com malha 32x1x1')
        plt.legend(('CPR-$2^a$ ordem', 'CPR-$3^a$ order', 'CPR-$4^a$ ordem', 'Semi-Analítica'))
        plt.savefig('results/compositional/FR/Nk_Burgers_32.png')
        import pdb; pdb.set_trace()
