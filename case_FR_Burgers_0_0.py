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
t = 0.5

def f(x, pos, t):
    return x - (1/(2*math.pi))*np.sin(2*math.pi*(pos-t*x))


datas = np.load('flying/results_Burger_10000_upw_2060.npy', allow_pickle=True)

for data in datas[1:]:
    Nk_FOU = data[12][0]

for arq in arquivos:
    if  arq.startswith(name):


        '---------------------------Flux Reconstruction------------------------'
        #datas = np.load('flying/results_Burger_8_FR2_600.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_8_FR2_RK3_nolim_3.npy', allow_pickle=True)
        for data in datas[1:]:
            t8_FR2 = data[2]
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
            Nk8_FR_avg = np.sum(data[12][0] * GL.weights,axis=-1)/np.sum(GL.weights)


        #datas = np.load('flying/results_Burger_16_FR2_600.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_16_FR2_RK3_nolim_4.npy', allow_pickle=True)

        for data in datas[1:]:
            t16_FR2 = data[2]
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
            Nk16_FR_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2



        #datas = np.load('flying/results_Burger_32_FR2_1000.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_32_FR2_RK3_nolim_6.npy', allow_pickle=True)

        for data in datas[1:]:
            t32_FR2 = data[2]
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

        #datas = np.load('flying/results_Burger_64_FR2_600.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_64_FR2_RK3_nolim_11.npy', allow_pickle=True)

        for data in datas[1:]:
            t64_FR2 = data[2]
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


        #datas = np.load('flying/results_Burger_128_FR2_750.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_128_FR2_RK3_nolim_21.npy', allow_pickle=True)

        for data in datas[1:]:
            t128_FR2 = data[2]
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
            Nk128_FR_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2


        #datas = np.load('flying/results_Burger_256_FR2_3001.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_256_FR2_RK3_nolim_42.npy', allow_pickle=True)

        for data in datas[1:]:
            t256_FR2 = data[2]
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

        #datas = np.load('flying/results_Burger_512_FR2_1501.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_512_FR2_RK3_nolim_80.npy', allow_pickle=True)

        for data in datas[1:]:
            t512_FR2 = data[2]
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
            Nk512_FR_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2


        '---------------------------Flux Reconstruction------------------------'
        #datas = np.load('flying/results_Burger_8_FR3_600.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_8_FR3_RK3_nolim_4.npy', allow_pickle=True)
        for data in datas[1:]:
            t8_FR3 = data[2]
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
            Nk8_FR3_avg = np.sum(data[12][0] * GL.weights,axis=-1)/np.sum(GL.weights)

        #datas = np.load('flying/results_Burger_16_FR3_1000.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_16_FR3_RK3_nolim_6.npy', allow_pickle=True)

        for data in datas[1:]:
            t16_FR3 = data[2]
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
            Nk16_FR3_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        #datas = np.load('flying/results_Burger_32_FR3_3001.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_32_FR3_RK3_nolim_10.npy', allow_pickle=True)

        for data in datas[1:]:
            t32_FR3 = data[2]
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

        #datas = np.load('flying/results_Burger_64_FR3_600.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_64_FR3_RK3_nolim_18.npy', allow_pickle=True)

        for data in datas[1:]:
            t64_FR3 = data[2]
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
            Nk64_FR3_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        #datas = np.load('flying/results_Burger_128_FR3_750.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_128_FR3_RK3_nolim_34.npy', allow_pickle=True)

        for data in datas[1:]:
            t128_FR3 = data[2]
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
            Nk128_FR3_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        #datas = np.load('flying/results_Burger_256_FR3_3001.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_256_FR3_RK3_nolim_68.npy', allow_pickle=True)

        for data in datas[1:]:
            t256_FR3 = data[2]
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

        #datas = np.load('flying/results_Burger_512_FR3_1501.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_512_FR3_RK3_nolim_133.npy', allow_pickle=True)

        for data in datas[1:]:
            t512_FR3 = data[2]
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
            Nk512_FR3_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2
        '''----------------------------FR4-------------------------'''

        #datas = np.load('flying/results_Burger_8_FR4_1000.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_8_FR4_RK3_nolim_4.npy', allow_pickle=True)
        #datas = np.load('flying/results_Burger_8_FR4_RK3_nolim_4.npy', allow_pickle=True)
        for data in datas[1:]:
            t8_FR4 = data[2]
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
            Nk8_FR4_avg = np.sum(data[12][0] * GL.weights,axis=-1)/np.sum(GL.weights)

        #datas = np.load('flying/results_Burger_16_FR4_3001.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_16_FR4_RK3_nolim_7.npy', allow_pickle=True)
        for data in datas[1:]:
            t16_FR4 = data[2]
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
            Nk16_FR4_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        #datas = np.load('flying/results_Burger_32_FR4_3001.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_32_FR4_RK3_nolim_13.npy', allow_pickle=True)

        for data in datas[1:]:
            t32_FR4 = data[2]
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
            Nk32_FR4_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        #datas = np.load('flying/results_Burger_64_FR4_600.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_64_FR4_RK3_nolim_24.npy', allow_pickle=True)

        for data in datas[1:]:
            t64_FR4 = data[2]
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
            Nk64_FR4_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        #datas = np.load('flying/results_Burger_128_FR4_750.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_128_FR4_RK3_nolim_47.npy', allow_pickle=True)

        for data in datas[1:]:
            t128_FR4 = data[2]
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
            Nk128_FR4_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        #datas = np.load('flying/results_Burger_256_FR4_3001.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_256_FR4_RK3_nolim_95.npy', allow_pickle=True)

        for data in datas[1:]:
            t256_FR4 = data[2]
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

        #datas = np.load('flying/results_Burger_512_FR4_1501.npy', allow_pickle=True)
        datas = np.load('flying/results_Burger_512_FR4_RK3_nolim_191.npy', allow_pickle=True)

        for data in datas[1:]:
            t512_FR4 = data[2]
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
            Nk512_FR4_avg = np.sum(data[12][0]*GL.weights,axis=-1)/2

        'FOU'

        datas = np.load('flying/results_Burger_8_FOU_ROE_RK3_2.npy', allow_pickle=True)
        for data in datas[1:]:
            t8_FOU = data[2]
            Nk8_FOU = data[12][0].flatten()
            n = 8
            x8 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk8_ans = np.empty_like(Nk8_FOU)
            for i in range(n):
                Nk8_ans[i] = root_scalar(f, args=(x8[i], t), method='toms748', bracket=[-1, 1]).root
            e8_L1_FOU = np.sum(abs(Nk8_ans - Nk8_FOU)) * 1 / n
            import pdb; pdb.set_trace()

        datas = np.load('flying/results_Burger_16_FOU_ROE_RK3_2.npy', allow_pickle=True)
        for data in datas[1:]:
            t16_FOU = data[2]
            Nk16_FOU = data[12][0].flatten()
            n = 16
            x16 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk16_ans = np.empty_like(Nk16_FOU)
            for i in range(n):
                Nk16_ans[i] = root_scalar(f, args=(x16[i], t), method='toms748', bracket=[-1, 1]).root
            e16_L1_FOU = np.sum(abs(Nk16_ans - Nk16_FOU)) * 1 / n

        datas = np.load('flying/results_Burger_32_FOU_ROE_RK3_3.npy', allow_pickle=True)
        for data in datas[1:]:
            t32_FOU = data[2]
            Nk32_FOU = data[12][0].flatten()
            n = 32
            x32 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk32_ans = np.empty_like(Nk32_FOU)
            for i in range(n):
                Nk32_ans[i] = root_scalar(f, args=(x32[i], t), method='toms748', bracket=[-1, 1]).root
            e32_L1_FOU = np.sum(abs(Nk32_ans - Nk32_FOU)) * 1 / n

        datas = np.load('flying/results_Burger_64_FOU_ROE_RK3_5.npy', allow_pickle=True)
        for data in datas[1:]:
            t64_FOU = data[2]
            Nk64_FOU = data[12][0].flatten()
            n = 64
            x64 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk64_ans = np.empty_like(Nk64_FOU)
            for i in range(n):
                Nk64_ans[i] = root_scalar(f, args=(x64[i], t), method='toms748', bracket=[-1, 1]).root
            e64_L1_FOU = np.sum(abs(Nk64_ans - Nk64_FOU)) * 1 / n

        datas = np.load('flying/results_Burger_128_FOU_ROE_RK3_8.npy', allow_pickle=True)
        for data in datas[1:]:
            t128_FOU = data[2]
            Nk128_FOU = data[12][0].flatten()
            n = 128
            x128 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk128_ans = np.empty_like(Nk128_FOU)
            for i in range(n):
                Nk128_ans[i] = root_scalar(f, args=(x128[i], t), method='toms748', bracket=[-1, 1]).root
            e128_L1_FOU = np.sum(abs(Nk128_ans - Nk128_FOU)) * 1 / n


        datas = np.load('flying/results_Burger_256_FOU_ROE_RK3_15.npy', allow_pickle=True)
        for data in datas[1:]:
            t256_FOU = data[2]
            Nk256_FOU = data[12][0].flatten()
            n = 256
            x256 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk256_ans = np.empty_like(Nk256_FOU)
            for i in range(n):
                Nk256_ans[i] = root_scalar(f, args=(x256[i], t), method='toms748', bracket=[-1, 1]).root
            e256_L1_FOU = np.sum(abs(Nk256_ans - Nk256_FOU)) * 1 / n

        datas = np.load('flying/results_Burger_512_FOU_ROE_RK3_28.npy', allow_pickle=True)
        for data in datas[1:]:
            t512_FOU = data[2]
            Nk512_FOU = data[12][0].flatten()
            n = 512
            x512 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk512_ans = np.empty_like(Nk512_FOU)
            for i in range(n):
                Nk512_ans[i] = root_scalar(f, args=(x512[i], t), method='toms748', bracket=[-1, 1]).root
            e512_L1_FOU = np.sum(abs(Nk512_ans - Nk512_FOU)) * 1 / n

        datas = np.load('flying/results_Burger_1024_FOU_ROE_RK3_56.npy', allow_pickle=True)
        for data in datas[1:]:
            t1024_FOU = data[2]
            Nk1024_FOU = data[12][0].flatten()
            n = 1024
            x1024 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk1024_ans = np.empty_like(Nk1024_FOU)
            for i in range(n):
                Nk1024_ans[i] = root_scalar(f, args=(x1024[i], t), method='toms748', bracket=[-1, 1]).root
            e1024_L1_FOU = np.sum(abs(Nk1024_ans - Nk1024_FOU)) * 1 / n

        datas = np.load('flying/results_Burger_2048_FOU_ROE_RK3_110.npy', allow_pickle=True)
        for data in datas[1:]:
            t2048_FOU = data[2]
            Nk2048_FOU = data[12][0].flatten()
            n = 2048
            x2048 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk2048_ans = np.empty_like(Nk2048_FOU)
            for i in range(n):
                Nk2048_ans[i] = root_scalar(f, args=(x2048[i], t), method='toms748', bracket=[-1, 1]).root
            e2048_L1_FOU = np.sum(abs(Nk2048_ans - Nk2048_FOU)) * 1 / n

        datas = np.load('flying/results_Burger_4096_FOU_ROE_RK3_218.npy', allow_pickle=True)
        for data in datas[1:]:
            t4096_FOU = data[2]
            Nk4096_FOU = data[12][0].flatten()
            n = 4096
            x4096 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk4096_ans = np.empty_like(Nk4096_FOU)
            for i in range(n):
                Nk4096_ans[i] = root_scalar(f, args=(x4096[i], t), method='toms748', bracket=[-1, 1]).root
            e4096_L1_FOU = np.sum(abs(Nk4096_ans - Nk4096_FOU)) * 1 / n

        datas = np.load('flying/results_Burger_8192_FOU_ROE_RK3_435.npy', allow_pickle=True)
        for data in datas[1:]:
            t8192_FOU = data[2]
            Nk8192_FOU = data[12][0].flatten()
            n = 8192
            x8192 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk8192_ans = np.empty_like(Nk8192_FOU)
            for i in range(n):
                Nk8192_ans[i] = root_scalar(f, args=(x8192[i], t), method='toms748', bracket=[-1, 1]).root
            e8192_L1_FOU = np.sum(abs(Nk8192_ans - Nk8192_FOU)) * 1 / n


        'MUSCL'

        datas = np.load('flying/results_Burger_8_MUSCL_LLF_RK3_nolim_2.npy', allow_pickle=True)
        for data in datas[1:]:
            t8_MUSCL = data[2]
            Nk8_MUSCL = data[12][0].flatten()
            n = 8
            Nk8_ans = np.empty_like(Nk8_MUSCL)
            for i in range(n):
                Nk8_ans[i] = root_scalar(f, args=(x8[i], t), method='toms748', bracket=[-1, 1]).root
            e8_L1_MUSCL = np.sum(abs(Nk8_ans - Nk8_MUSCL)) * 1 / n


        datas = np.load('flying/results_Burger_16_MUSCL_LLF_RK3_nolim_4.npy', allow_pickle=True)
        for data in datas[1:]:
            t16_MUSCL = data[2]
            Nk16_MUSCL = data[12][0].flatten()
            n = 16
            x16 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk16_ans = np.empty_like(Nk16_MUSCL)
            for i in range(n):
                Nk16_ans[i] = root_scalar(f, args=(x16[i], t), method='toms748', bracket=[-1, 1]).root
            e16_L1_MUSCL = np.sum(abs(Nk16_ans - Nk16_MUSCL)) * 1 / n
            R16_L1_MUSCL = math.log(e8_L1_MUSCL/e16_L1_MUSCL,2)


        datas = np.load('flying/results_Burger_32_MUSCL_LLF_RK3_nolim_6.npy', allow_pickle=True)
        for data in datas[1:]:
            t32_MUSCL = data[2]
            Nk32_MUSCL = data[12][0].flatten()
            n = 32
            x32 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk32_ans = np.empty_like(Nk32_MUSCL)
            for i in range(n):
                Nk32_ans[i] = root_scalar(f, args=(x32[i], t), method='toms748', bracket=[-1, 1]).root
            e32_L1_MUSCL = np.sum(abs(Nk32_ans - Nk32_MUSCL)) * 1 / n
            R32_L1_MUSCL = math.log(e16_L1_MUSCL/e32_L1_MUSCL,2)

        datas = np.load('flying/results_Burger_64_MUSCL_LLF_RK3_nolim_11.npy', allow_pickle=True)
        for data in datas[1:]:
            t64_MUSCL = data[2]
            Nk64_MUSCL = data[12][0].flatten()
            n = 64
            x64 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk64_ans = np.empty_like(Nk64_MUSCL)
            for i in range(n):
                Nk64_ans[i] = root_scalar(f, args=(x64[i], t), method='toms748', bracket=[-1, 1]).root
            e64_L1_MUSCL = np.sum(abs(Nk64_ans - Nk64_MUSCL)) * 1 / n
            R64_L1_MUSCL = math.log(e32_L1_MUSCL/e64_L1_MUSCL,2)

        datas = np.load('flying/results_Burger_128_MUSCL_LLF_RK3_nolim_21.npy', allow_pickle=True)
        for data in datas[1:]:
            t128_MUSCL = data[2]
            Nk128_MUSCL = data[12][0].flatten()
            n = 128
            x128 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk128_ans = np.empty_like(Nk128_MUSCL)
            for i in range(n):
                Nk128_ans[i] = root_scalar(f, args=(x128[i], t), method='toms748', bracket=[-1, 1]).root
            e128_L1_MUSCL = np.sum(abs(Nk128_ans - Nk128_MUSCL)) * 1 / n
            R128_L1_MUSCL = math.log(e64_L1_MUSCL/e128_L1_MUSCL,2)

        datas = np.load('flying/results_Burger_256_MUSCL_LLF_RK3_nolim_42.npy', allow_pickle=True)
        for data in datas[1:]:
            t256_MUSCL = data[2]
            Nk256_MUSCL = data[12][0].flatten()
            n = 256
            x256 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk256_ans = np.empty_like(Nk256_MUSCL)
            for i in range(n):
                Nk256_ans[i] = root_scalar(f, args=(x256[i], t), method='toms748', bracket=[-1, 1]).root
            e256_L1_MUSCL = np.sum(abs(Nk256_ans - Nk256_MUSCL)) * 1 / n
            R256_L1_MUSCL = math.log(e128_L1_MUSCL/e256_L1_MUSCL,2)


        datas = np.load('flying/results_Burger_512_MUSCL_LLF_RK3_nolim_80.npy', allow_pickle=True)
        for data in datas[1:]:
            t512_MUSCL = data[2]
            Nk512_MUSCL = data[12][0].flatten()
            n = 512
            x512 = np.linspace(0+(1)*1/(2*n),1-(1)*1/(2*n),n)
            Nk512_ans = np.empty_like(Nk512_MUSCL)
            for i in range(n):
                Nk512_ans[i] = root_scalar(f, args=(x512[i], t), method='toms748', bracket=[-1, 1]).root
            e512_L1_MUSCL = np.sum(abs(Nk512_ans - Nk512_MUSCL)) * 1 / n
            R512_L1_MUSCL = math.log(e256_L1_MUSCL/e512_L1_MUSCL,2)


        plt.figure(1)
        x = np.linspace(0+1/20000, 1-1/20000, 10000)
        plt.plot(x, Nk_FOU, 'k', x512_4, Nk512_ans_4, 'b')
        plt.legend(('FOU', 'Analytical'))
        plt.savefig('results/compositional/FR/Nk_Burgers_sol_03t')

        plt.figure(2)
        x = np.linspace(0+1/20000, 1-1/20000, 10000)
        plt.plot(x512_4, Nk512_ans_4, 'k')
        plt.plot(x8_2, Nk8_FR, '-g*', x8_3, Nk8_FR3, '-bo', x8_4, Nk8_FR4, '-ys')
        plt.legend(('Analytical', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.savefig('results/compositional/FR/Nk_Burgers_8.png')

        plt.figure(3)
        x = np.log10(np.array([8,16,32,64,128,256,512]))
        y = np.log10(np.array([e8_L1, e16_L1, e32_L1, e64_L1, e128_L1, e256_L1, e512_L1]))
        y_FR3 = np.log10(np.array([e8_L1_3, e16_L1_3, e32_L1_3, e64_L1_3, e128_L1_3, e256_L1_3, e512_L1_3]))
        y_FR4 = np.log10(np.array([e8_L1_4, e16_L1_4, e32_L1_4, e64_L1_4, e128_L1_4, e256_L1_4, e512_L1_4]))
        #y_FR5 = np.log10(np.array([e8_L1_5, e16_L1_5, e32_L1_5, e64_L1_5, e128_L1_5, e256_L1_5]))
        y_MUSCL = np.log10(np.array([e8_L1_MUSCL, e16_L1_MUSCL, e32_L1_MUSCL, e64_L1_MUSCL, e128_L1_MUSCL, e256_L1_MUSCL, e512_L1_MUSCL]))
        y_ref = -2*x +2
        plt.plot(x, y, '-ro', x, y_FR3, '-g^', x, y_FR4, '-ys', x, y_MUSCL, 'm*', mfc='none')
        plt.grid()
        plt.ylabel('$log_{10}({E}_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.legend(('FR - $\mathcal{P}_1$', 'FR - $\mathcal{P}_2$', 'FR - $\mathcal{P}_3$', 'MUSCL'))
        plt.savefig('results/compositional/FR/Nk_Burgers_convergence.png')

        plt.figure(4)
        x32 = np.linspace(0+1/64,1-1/64,32)
        plt.plot(x32, Nk32_FR_avg, '-r^', x32, Nk32_FR3_avg, '-go', x32, Nk32_FR4_avg, '-ys', mfc='none')
        plt.plot(x32, Nk32_MUSCL, '-m*', x512_4, Nk512_ans_4, 'k', mfc='none')
        plt.grid()
        plt.xlim(0.18, 0.41)
        plt.ylim(0.1, 0.17)
        plt.ylabel('$N_k$')
        plt.xlabel('Distance')
        plt.title('Results for t=0.3 with mesh 32x1x1')
        plt.legend(('FR-P1', 'FR-P2', 'FR-P3','MUSCL', 'Semi-Analytical'))
        plt.savefig('results/compositional/FR/Nk_Burgers_32_ing.png')

        plt.figure(6)
        plt.plot(x16, Nk16_FR_avg, '-r^', x16, Nk16_FR3_avg, '-go', x16, Nk16_FR4_avg, '-ys', mfc='none')
        plt.plot(x16, Nk16_FOU, 'b<', x16, Nk16_MUSCL, '-m*', mfc='none')
        plt.plot(x512_4, Nk512_ans_4, 'k')
        plt.grid()
        plt.ylabel('$N_k$')
        plt.xlabel('Distance')
        plt.title('Results for t=0.3 with mesh 16x1x1')
        plt.legend(('FR-P1', 'FR-P2', 'FR-P3', 'FOU',  'MUSCL', 'Semi-Analytical'))
        plt.savefig('results/compositional/FR/Nk_Burgers_16_ing.png')

        plt.figure(7)
        x16 = np.linspace(0+1/32,1-1/32,16)
        plt.plot(x8, Nk8_FR_avg, '-r^', x8, Nk8_FR3_avg, '-go', x8,\
            Nk8_FR4_avg, '-ys', x8, Nk8_FOU, '-b<', x8, Nk8_MUSCL, '-m*', mfc='none')
        plt.plot(x512_4, Nk512_ans_4, 'k')
        plt.grid()
        #plt.xlim(0.18, 0.41)
        #plt.ylim(0.1, 0.17)
        plt.ylabel('$N_k$')
        plt.xlabel('Distance')
        plt.title('Results for t=0.3 with mesh 8x1x1')
        plt.legend(('FR-P1', 'FR-P2', 'FR-P3', 'FOU','MUSCL', 'Semi-Analytical'))
        plt.savefig('results/compositional/FR/Nk_Burgers_8_ing.png')

        plt.figure(8)
        x16 = np.linspace(0+1/32,1-1/32,16)
        plt.plot(x128, Nk128_FR_avg, '-r^', x128, Nk128_FR3_avg, '-go', x128,\
            Nk128_FR4_avg, '-ys', x128, Nk128_FOU, '-b<', x128, Nk128_MUSCL, '-m*', mfc='none')
        plt.plot(x512_4, Nk512_ans_4, 'k')
        plt.grid()
        #plt.xlim(0.18, 0.41)
        #plt.ylim(0.1, 0.17)
        plt.ylabel('$N_k$')
        plt.xlabel('Distance')
        plt.title('Results for t=0.3 with mesh 128x1x1')
        plt.legend(('FR-P1', 'FR-P2', 'FR-P3', 'FOU', 'Semi-Analytical'))
        plt.savefig('results/compositional/FR/Nk_Burgers_128_ing.png')

        plt.figure(9)
        e_FOU = np.array([e8_L1_FOU, e16_L1_FOU, e32_L1_FOU, e64_L1_FOU, e128_L1_FOU, \
            e256_L1_FOU, e512_L1_FOU, e1024_L1_FOU, e2048_L1_FOU, e4096_L1_FOU, \
            ])#e8192_L1_FOU, e16384_L1_FOU])
        e_MUSCL = np.array([e8_L1_MUSCL, e16_L1_MUSCL, e32_L1_MUSCL, e64_L1_MUSCL, \
            e128_L1_MUSCL, e256_L1_MUSCL, e512_L1_MUSCL])
        e_FR2 = np.array([e8_L1, e16_L1, e32_L1, e64_L1, e128_L1, e256_L1, e512_L1])
        e_FR3 = np.array([e8_L1_3, e16_L1_3, e32_L1_3, e64_L1_3, e128_L1_3, e256_L1_3, e512_L1_3])
        e_FR4 = np.array([e8_L1_4, e16_L1_4, e32_L1_4, e64_L1_4, e128_L1_4, e256_L1_4, e512_L1_4])
        t_FOU = np.array([t8_FOU, t16_FOU, t32_FOU, t64_FOU, t128_FOU, t256_FOU, \
            t512_FOU, t1024_FOU, t2048_FOU, t4096_FOU])#t8192_FOU, t16384_FOU])
        t_FR2 = np.array([t8_FR2, t16_FR2, t32_FR2, t64_FR2, t128_FR2, t256_FR2, t512_FR2])
        t_FR3 = np.array([t8_FR3, t16_FR3, t32_FR3, t64_FR3, t128_FR3, t256_FR3, t512_FR3])
        t_FR4 = np.array([t8_FR4, t16_FR4, t32_FR4, t64_FR4, t128_FR4, t256_FR4, t512_FR4])
        t_MUSCL = np.array([t8_MUSCL, t16_MUSCL, t32_MUSCL, t64_MUSCL, t128_MUSCL, t256_MUSCL, t512_MUSCL])
        plt.plot(t_FOU, e_FOU, '-bo', t_MUSCL, e_MUSCL, '-mv', t_FR2, e_FR2, '-gs', t_FR3,\
            e_FR3, '-y<', t_FR4, e_FR4, '-r*', mfc='none')
        plt.grid()
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e-10, 1e-1)
        minor_ticks = np.array([1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1])
        major_ticks = np.array([1e-10, 1e-8, 1e-6, 1e-4, 1e-2])
        plt.yticks(minor_ticks)
        #plt.grid(which='minor', axis='y')
        plt.ylabel('$|E_{L1}|$')
        plt.xlabel('CPU time [s]')
        plt.title('Results for t=0.3')
        plt.legend(('FOU', 'MUSCL', 'FR-P1', 'FR-P2', 'FR-P3'))
        plt.savefig('results/compositional/FR/eL1_time_Burgers_FR_FOU.png')

        plt.figure(10)
        plt.plot(x1024, Nk1024_FOU, '-b<', mfc='none')
        plt.plot(x512_4, Nk512_ans_4, 'k')
        plt.grid()
        plt.ylabel('$N_k$')
        plt.xlabel('Distance')
        plt.title('Results for t=0.3 with mesh 1024x1x1')
        plt.legend(('FOU', 'Semi-Analytical'))
        plt.savefig('results/compositional/FR/Nk_Burgers_FOU_1024_ing.png')

        plt.figure(11)
        plt.plot(x512, Nk512_FOU, '-b<', mfc='none')
        plt.plot(x512_4, Nk512_ans_4, 'k')
        plt.grid()
        plt.ylabel('$N_k$')
        plt.xlabel('Distance')
        plt.title('Results for t=0.3 with mesh 512x1x1')
        plt.legend(('FOU', 'Semi-Analytical'))
        plt.savefig('results/compositional/FR/Nk_Burgers_FOU_512_ing.png')

        plt.figure(12)
        plt.plot(x512, Nk512_FR_avg, '-r^', x512, Nk512_FR3_avg, '-go', x512,\
            Nk512_FR4_avg, '-ys', x512, Nk512_FOU, '-b<', x512, Nk512_MUSCL, '-m*', mfc='none')
        plt.plot(x512_4, Nk512_ans_4, 'k')
        plt.grid()
        #plt.xlim(0.18, 0.41)
        #plt.ylim(0.1, 0.17)
        plt.ylabel('$N_k$')
        plt.xlabel('Distance')
        plt.title('Results for t=0.3 with mesh 512x1x1')
        plt.legend(('FR-P1', 'FR-P2', 'FR-P3', 'FOU', 'Semi-Analytical'))
        plt.savefig('results/compositional/FR/Nk_Burgers_512_ing.png')
        import pdb; pdb.set_trace()
