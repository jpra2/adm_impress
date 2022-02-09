import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
import quadpy

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
t = 1
x = np.linspace(-1,1,10000)
Nk_ans = np.sin(math.pi*(t-x))
f = interp1d(x,Nk_ans)

for arq in arquivos:
    if  arq.startswith(name):
        '------------------------------MUSCL-----------------------------------'
        datas = np.load('flying/results_Linear_adv_32_MUSCL_1000.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk32_MUSCL = data[12][0]
            n = 32
            x32 = np.linspace(-1+1/n,1-1/n,n)
            Nk32_ans = np.sin(math.pi*(t-x32))
            e32_L1_MUSCL = np.sum(abs(Nk32_ans - Nk32_MUSCL)) * 2 / n

        datas = np.load('flying/results_Linear_adv_64_MUSCL_1000.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk64_MUSCL = data[12][0]
            n = 64
            x64 = np.linspace(-1+1/n,1-1/n,n)
            Nk64_ans = np.sin(math.pi*(t-x64))
            e64_L1_MUSCL = np.sum(abs(Nk64_ans - Nk64_MUSCL)) * 2 / n
            R64_MUSCL = math.log(e32_L1_MUSCL/e64_L1_MUSCL,2)

        datas = np.load('flying/results_Linear_adv_128_MUSCL_1000.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk128_MUSCL = data[12][0]
            n = 128
            x128 = np.linspace(-1+1/n,1-1/n,n)
            Nk128_ans = np.sin(math.pi*(t-x128))
            e128_L1_MUSCL = np.sum(abs(Nk128_ans - Nk128_MUSCL)) * 2 / n
            R128_MUSCL = math.log(e64_L1_MUSCL/e128_L1_MUSCL,2)

        datas = np.load('flying/results_Linear_adv_256_MUSCL_10001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk256_MUSCL = data[12][0]
            n = 256
            x256 = np.linspace(-1+1/n,1-1/n,n)
            Nk256_ans = np.sin(math.pi*(t-x256))
            e256_L1_MUSCL = np.sum(abs(Nk256_ans - Nk256_MUSCL)) * 2 / n
            R256_MUSCL = math.log(e128_L1_MUSCL/e256_L1_MUSCL,2)

        datas = np.load('flying/results_Linear_adv_512_MUSCL_10001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk512_MUSCL = data[12][0]
            n = 512
            x512 = np.linspace(-1+1/n,1-1/n,n)
            Nk512_ans = np.sin(math.pi*(t-x512))
            e512_L1_MUSCL = np.sum(abs(Nk512_ans - Nk512_MUSCL)) * 2 / n
            R512_MUSCL = (np.log(e256_L1_MUSCL) - np.log(e512_L1_MUSCL))/(np.log(2/256) - np.log(2/512))



        '---------------------------Flux Reconstruction------------------------'
        datas = np.load('flying/results_Linear_adv_8_FR2_1000.npy', allow_pickle=True)

        for data in datas[1:]:
            xx = data[12][0]
            Nk8_FR = data[12][0].flatten()
            n = 8
            x8_2 = np.empty((n,2))
            GL = quadpy.c1.gauss_lobatto(2)
            points = GL.points
            for i in range(len(points)):
                x8_2[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x8_2 = x8_2.flatten()
            Nk8_ans_2 = np.sin(math.pi*(t-x8_2))
            e8_L1 = np.sum(abs(Nk8_ans_2 - Nk8_FR)) * 2 / n
            e8_L2 = np.sqrt(np.sum((Nk8_ans_2 - Nk8_FR)**2) * 2 / n)

            Nk8_FR_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            x8_cell = np.linspace(-1+1/n,1-1/n, n)
            Nk8_ans_2_cell = np.sin(math.pi*(t-x8_cell))
            e8_L1_cell = np.sum(abs(Nk8_ans_2_cell - Nk8_FR_cell)) * 2 / n

        datas = np.load('flying/results_Linear_adv_16_FR2_1000.npy', allow_pickle=True)

        for data in datas[1:]:
            xx = data[12][0]
            Nk16_FR = data[12][0].flatten()
            n = 16
            x16_2 = np.empty((n,2))
            for i in range(len(points)):
                x16_2[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x16_2 = x16_2.flatten()
            Nk16_ans_2 = np.sin(math.pi*(t-x16_2))
            e16_L1 = np.sum(abs(Nk16_ans_2 - Nk16_FR)) * 2 / n
            e16_L2 = np.sqrt(np.sum((Nk16_ans_2 - Nk16_FR)**2) * 2 / n)
            R16_L1 = math.log(e8_L1/e16_L1,2)
            R16_L2 = math.log(e8_L2/e16_L2,2)

            Nk16_FR_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            x16_cell = np.linspace(-1+1/n,1-1/n, n)
            Nk16_ans_2_cell = np.sin(math.pi*(t-x16_cell))
            e16_L1_cell = np.sum(abs(Nk16_ans_2_cell - Nk16_FR_cell)) * 2 / n

        datas = np.load('flying/results_Linear_adv_32_FR2_1000.npy', allow_pickle=True)

        for data in datas[1:]:
            xx = data[12][0]
            Nk32_FR = data[12][0].flatten()
            n = 32
            x32_2 = np.empty((n,2))
            for i in range(len(points)):
                x32_2[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x32_2 = x32_2.flatten()
            Nk32_ans_2 = np.sin(math.pi*(t-x32_2))
            e32_L1 = np.sum(abs(Nk32_ans_2 - Nk32_FR)) * 2 / n
            e32_L2 = np.sqrt(np.sum((Nk32_ans_2 - Nk32_FR)**2) * 2 / n)
            R32_L1 = math.log(e16_L1/e32_L1,2)
            R32_L2 =  math.log(e16_L2/e32_L2,2)

            Nk32_FR_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            x32_cell = np.linspace(-1+1/n,1-1/n, n)
            Nk32_ans_2_cell = np.sin(math.pi*(t-x32_cell))
            e32_L1_cell = np.sum(abs(Nk32_ans_2_cell - Nk32_FR_cell)) * 2 / n

        datas = np.load('flying/results_Linear_adv_64_FR2_1000.npy', allow_pickle=True)

        for data in datas[1:]:
            xx=data[12][0]
            Nk64_FR = data[12][0].flatten()
            n = 64
            x64_2 = np.empty((n,2))
            for i in range(len(points)):
                x64_2[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x64_2 = x64_2.flatten()
            Nk64_ans_2 = np.sin(math.pi*(t-x64_2))
            e64_L1 = np.sum(abs(Nk64_ans_2 - Nk64_FR)) * 2 / n
            e64_L2 = np.sqrt(np.sum((Nk64_ans_2 - Nk64_FR)**2) * 2 / n)
            R64_L1 = math.log(e32_L1/e64_L1,2)
            R64_L2 = math.log(e32_L2/e64_L2,2)

            Nk64_FR_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            x64_cell = np.linspace(-1+1/n,1-1/n, n)
            Nk64_ans_2_cell = np.sin(math.pi*(t-x64_cell))
            e64_L1_cell = np.sum(abs(Nk64_ans_2_cell - Nk64_FR_cell)) * 2 / n
            import pdb; pdb.set_trace()
            
        datas = np.load('flying/results_Linear_adv_128_FR2_1000.npy', allow_pickle=True)

        for data in datas[1:]:
            xx=data[12][0]
            Nk128_FR = data[12][0].flatten()
            n = 128
            x128_2 = np.empty((n,2))
            for i in range(len(points)):
                x128_2[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x128_2 = x128_2.flatten()
            Nk128_ans_2 = np.sin(math.pi*(t-x128_2))
            e128_L1 = np.sum(abs(Nk128_ans_2 - Nk128_FR)) * 2 / n
            e128_L2 = np.sqrt(np.sum((Nk128_ans_2 - Nk128_FR)**2) * 2 / n)
            R128_L1 = math.log(e64_L1/e128_L1,2)
            R128_L2 = math.log(e64_L2/e128_L2,2)

            Nk128_FR_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            x128_cell = np.linspace(-1+1/n,1-1/n, n)
            Nk128_ans_2_cell = np.sin(math.pi*(t-x128_cell))
            e128_L1_cell = np.sum(abs(Nk128_ans_2_cell - Nk128_FR_cell)) * 2 / n


        datas = np.load('flying/results_Linear_adv_256_FR2_1000.npy', allow_pickle=True)

        for data in datas[1:]:
            xx = data[12][0]
            Nk256_FR = data[12][0].flatten()
            n = 256
            x256_2 = np.empty((n,2))
            for i in range(len(points)):
                x256_2[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x256_2 = x256_2.flatten()
            Nk256_ans_2 = np.sin(math.pi*(t-x256_2))
            e256_L1 = np.sum(abs(Nk256_ans_2 - Nk256_FR)) * 2 / n
            e256_L2 = np.sqrt(np.sum((Nk256_ans_2 - Nk256_FR)**2) * 2 / n)
            R256_L1 = math.log(e128_L1/e256_L1,2)
            R256_L2 = math.log(e128_L2/e256_L2,2)

            Nk256_FR_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            x256_cell = np.linspace(-1+1/n,1-1/n, n)
            Nk256_ans_2_cell = np.sin(math.pi*(t-x256_cell))
            e256_L1_cell = np.sum(abs(Nk256_ans_2_cell - Nk256_FR_cell)) * 2 / n


        datas = np.load('flying/results_Linear_adv_512_FR2_1000.npy', allow_pickle=True)

        for data in datas[1:]:
            xx = data[12][0]
            Nk512_FR = data[12][0].flatten()
            n = 512
            x512_2 = np.empty((n,2))
            for i in range(len(points)):
                x512_2[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)

            x512_2 = x512_2.flatten()
            Nk512_ans_2 = np.sin(math.pi*(t-x512_2))
            e512_L1 = np.sum(abs(Nk512_ans_2 - Nk512_FR)) * 2 / n
            e512_L2 = np.sqrt(np.sum((Nk512_ans_2 - Nk512_FR)**2) * 2 / n)
            R512_L1 = (np.log(e256_L1) - np.log(e512_L1))/(np.log(2/256) - np.log(2/512))
            R512_L2 = math.log(e256_L2/e512_L2,2)

            Nk512_FR_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            x512_cell = np.linspace(-1+1/n,1-1/n, n)
            Nk512_ans_2_cell = np.sin(math.pi*(t-x512_cell))
            e512_L1_cell = np.sum(abs(Nk512_ans_2_cell - Nk512_FR_cell)) * 2 / n

        #    p_resp = np.linspace(0.623843,0,100)

        '---------------------------Flux Reconstruction------------------------'
        datas = np.load('flying/results_Linear_adv_8_FR3_2001.npy', allow_pickle=True)

        for data in datas[1:]:
            xx = data[12][0]

            Nk8_FR3 = data[12][0].flatten()
            n = 8
            x8_3 = np.empty((n,3))
            GL = quadpy.c1.gauss_lobatto(3)
            points = GL.points
            for i in range(3):
                x8_3[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x8_3 = x8_3.flatten()
            Nk8_ans_3 = np.sin(math.pi*(t-x8_3))
            e8_L1_3 = np.sum(abs(Nk8_ans_3 - Nk8_FR3)) * 2 / n

            Nk8_FR3_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            Nk8_ans_3_cell = np.sin(math.pi*(t-x8_cell))
            e8_L1_3_cell = np.sum(abs(Nk8_ans_3_cell - Nk8_FR3_cell)) * 2 / n


        datas = np.load('flying/results_Linear_adv_16_FR3_2001.npy', allow_pickle=True)

        for data in datas[1:]:
            xx = data[12][0]
            Nk16_FR3 = data[12][0].flatten()
            n = 16
            x16_3 = np.empty((n,3))
            for i in range(3):
                x16_3[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x16_3 = x16_3.flatten()
            Nk16_ans_3 = np.sin(math.pi*(t-x16_3))
            e16_L1_3 = np.sum(abs(Nk16_ans_3 - Nk16_FR3)) * 2 / n
            R16_3 = math.log(e8_L1_3/e16_L1_3,2)

            Nk16_FR3_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            Nk16_ans_3_cell = np.sin(math.pi*(t-x16_cell))
            e16_L1_3_cell = np.sum(abs(Nk16_ans_3_cell - Nk16_FR3_cell)) * 2 / n


        datas = np.load('flying/results_Linear_adv_32_FR3_2001.npy', allow_pickle=True)

        for data in datas[1:]:
            xx = data[12][0]
            Nk32_FR3 = data[12][0].flatten()
            n = 32
            x32_3 = np.empty((n,3))
            for i in range(3):
                x32_3[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x32_3 = x32_3.flatten()
            Nk32_ans_3 = np.sin(math.pi*(t-x32_3))
            e32_L1_3 = np.sum(abs(Nk32_ans_3 - Nk32_FR3)) * 2 / n
            R32_3 = math.log(e16_L1_3/e32_L1_3,2)

            Nk32_FR3_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            Nk32_ans_3_cell = np.sin(math.pi*(t-x32_cell))
            e32_L1_3_cell = np.sum(abs(Nk32_ans_3_cell - Nk32_FR3_cell)) * 2 / n


        datas = np.load('flying/results_Linear_adv_64_FR3_2001.npy', allow_pickle=True)

        for data in datas[1:]:
            xx = data[12][0]
            Nk64_FR3 = data[12][0].flatten()
            n = 64
            x64_3 = np.empty((n,3))
            for i in range(3):
                x64_3[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x64_3 = x64_3.flatten()
            Nk64_ans_3 = np.sin(math.pi*(t-x64_3))
            e64_L1_3 = np.sum(abs(Nk64_ans_3 - Nk64_FR3)) * 2 / n
            R64_3 = math.log(e32_L1_3/e64_L1_3,2)

            Nk64_FR3_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            Nk64_ans_3_cell = np.sin(math.pi*(t-x64_cell))
            e64_L1_3_cell = np.sum(abs(Nk64_ans_3_cell - Nk64_FR3_cell)) * 2 / n


        datas = np.load('flying/results_Linear_adv_128_FR3_2001.npy', allow_pickle=True)

        for data in datas[1:]:
            xx = data[12][0]
            Nk128_FR3 = data[12][0].flatten()
            n = 128
            x128_3 = np.empty((n,3))
            for i in range(3):
                x128_3[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x128_3 = x128_3.flatten()
            Nk128_ans_3 = np.sin(math.pi*(t-x128_3))
            e128_L1_3 = np.sum(abs(Nk128_ans_3 - Nk128_FR3)) * 2 / n
            R128_3 = math.log(e64_L1_3/e128_L1_3,2)

            Nk128_FR3_cell = 1/2*np.sum(GL.weights * xx, axis=-1)
            Nk128_ans_3_cell = np.sin(math.pi*(t-x128_cell))
            e128_L1_3_cell = np.sum(abs(Nk128_ans_3_cell - Nk128_FR3_cell)) * 2 / n


        datas = np.load('flying/results_Linear_adv_256_FR3_2001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk256_FR3 = data[12][0].flatten()
            n = 256
            x256_3 = np.empty((n,3))
            for i in range(3):
                x256_3[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x256_3 = x256_3.flatten()
            Nk256_ans_3 = np.sin(math.pi*(t-x256_3))
            e256_L1_3 = np.sum(abs(Nk256_ans_3 - Nk256_FR3)) * 2 / n
            R256_3 = math.log(e128_L1_3/e256_L1_3,2)

        datas = np.load('flying/results_Linear_adv_512_FR3_2001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk512_FR3 = data[12][0].flatten()
            n = 512
            x512_3 = np.empty((n,3))
            for i in range(3):
                x512_3[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x512_3 = x512_3.flatten()
            Nk512_ans_3 = np.sin(math.pi*(t-x512_3))
            e512_L1_3 = np.sum(abs(Nk512_ans_3 - Nk512_FR3)) * 2 / n
            R512_3 = (np.log(e256_L1_3) - np.log(e512_L1_3))/(np.log(2/256) - np.log(2/512))

        '''----------------------------FR4-------------------------'''

        datas = np.load('flying/results_Linear_adv_8_FR4_3334.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk8_FR4 = data[12][0].flatten()
            n = 8
            x8_4 = np.empty((n,4))
            GL = quadpy.c1.gauss_lobatto(4)
            points = GL.points
            for i in range(4):
                x8_4[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x8_4 = x8_4.flatten()
            Nk8_ans_4 = np.sin(math.pi*(t-x8_4))
            e8_L1_4 = np.sum(abs(Nk8_ans_4 - Nk8_FR4)) * 2 / n

        datas = np.load('flying/results_Linear_adv_16_FR4_3334.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk16_FR4 = data[12][0].flatten()
            n = 16
            x16_4 = np.empty((n,4))
            for i in range(4):
                x16_4[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x16_4 = x16_4.flatten()
            Nk16_ans_4 = np.sin(math.pi*(t-x16_4))
            e16_L1_4 = np.sum(abs(Nk16_ans_4 - Nk16_FR4)) * 2 / n
            R16_4 = math.log(e8_L1_4/e16_L1_4, 2)

        datas = np.load('flying/results_Linear_adv_32_FR4_3334.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk32_FR4 = data[12][0].flatten()
            n = 32
            x32_4 = np.empty((n,4))
            for i in range(4):
                x32_4[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x32_4 = x32_4.flatten()
            Nk32_ans_4 = np.sin(math.pi*(t-x32_4))
            e32_L1_4 = np.sum(abs(Nk32_ans_4 - Nk32_FR4)) * 2 / n
            R32_4 = math.log(e16_L1_4/e32_L1_4,2)

        datas = np.load('flying/results_Linear_adv_64_FR4_3334.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk64_FR4 = data[12][0].flatten()
            n = 64
            x64_4 = np.empty((n,4))
            for i in range(4):
                x64_4[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x64_4 = x64_4.flatten()
            Nk64_ans_4 = np.sin(math.pi*(t-x64_4))
            e64_L1_4 = np.sum(abs(Nk64_ans_4 - Nk64_FR4)) * 2 / n
            R64_4 = math.log(e32_L1_4/e64_L1_4,2)

        datas = np.load('flying/results_Linear_adv_128_FR4_5001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk128_FR4 = data[12][0].flatten()
            n = 128
            x128_4 = np.empty((n,4))
            for i in range(4):
                x128_4[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x128_4 = x128_4.flatten()
            Nk128_ans_4 = np.sin(math.pi*(t-x128_4))
            e128_L1_4 = np.sum(abs(Nk128_ans_4 - Nk128_FR4)) * 2 / n
            R128_4 = math.log(e64_L1_4/e128_L1_4,2)

        datas = np.load('flying/results_Linear_adv_256_FR4_20001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk256_FR4 = data[12][0].flatten()
            n = 256
            x256_4 = np.empty((n,4))
            for i in range(4):
                x256_4[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x256_4 = x256_4.flatten()
            Nk256_ans_4 = np.sin(math.pi*(t-x256_4))
            e256_L1_4 = np.sum(abs(Nk256_ans_4 - Nk256_FR4)) * 2 / n
            R256_4 = math.log(e128_L1_4/e256_L1_4,2)

        datas = np.load('flying/results_Linear_adv_512_FR4_20001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk512_FR4 = data[12][0].flatten()
            n = 512
            x512_4 = np.empty((n,4))
            for i in range(4):
                x512_4[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x512_4 = x512_4.flatten()
            Nk512_ans_4 = np.sin(math.pi*(t-x512_4))
            e512_L1_4 = np.sum(abs(Nk512_ans_4 - Nk512_FR4)) * 2 / n
            R512_4 = math.log(e256_L1_4/e512_L1_4,2)

        '''----------------------------FR5-------------------------'''

        datas = np.load('flying/results_Linear_adv_8_FR5_10001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk8_FR5 = data[12][0].flatten()
            n = 8
            x8_5 = np.empty((n,5))
            GL = quadpy.c1.gauss_lobatto(5)
            points = GL.points
            for i in range(5):
                x8_5[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x8_5 = x8_5.flatten()
            Nk8_ans_5 = np.sin(math.pi*(t-x8_5))
            e8_L1_5 = np.sum(abs(Nk8_ans_5 - Nk8_FR5)) * 2 / n

        datas = np.load('flying/results_Linear_adv_16_FR5_10001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk16_FR5 = data[12][0].flatten()
            n = 16
            x16_5 = np.empty((n,5))
            for i in range(5):
                x16_5[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x16_5 = x16_5.flatten()
            Nk16_ans_5 = np.sin(math.pi*(t-x16_5))
            e16_L1_5 = np.sum(abs(Nk16_ans_5 - Nk16_FR5)) * 2 / n
            R16_5 = math.log(e8_L1_5/e16_L1_5, 2)


        datas = np.load('flying/results_Linear_adv_32_FR5_12501.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk32_FR5 = data[12][0].flatten()
            n = 32
            x32_5 = np.empty((n,5))
            for i in range(5):
                x32_5[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x32_5 = x32_5.flatten()
            Nk32_ans_5 = np.sin(math.pi*(t-x32_5))
            e32_L1_5 = np.sum(abs(Nk32_ans_5 - Nk32_FR5)) * 2 / n
            R32_5 = math.log(e16_L1_5/e32_L1_5,2)


        datas = np.load('flying/results_Linear_adv_64_FR5_12501.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk64_FR5 = data[12][0].flatten()
            n = 64
            x64_5 = np.empty((n,5))
            for i in range(5):
                x64_5[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x64_5 = x64_5.flatten()
            Nk64_ans_5 = np.sin(math.pi*(t-x64_5))
            e64_L1_5 = np.sum(abs(Nk64_ans_5 - Nk64_FR5)) * 2 / n
            R64_5 = math.log(e32_L1_5/e64_L1_5,2)

        datas = np.load('flying/results_Linear_adv_128_FR5_20001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk128_FR5 = data[12][0].flatten()
            n = 128
            x128_5 = np.empty((n,5))
            for i in range(5):
                x128_5[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x128_5 = x128_5.flatten()
            Nk128_ans_5 = np.sin(math.pi*(t-x128_5))
            e128_L1_5 = np.sum(abs(Nk128_ans_5 - Nk128_FR5)) * 2 / n
            R128_5 = math.log(e64_L1_5/e128_L1_5,2)

        datas = np.load('flying/results_Linear_adv_256_FR5_20001.npy', allow_pickle=True)

        for data in datas[1:]:
            Nk256_FR5 = data[12][0].flatten()
            n = 256
            x256_5 = np.empty((n,5))
            for i in range(5):
                x256_5[:,i] = np.linspace(-1+(points[i] + 1)*1/n,1-(1-points[i])*1/n,n)
            x256_5 = x256_5.flatten()
            Nk256_ans_5 = np.sin(math.pi*(t-x256_5))
            e256_L1_5 = np.sum(abs(Nk256_ans_5 - Nk256_FR5)) * 2 / n
            R256_5 = math.log(e128_L1_5/e256_L1_5,2)



        plt.figure(30)
        plt.plot(x8_2, Nk8_FR, 'r', x, Nk_ans, 'k', x8_3, Nk8_FR3, 'y')
        plt.grid()
        plt.legend(('FR-2nd order', 'Analytical Solution', 'MUSCL', 'FR-3rd order'))
        plt.ylabel('Nk (moles)')
        plt.xlabel('Dimensionless distance')
        plt.title('Water horizontal displacement')
        plt.savefig('results/compositional/FR_paper/sinoidal_pulse/Nk_seno_8_func.eps', format='eps')

        plt.figure(31)

        plt.plot(x16_2, Nk16_FR, 'r', x, Nk_ans, 'k', x16_3, Nk16_FR3, 'y')
        plt.grid()
        plt.legend(('FR-2nd order', 'Analytical Solution', 'MUSCL', 'FR-3rd order'))
        plt.ylabel('Nk (moles)')
        plt.xlabel('Dimensionless distance')
        plt.title('Water horizontal displacement')
        plt.savefig('results/compositional/FR_paper/sinoidal_pulse/Nk_seno_16_func.eps', format='eps')


        plt.figure(1)

        plt.plot(x32_2, Nk32_FR, 'r', x, Nk_ans, 'k', x32, Nk32_MUSCL, 'g', x32_3, Nk32_FR3, 'y')
        plt.grid()
        plt.legend(('FR-2nd order', 'Analytical Solution', 'MUSCL', 'FR-3rd order'))
        plt.ylabel('Nk (moles)')
        plt.xlabel('Dimensionless distance')
        plt.title('Water horizontal displacement')
        plt.savefig('results/compositional/FR_paper/sinoidal_pulse/Nk_seno_32_func.eps', format='eps')

        plt.figure(2)

        plt.plot(x64_2, Nk64_FR, 'r', x, Nk_ans, 'k', x64, Nk64_MUSCL, 'g', x64_3, Nk64_FR3, 'y')
        plt.grid()
        plt.legend(('FR-2nd order', 'Analytical Solution', 'MUSCL', 'FR-3rd order'))
        plt.ylabel('Nk (moles)')
        plt.xlabel('Dimensionless distance')
        plt.title('Water horizontal displacement')
        plt.savefig('results/compositional/FR/Nk_seno_64_func.png')

        plt.figure(3)

        plt.plot(x128_2, Nk128_FR, 'r', x, Nk_ans, 'k', x128, Nk128_MUSCL, 'g', x128_3, Nk128_FR3, 'y')
        plt.grid()
        plt.legend(('FR-2nd order', 'Analytical Solution', 'MUSCL', 'FR-3rd order'))
        plt.ylabel('Nk (moles)')
        plt.xlabel('Dimensionless distance')
        plt.title('Water horizontal displacement')
        plt.savefig('results/compositional/FR_paper/sinoidal_pulse/Nk_seno_128_func.eps', format='eps')

        plt.figure(4)

        plt.plot(x256_2, Nk256_FR, 'r', x, Nk_ans, 'k', x256, Nk256_MUSCL, 'g', x256_3, Nk256_FR3, 'y')
        plt.grid()
        plt.legend(('FR-2nd order', 'Analytical Solution', 'MUSCL', 'FR-3rd order'))
        plt.ylabel('Nk (moles)')
        plt.xlabel('Dimensionless distance')
        plt.title('Water horizontal displacement')
        plt.savefig('results/compositional/FR_paper/sinoidal_pulse/Nk_seno_256_func.eps', format='eps')


        plt.figure(5)

        plt.plot(x, Nk_ans, 'k', x512, Nk512_MUSCL, 'g', x512_2, Nk512_FR, 'r', x512_3, Nk512_FR3, 'y')#, x512, Nk512_FR4, 'm')
        plt.grid()
        plt.legend(('Analytical Solution', 'MUSCL', 'FR-2nd order', 'FR-3rd order'))
        plt.ylabel('Nk (moles)')
        plt.xlabel('Dimensionless distance')
        plt.title('Water horizontal displacement')
        plt.savefig('results/compositional/FR_paper/sinoidal_pulse/Nk_seno_512_func.eps', format='eps')

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
        plt.legend(('CPR-$2^a$ ordem', 'CPR-3$^a$ ordem', 'CPR-4$^a$ ordem'))#, 'FR-5th order'))
        plt.savefig('results/compositional/FR/Nk_seno_convergencia.png')
        import pdb; pdb.set_trace()
