import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline
import quadpy
from case_plots.BL_exact_Gustavo import BL_exact

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)


fi, A, Swr, Sor, ftime, L = 0.2, 0.02592, 0, 0, 1500, 300
for  arq in arquivos:
    if  arq.startswith(name):
        '''-------------------------FR2 RESULTS------------------------'''
        datas = np.load('flying/results_BL_8_FR2_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw8_avg = data[5]
            n=8
            GL = quadpy.c1.gauss_lobatto(2)
            points = GL.points
            Sw8 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            GL = quadpy.c1.gauss_lobatto(2)
            points = GL.points
            x8 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x8[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x8 = x8.flatten()
            x8_avg = np.linspace(0+L/(2*n),L-L/(2*n),n)
            #x8 = np.linspace(0,1,n)
            S0 =   np.ones_like(x8)*0
            f8 = BL_exact(S0, x8, fi, A, Swr, Sor, ftime)
            f8[0] = 1
            e8_L1_FR2 = (sum(abs(f8-Sw8))*(300 / n))
            #e8_L1_FR2_avg = (sum(abs(f(x8_avg)-Sw8_avg))*(300 / n))

        datas = np.load('flying/results_BL_16_FR2_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw16_avg = data[5]
            n = 16
            Sw16 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x16 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x16[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x16 = x16.flatten()
            x16_avg = np.linspace(0+L/(2*n),L-L/(2*n),n)
            S0 =   np.ones_like(x16)*0
            f16 = BL_exact(S0, x16, fi, A, Swr, Sor, ftime)
            f16[0] = 1

            e16_L1_FR2 = (sum(abs(f16-Sw16))*(300 / n))
            #e16_L1_FR2_avg = (sum(abs(f(x16_avg)-Sw16_avg))*(300 / n))
            R16_L1_FR2 = math.log(e8_L1_FR2/e16_L1_FR2,2)
            #R16_L1_FR2_avg = math.log(e8_L1_FR2_avg/e16_L1_FR2_avg,2)

        datas = np.load('flying/results_BL_32_FR2_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw32_avg = data[5]
            n = 32
            Sw32 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x32 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x32[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x32 = x32.flatten()
            x32_avg = np.linspace(0+L/(2*n),L-L/(2*n),n)

            S0 =   np.ones_like(x32)*0
            f32 = BL_exact(S0, x32, fi, A, Swr, Sor, ftime)
            f32[0] = 1
            #x32= np.linspace(0,1,n)
            e32_L1_FR2 = (sum(abs(f32-Sw32))*(300 / n))
            R32_L1_FR2 = math.log(e16_L1_FR2/e32_L1_FR2,2)
            #e32_L1_FR2_avg = (sum(abs(f(x32_avg)-Sw32_avg))*(300 / n))
            #R32_L1_FR2_avg = math.log(e16_L1_FR2_avg/e32_L1_FR2_avg,2)

        datas = np.load('flying/results_BL_64_FR2_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw64_avg = data[5]
            n = 64
            Sw64 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x64 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x64[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x64 = x64.flatten()
            x64_avg = np.linspace(0+L/(2*n),L-L/(2*n),n)

            S0 =   np.ones_like(x64)*0
            f64 = BL_exact(S0, x64, fi, A, Swr, Sor, ftime)
            f64[0] = 1
            e64_L1_FR2 = (sum(abs(f64-Sw64))*(300 / n))
            R64_L1_FR2 = math.log(e32_L1_FR2/e64_L1_FR2,2)
            #e64_L1_FR2_avg = (sum(abs(f(x64_avg)-Sw64_avg))*(300 / n))
            #R64_L1_FR2_avg = math.log(e32_L1_FR2_avg/e64_L1_FR2_avg,2)

        datas = np.load('flying/results_BL_128_FR2_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw128_avg = data[5]
            n = 128
            Sw128 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x128 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x128[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x128 = x128.flatten()
            #x128 = np.linspace(0,1,n)
            x128_avg = np.linspace(0+L/(2*n),L-L/(2*n),n)

            S0 =   np.ones_like(x128)*0
            f128 = BL_exact(S0, x128, fi, A, Swr, Sor, ftime)
            f128[0] = 1
            e128_L1_FR2 = (sum(abs(f128-Sw128))*(300 / n))
            R128_L1_FR2 = math.log(e64_L1_FR2/e128_L1_FR2,2)
            #e128_L1_FR2_avg = (sum(abs(f(x128_avg)-Sw128_avg))*(300 / n))
            #R128_L1_FR2_avg = math.log(e64_L1_FR2_avg/e128_L1_FR2_avg,2)

        datas = np.load('flying/results_BL_256_FR2_6000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw256_avg = data[5]
            n = 256
            Sw256 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x256 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x256[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x256 = x256.flatten()
            #x256 = np.linspace(0,1,n)
            x256_avg = np.linspace(0+L/(2*n),L-L/(2*n),n)

            S0 =   np.ones_like(x256)*0
            f256 = BL_exact(S0, x256, fi, A, Swr, Sor, ftime)
            f256[0] = 1

            e256_L1_FR2 = (sum(abs(f256-Sw256))*(300 / n))
            R256_L1_FR2 = math.log(e128_L1_FR2/e256_L1_FR2,2)
            #e256_L1_FR2_avg = (sum(abs(f(x256_avg)-Sw256_avg))*(300 / n))
            #R256_L1_FR2_avg = math.log(e128_L1_FR2_avg/e256_L1_FR2_avg,2)


        datas = np.load('flying/results_BL_512_FR2_6000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw512_avg = data[5]
            n = 512
            Sw512 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x512 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x512[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x512 = x512.flatten()
            #x512 = np.linspace(0,1,n)

            S0 =   np.ones_like(x512)*0
            f512 = BL_exact(S0, x512, fi, A, Swr, Sor, ftime)
            f512[0] = 1

            x512_avg = np.linspace(0+L/(2*n),L-L/(2*n),n)
            e512_L1_FR2 = (sum(abs(f512-Sw512))*(300 / n))
            R512_L1_FR2 = math.log(e256_L1_FR2/e512_L1_FR2,2)
            #e512_L1_FR2_avg = (sum(abs(f(x512_avg)-Sw512_avg))*(300 / n))
            #R512_L1_FR2_avg = math.log(e256_L1_FR2_avg/e512_L1_FR2_avg,2)

        '''-------------------------FR3 RESULTS------------------------'''

        datas = np.load('flying/results_BL_8_FR3_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw8_avg_FR3 = data[5]
            n=8
            Sw8_FR3 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            GL = quadpy.c1.gauss_lobatto(3)
            points = GL.points
            x8_FR3 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x8_FR3[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x8_FR3 = x8_FR3.flatten()

            S0 =   np.ones_like(x8_FR3)*0
            f8_FR3 = BL_exact(S0, x8_FR3, fi, A, Swr, Sor, ftime)
            f8_FR3[0] = 1
            #x8 = np.linspace(0,1,n)
            e8_L1_FR3 = (sum(abs(f8_FR3-Sw8_FR3))*(300 / n))
            #e8_L1_FR3_avg = (sum(abs(f(x8_avg)-Sw8_avg_FR3))*(300 / n))

        datas = np.load('flying/results_BL_16_FR3_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw16_avg_FR3 = data[5]
            n = 16
            Sw16_FR3 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x16_FR3 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x16_FR3[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x16_FR3 = x16_FR3.flatten()

            S0 =   np.ones_like(x16_FR3)*0
            f16_FR3 = BL_exact(S0, x16_FR3, fi, A, Swr, Sor, ftime)
            f16_FR3[0] = 1

            e16_L1_FR3 = (sum(abs(f16_FR3-Sw16_FR3))*(300 / n))
            #e16_L1_FR3_avg = (sum(abs(f(x16_avg)-Sw16_avg_FR3))*(300 / n))
            R16_L1_FR3 = math.log(e8_L1_FR3/e16_L1_FR3,2)
            #R16_L1_FR3_avg = math.log(e8_L1_FR3_avg/e16_L1_FR3_avg,2)

        datas = np.load('flying/results_BL_32_FR3_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw32_avg_FR3 = data[5]
            n = 32
            Sw32_FR3 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x32_FR3 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x32_FR3[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x32_FR3 = x32_FR3.flatten()

            S0 =   np.ones_like(x32_FR3)*0
            f32_FR3 = BL_exact(S0, x32_FR3, fi, A, Swr, Sor, ftime)
            f32_FR3[0] = 1

            e32_L1_FR3 = (sum(abs(f32_FR3-Sw32_FR3))*(300 / n))
            R32_L1_FR3 = math.log(e16_L1_FR3/e32_L1_FR3,2)
            #e32_L1_FR3_avg = (sum(abs(f(x32_avg)-Sw32_avg_FR3))*(300 / n))
            #R32_L1_FR3_avg = math.log(e16_L1_FR3_avg/e32_L1_FR3_avg,2)

        datas = np.load('flying/results_BL_64_FR3_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw64_avg_FR3 = data[5]
            n = 64
            Sw64_FR3 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x64_FR3 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x64_FR3[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x64_FR3 = x64_FR3.flatten()

            S0 =   np.ones_like(x64_FR3)*0
            f64_FR3 = BL_exact(S0, x64_FR3, fi, A, Swr, Sor, ftime)
            f64_FR3[0] = 1

            e64_L1_FR3 = (sum(abs(f64_FR3-Sw64_FR3))*(300 / n))
            R64_L1_FR3 = math.log(e32_L1_FR3/e64_L1_FR3,2)
            #e64_L1_FR3_avg = (sum(abs(f(x64_avg)-Sw64_avg_FR3))*(300 / n))
            #R64_L1_FR3_avg = math.log(e32_L1_FR3_avg/e64_L1_FR3_avg,2)

        datas = np.load('flying/results_BL_128_FR3_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw128_avg_FR3 = data[5]
            n = 128
            Sw128_FR3 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x128_FR3 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x128_FR3[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x128_FR3 = x128_FR3.flatten()

            S0 =   np.ones_like(x128_FR3)*0
            f128_FR3 = BL_exact(S0, x128_FR3, fi, A, Swr, Sor, ftime)
            f128_FR3[0] = 1

            e128_L1_FR3 = (sum(abs(f128_FR3-Sw128_FR3))*(300 / n))
            R128_L1_FR3 = math.log(e64_L1_FR3/e128_L1_FR3,2)
            #e128_L1_FR3_avg = (sum(abs(f(x128_avg)-Sw128_avg_FR3))*(300 / n))
            #R128_L1_FR3_avg = math.log(e64_L1_FR3_avg/e128_L1_FR3_avg,2)

        datas = np.load('flying/results_BL_256_FR3_6000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw256_avg_FR3 = data[5]
            n = 256
            Sw256_FR3 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x256_FR3 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x256_FR3[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x256_FR3 = x256_FR3.flatten()

            S0 =   np.ones_like(x256_FR3)*0
            f256_FR3 = BL_exact(S0, x256_FR3, fi, A, Swr, Sor, ftime)
            f256_FR3[0] = 1

            e256_L1_FR3 = (sum(abs(f256_FR3-Sw256_FR3))*(300 / n))
            R256_L1_FR3 = math.log(e128_L1_FR3/e256_L1_FR3,2)
            #e256_L1_FR3_avg = (sum(abs(f(x256_avg)-Sw256_avg_FR3))*(300 / n))
            #R256_L1_FR3_avg = math.log(e128_L1_FR3_avg/e256_L1_FR3_avg,2)

        datas = np.load('flying/results_BL_512_FR3_12000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw512_avg_FR3 = data[5]
            n = 512
            Sw512_FR3 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x512_FR3 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x512_FR3[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x512_FR3 = x512_FR3.flatten()

            S0 =   np.ones_like(x512_FR3)*0
            f512_FR3 = BL_exact(S0, x512_FR3, fi, A, Swr, Sor, ftime)
            f512_FR3[0] = 1
            e512_L1_FR3 = (sum(abs(f512_FR3-Sw512_FR3))*(300 / n))
            #e512_L1_FR3_avg = (sum(abs(f(x512_avg)-Sw512_avg_FR3))*(300 / n))
            R512_L1_FR3 = math.log(e256_L1_FR3/e512_L1_FR3,2)
            #R512_L1_FR3_avg = math.log(e256_L1_FR3_avg/e512_L1_FR3_avg,2)


        '''-------------------------FR3 RESULTS------------------------'''

        datas = np.load('flying/results_BL_8_FR4_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw8_avg_FR4 = data[5]
            n=8
            Sw8_FR4 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            GL = quadpy.c1.gauss_lobatto(4)
            points = GL.points
            x8_FR4 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x8_FR4[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x8_FR4 = x8_FR4.flatten()

            S0 =   np.ones_like(x8_FR4)*0
            f8_FR4 = BL_exact(S0, x8_FR4, fi, A, Swr, Sor, ftime)
            f8_FR4[0] = 1
            #x8 = np.linspace(0,1,n)
            e8_L1_FR4 = (sum(abs(f8_FR4-Sw8_FR4))*(300 / n))
            #e8_L1_FR4_avg = (sum(abs(f(x8_avg)-Sw8_avg_FR4))*(300 / n))

        datas = np.load('flying/results_BL_16_FR4_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw16_avg_FR4 = data[5]
            n=16
            Sw16_FR4 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x16_FR4 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x16_FR4[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x16_FR4 = x16_FR4.flatten()

            S0 =   np.ones_like(x16_FR4)*0
            f16_FR4 = BL_exact(S0, x16_FR4, fi, A, Swr, Sor, ftime)
            f16_FR4[0] = 1
            e16_L1_FR4 = (sum(abs(f16_FR4-Sw16_FR4))*(300 / n))
            #e16_L1_FR4_avg = (sum(abs(f(x16_avg)-Sw16_avg_FR4))*(300 / n))
            R16_L1_FR4 = math.log(e8_L1_FR4/e16_L1_FR4,2)
            #R16_L1_FR4_avg = math.log(e8_L1_FR4_avg/e16_L1_FR4_avg,2)

        datas = np.load('flying/results_BL_32_FR4_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw32_avg_FR4 = data[5]
            n=32
            Sw32_FR4 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x32_FR4 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x32_FR4[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x32_FR4 = x32_FR4.flatten()

            S0 =   np.ones_like(x32_FR4)*0
            f32_FR4 = BL_exact(S0, x32_FR4, fi, A, Swr, Sor, ftime)
            f32_FR4[0] = 1
            e32_L1_FR4 = (sum(abs(f32_FR4-Sw32_FR4))*(300 / n))
            #e32_L1_FR4_avg = (sum(abs(f(x32_avg)-Sw32_avg_FR4))*(300 / n))
            R32_L1_FR4 = math.log(e16_L1_FR4/e32_L1_FR4,2)
            #R32_L1_FR4_avg = math.log(e16_L1_FR4_avg/e32_L1_FR4_avg,2)

        datas = np.load('flying/results_BL_64_FR4_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw64_avg_FR4 = data[5]
            n=64
            Sw64_FR4 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x64_FR4 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x64_FR4[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x64_FR4 = x64_FR4.flatten()

            S0 =   np.ones_like(x64_FR4)*0
            f64_FR4 = BL_exact(S0, x64_FR4, fi, A, Swr, Sor, ftime)
            f64_FR4[0] = 1
            e64_L1_FR4 = (sum(abs(f64_FR4-Sw64_FR4))*(300 / n))
            #e64_L1_FR4_avg = (sum(abs(f(x64_avg)-Sw64_avg_FR4))*(300 / n))
            R64_L1_FR4 = math.log(e32_L1_FR4/e64_L1_FR4,2)
            #R64_L1_FR4_avg = math.log(e32_L1_FR4_avg/e64_L1_FR4_avg,2)

        datas = np.load('flying/results_BL_128_FR4_3000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw128_avg_FR4 = data[5]
            n=128
            Sw128_FR4 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x128_FR4 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x128_FR4[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x128_FR4 = x128_FR4.flatten()

            S0 =   np.ones_like(x128_FR4)*0
            f128_FR4 = BL_exact(S0, x128_FR4, fi, A, Swr, Sor, ftime)
            f128_FR4[0] = 1
            e128_L1_FR4 = (sum(abs(f128_FR4-Sw128_FR4))*(300 / n))
            #e128_L1_FR4_avg = (sum(abs(f(x128_avg)-Sw128_avg_FR4))*(300 / n))
            R128_L1_FR4 = math.log(e64_L1_FR4/e128_L1_FR4,2)
            #R128_L1_FR4_avg = math.log(e64_L1_FR4_avg/e128_L1_FR4_avg,2)

        datas = np.load('flying/results_BL_256_FR4_6000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw256_avg_FR4 = data[5]
            n=256
            Sw256_FR4 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x256_FR4 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x256_FR4[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x256_FR4 = x256_FR4.flatten()

            S0 =   np.ones_like(x256_FR4)*0
            f256_FR4 = BL_exact(S0, x256_FR4, fi, A, Swr, Sor, ftime)
            f256_FR4[0] = 1
            e256_L1_FR4 = (sum(abs(f256_FR4-Sw256_FR4))*(300 / n))
            #e256_L1_FR4_avg = (sum(abs(f(x256_avg)-Sw256_avg_FR4))*(300 / n))
            R256_L1_FR4 = math.log(e128_L1_FR4/e256_L1_FR4,2)
            #R256_L1_FR4_avg = math.log(e128_L1_FR4_avg/e256_L1_FR4_avg,2)

        datas = np.load('flying/results_BL_512_FR4_12000.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw512_avg_FR4 = data[5]
            n=512
            Sw512_FR4 = data[12][-1].flatten() * (1/55509.29780738) / (300/n*0.2)
            x512_FR4 = np.empty_like(data[12][-1])
            for i in range(len(points)):
                x512_FR4[:,i] = np.linspace(0+(points[i] + 1)*L/(2*n),L-(1-points[i])*L/(2*n),n)
            x512_FR4 = x512_FR4.flatten()

            S0 =   np.ones_like(x512_FR4)*0
            f512_FR4 = BL_exact(S0, x512_FR4, fi, A, Swr, Sor, ftime)
            f512_FR4[0] = 1
            e512_L1_FR4 = (sum(abs(f512_FR4-Sw512_FR4))*(300 / n))
            #e512_L1_FR4_avg = (sum(abs(f(x512_avg)-Sw512_avg_FR4))*(300 / n))
            R512_L1_FR4 = math.log(e256_L1_FR4/e512_L1_FR4,2)
            #R512_L1_FR4_avg = math.log(e256_L1_FR4_avg/e512_L1_FR4_avg,2)


        plt.figure(1)
        plt.plot(x8_avg, Sw8_avg, 'tab:pink', x16_avg, Sw16_avg, 'b', x32_avg, Sw32_avg, 'y',
            x64_avg, Sw64_avg, 'g', x128_avg, Sw128_avg, 'c', x256_avg, Sw256_avg, 'm',
            x512_avg, Sw512_avg, 'r')
    #    plt.plot(xD, SwD, 'k')
        plt.grid()
        plt.legend(('8 elements','16 elements', '32 elements', '64 elements', '128 elements',
                    '256 elements', '512 elements'))
        plt.title('Buckley-Leverett Solution Example Bastian - FR 2nd order')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Bastian_FR2_meshes.eps', format='eps')

        plt.figure(2)
        plt.plot(x8_avg, Sw8_avg_FR3, 'tab:pink', x16_avg, Sw16_avg_FR3, 'b', x32_avg, Sw32_avg_FR3, 'y',
            x64_avg, Sw64_avg_FR3, 'g', x128_avg, Sw128_avg_FR3, 'c', x256_avg, Sw256_avg_FR3, 'm',
            x512_avg, Sw512_avg_FR3, 'r')
        #plt.plot(xD, SwD, 'k')
        plt.grid()
        plt.legend(('8 elements','16 elements', '32 elements', '64 elements', '128 elements',
                    '256 elements', '512 elements', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example Bastian - FR 3rd order')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Bastian_FR3_meshes.eps', format='eps')

        plt.figure(3)
        plt.plot(x8_avg, Sw8_avg_FR4, 'tab:pink', x16_avg, Sw16_avg_FR4, 'b', x32_avg, Sw32_avg_FR4, 'y',
            x64_avg, Sw64_avg_FR4, 'g', x128_avg, Sw128_avg_FR4, 'c', x256_avg, Sw256_avg_FR4, 'm',
            x512_avg, Sw512_avg_FR4, 'r')
        #plt.plot(xD, SwD, 'k')
        plt.grid()
        plt.legend(('8 elements','16 elements', '32 elements', '64 elements', '128 elements',
                    '256 elements', '512 elements', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example Bastian - FR 4th order')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Bastian_FR4_meshes.eps', format='eps')

        plt.figure(4)
        x = np.log10([32, 64, 128, 256, 512])
        y = np.log10([e32_L1_FR2, e64_L1_FR2, e128_L1_FR2, e256_L1_FR2, e512_L1_FR2])
        R = np.array([R64_L1_FR2, R128_L1_FR2, R256_L1_FR2, R512_L1_FR2])

        for i in range(len(y)-1):
            plt.plot(x[i:i+2], y[i:i+2])
        plt.legend(('{0:.2f}'.format(R[0]), '{0:.2f}'.format(R[1]), '{0:.2f}'.format(R[2]),
        '{0:.2f}'.format(R[3])))
        plt.plot(x,-x+1,'k')
        plt.grid()
        #plt.legend(('FR-2nd order', 'Reference Line'))
        plt.title('Buckley-Leverett Solution Example Bastian - FR 2nd order')
        plt.ylabel('log$_{10}(E_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Bastian_FR2_converg.eps', format='eps')

        plt.figure(5)
        y3 = np.log10([e32_L1_FR3, e64_L1_FR3, e128_L1_FR3, e256_L1_FR3, e512_L1_FR3])
        R3 = np.array([R64_L1_FR3, R128_L1_FR3, R256_L1_FR3, R512_L1_FR3])
        for i in range(len(y3)-1):
            plt.plot(x[i:i+2], y3[i:i+2])
        plt.legend(('{0:.2f}'.format(R3[0]), '{0:.2f}'.format(R3[1]), '{0:.2f}'.format(R3[2]),
        '{0:.2f}'.format(R3[3])))
        plt.plot(x,-x+1,'k')
        plt.grid()
        #plt.legend(('FR-2nd order', 'Reference Line'))
        plt.title('Buckley-Leverett Solution Example Bastian - FR 3rd order')
        plt.ylabel('log$_{10}(E_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Bastian_FR3_converg.eps', format='eps')

        plt.figure(6)
        y4 = np.log10([e32_L1_FR4, e64_L1_FR4, e128_L1_FR4, e256_L1_FR4, e512_L1_FR4])
        R4 = np.array([R64_L1_FR4, R128_L1_FR4, R256_L1_FR4, R512_L1_FR4])
        for i in range(len(y3)-1):
            plt.plot(x[i:i+2], y4[i:i+2])
        plt.legend(('{0:.2f}'.format(R4[0]), '{0:.2f}'.format(R4[1]), '{0:.2f}'.format(R4[2]),
        '{0:.2f}'.format(R4[3])))
        plt.plot(x,-x+1,'k')
        plt.grid()
        #plt.legend(('FR-2nd order', 'Reference Line'))
        plt.title('Buckley-Leverett Solution Example Bastian - FR 4th order')
        plt.ylabel('log$_{10}(E_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Bastian_FR4_converg.eps', format='eps')

        plt.figure(7)
        plt.plot(x,y,'-bo',x,y3,'-rP', x, y4, '-gv')
        plt.plot(x,-x+1,'k')
        plt.grid()
        plt.legend(('CPR-P1', 'CPR-P2', 'CPR-P3', 'Referência (1$^a$ ordem)'))
        #plt.title('Solução Problema de Bastian - FR convergence')
        plt.ylabel('log$_{10}(E_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.savefig('results/compositional/FR/saturation_BL_Bastian_FR_convergencia.png')

        fig = plt.figure(8)
        plt.plot(x32_avg, Sw32_avg, '-ro', x32_avg, Sw32_avg_FR3, '-gs', x32_avg, Sw32_avg_FR4, '-bv',
            xD, SwD, '-k')
        plt.legend(('FR 2nd order','FR 3rd order', 'FR 4th order', 'Analytical Solution'))

        plt.title('Buckley-Leverett Solution Example Bastian - 32 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Bastian_32_comparison.eps', format='eps')

        fig = plt.figure(9)
        plt.plot(x32_avg, Sw32_avg, '-ro', x32_avg, Sw32_avg_FR3, '-gs', x32_avg, Sw32_avg_FR4, '-bv',
            xD, SwD, '-k')
        plt.legend(('FR 2nd order','FR 3rd order', 'FR 4th order', 'Analytical Solution'))
        plt.xlim(0.55,0.8)
        plt.ylim(0.6,0.8)
        plt.title('Buckley-Leverett Solution Example Bastian - 32 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Bastian_32_comparison_zoom.eps', format='eps')
        import pdb; pdb.set_trace()

        plt.figure(3)
        plt.plot(x64, Sw64, 'r', x64, Sw64_upw, 'g', xD, SwD, 'k', x64, Sw64_FR, 'b')
        plt.legend(('MUSCL','FOUM', 'Analytical Solution', 'FR'))
        plt.title('Buckley-Leverett Solution Example - 64 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_64_comparison.eps', format='eps')

        plt.figure(4)
        plt.plot(x128, Sw128, 'r', x128, Sw128_upw, 'g', xD, SwD, 'k', x128, Sw128_FR, 'b')
        plt.legend(('MUSCL','FOUM', 'Analytical Solution', 'FR'))
        plt.title('Buckley-Leverett Solution Example - 128 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_128_comparison.eps', format='eps')

        plt.figure(5)
        plt.plot( x256, Sw256_upw, 'g', xD, SwD, 'k', x256, Sw256_FR, 'b', x256, Sw256, 'r')
        plt.legend(('FOUM', 'Analytical Solution', 'FR', 'MUSCL'))
        plt.title('Buckley-Leverett Solution Example - 256 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_256_comparison.eps', format='eps')

        plt.figure(6)
        plt.plot(x512, Sw512_upw, 'g', xD, SwD, 'k', x512, Sw512_FR, 'b', x512, Sw512, 'r', x512, Sw512_FR3, 'y')
        plt.legend(('FOUM', 'Analytical Solution', 'FR', 'MUSCL', 'FR (3rd order)'))
        plt.title('Buckley-Leverett Solution Example - 512 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_512_comparison.eps', format='eps')

        plt.figure(7)
        plt.plot(x8, Sw8, 'r', x8, Sw8_upw, 'g', xD, SwD, 'k', x8, Sw8, 'b')
        plt.legend(('MUSCL','FOUM', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - 8 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_8_comparison.eps', format='eps')

        plt.figure(8)
        plt.plot(x16, Sw16, 'r', x16, Sw16_upw, 'g', xD, SwD, 'k', x16, Sw16_FR, 'b')
        plt.legend(('MUSCL','FOUM', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - 16 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_16_comparison.eps', format='eps')

        plt.figure(15)
        x = np.log10(np.array([8,16,32,64,128,256, 512]))
        y_FR = np.log10(np.array([e8_L1_FR2, e16_L1_FR2, e32_L1_FR2, e64_L1_FR2, e128_L1_FR2,
            e256_L1_FR2, e512_L1_FR2]))
        y_MUSCL = np.log10(np.array([e8_L1_MUSCL, e16_L1_MUSCL, e32_L1_MUSCL, e64_L1_MUSCL,
            e128_L1_MUSCL, e256_L1_MUSCL, e512_L1_MUSCL]))
        y_upw = np.log10(np.array([e8_L1_upw, e16_L1_upw, e32_L1_upw, e64_L1_upw, e128_L1_upw,
            e256_L1_upw, e512_L1_upw]))

        y_ref = -x-0.1
        plt.plot(x, y_MUSCL, '-g^', x, y_ref, 'b', x, y_upw, 'y', x, y_FR, '-ro')
        plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{10}({E}_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.legend(('MUSCL-2nd order', 'Reference Line', 'FOU', 'FR'))
        plt.grid()
        plt.savefig('results/compositional/FR_paper/BL_Darlan_L1_convergence_order2' +'.eps', format='eps')

        import pdb; pdb.set_trace()
