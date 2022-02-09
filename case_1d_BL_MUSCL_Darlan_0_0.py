import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)


xD = np.loadtxt('x_BL_Darlan_semi_analytical.txt')
SwD = np.loadtxt('Sw_BL_Darlan_semi_analytical.txt')
f = interp1d(xD,SwD)

for  arq in arquivos:
    if  arq.startswith(name):
        '''-------------------------MUSCL LLF RESULTS------------------------'''
        #datas = np.load('flying/results_BL_Darlan_8t_875.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_8_MUSCL_43.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw8 = data[5]
            x8 = np.linspace(xD[-1],1,8)
            n = 8
            e8_L1_MUSCL = (sum(abs(f(x8)-Sw8))*(1/n))

        #datas = np.load('flying/results_BL_Darlan_16t_1658.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_16_MUSCL_83.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw16 = data[5]
            x16 = np.linspace(xD[-1],1,16)
            n = 16
            e16_L1_MUSCL = (sum(abs(f(x16)-Sw16))*(1/n))
            R16_L1_MUSCL = math.log(e8_L1_MUSCL/e16_L1_MUSCL,2)

        #datas = np.load('flying/results_BL_Darlan_32t_3213.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_32_MUSCL_175.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw32 = data[5]
            x32 = np.linspace(xD[-1]+1/64,1-1/64,32)
            n = 32
            e32_L1_MUSCL = (sum(abs(f(x32)-Sw32))*(1/n))
            R32_L1_MUSCL = math.log(e16_L1_MUSCL/e32_L1_MUSCL,2)

        #datas = np.load('flying/results_BL_Darlan_64t_1255.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_64_MUSCL_349.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw64 = data[5]
            x64 = np.linspace(xD[-1],1,64)
            n = 64
            e64_L1_MUSCL = (sum(abs(f(x64)-Sw64))*(1/n))
            R64_L1_MUSCL = math.log(e32_L1_MUSCL/e64_L1_MUSCL,2)

        datas = np.load('flying/results_BL_Darlan_128_MUSCL_692.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw128 = data[5]
            x128 = np.linspace(xD[-1],1,128)
            n = 128
            e128_L1_MUSCL = (sum(abs(f(x128)-Sw128))*(1/n))
            R128_L1_MUSCL = math.log(e64_L1_MUSCL/e128_L1_MUSCL,2)

        datas = np.load('flying/results_BL_Darlan_256_MUSCL_1389.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw256 = data[5]
            x256 = np.linspace(xD[-1],1,256)
            n = 256
            e256_L1_MUSCL = (sum(abs(f(x256)-Sw256))*(1/n))
            R256_L1_MUSCL = math.log(e128_L1_MUSCL/e256_L1_MUSCL,2)

        #datas = np.load('flying/results_BL_Darlan_512_MUSCL_2902.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_512_MUSCL_2792.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw512 = data[5]
            P512_MUSCL = data[4]
            Nk512_MUSCL = data[12][0]
            x512 = np.linspace(xD[-1],1,512)
            n = 512
            e512_L1_MUSCL = (sum(abs(f(x512)-Sw512))*(1/n))
            R512_L1_MUSCL = math.log(e256_L1_MUSCL/e512_L1_MUSCL,2)


        '''----------------------------UPWIND RESULTS------------------------'''

        #datas = np.load('flying/results_BL_Darlan_8_upw_1950.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_8_upw_40.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw8_upw = data[5]
            n = 8
            e8_L1_upw = (sum(abs(f(x8)-Sw8_upw))*(1/n))

        #datas = np.load('flying/results_BL_Darlan_16_upw_3492.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_16_upw_82.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw16_upw = data[5]
            n = 16
            e16_L1_upw = (sum(abs(f(x16)-Sw16_upw))*(1/n))
            R16_L1_upw = math.log(e8_L1_upw/e16_L1_upw,2)

        #datas = np.load('flying/results_BL_Darlan_32_upw_6488.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_32_upw_169.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw32_upw = data[5]
            n = 32
            e32_L1_upw = (sum(abs(f(x32)-Sw32_upw))*(1/n))
            R32_L1_upw = math.log(e16_L1_upw/e32_L1_upw,2)

        #datas = np.load('flying/results_BL_Darlan_64_upw_1267.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_64_upw_333.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw64_upw = data[5]
            n = 64
            e64_L1_upw = (sum(abs(f(x64)-Sw64_upw))*(1/n))
            R64_L1_upw = math.log(e32_L1_upw/e64_L1_upw,2)

        #datas = np.load('flying/results_BL_Darlan_128_upw_477.npy', allow_pickle=True)
        #datas = np.load('flying/results_BL_Darlan_128_upw_964.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_128_upw_660.npy', allow_pickle=True)
        #datas = np.load('flying/results_BL_Darlan_128_upw_1930.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw128_upw = data[5]
            n = 128
            e128_L1_upw = (sum(abs(f(x128)-Sw128_upw))*(1/n))
            R128_L1_upw = math.log(e64_L1_upw/e128_L1_upw,2)

        #datas = np.load('flying/results_BL_Darlan_256_upw_1074.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_256_upw_1319.npy', allow_pickle=True)
        #datas = np.load('flying/results_BL_Darlan_256_upw_4837.npy', allow_pickle=True)
        b=0
        for data in datas[datas.shape[0]-1:]:
            Sw256_upw = data[5]
            Nk256_upw = data[12][0]
            n = 256
            e256_L1_upw = (sum(abs(f(x256)-Sw256_upw))*(1/n))
            R256_L1_upw = math.log(e128_L1_upw/e256_L1_upw,2)

        #datas = np.load('flying/results_BL_Darlan_512_upw_2042.npy', allow_pickle=True)
        #datas = np.load('flying/results_BL_Darlan_512_upw_3775.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_512_upw_2633.npy', allow_pickle=True)
        #datas = np.load('flying/results_BL_Darlan_512_upw_37677.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw512_upw = data[5]
            P512_upw = data[4]
            Nk512_upw = data[12][0]
            n = 512
            e512_L1_upw = (sum(abs(f(x512)-Sw512_upw))*(1/n))
            R512_L1_upw = math.log(e256_L1_upw/e512_L1_upw,2)

        '''-------------------FLUX RECONSTRUCTION RESULTS (2nd order)--------------------'''
        #datas = np.load('flying/results_BL_Darlan_8_FR_88.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_8_FR_105.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw8_FR = data[5]
            Nk8_FR = data[12][0]
            n = 8
            e8_L1_FR2 = (sum(abs(f(x8)-Sw8_FR))*(1/n))

        #datas = np.load('flying/results_BL_Darlan_16_FR_176.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_16_FR2_205.npy', allow_pickle=True) #271

        for data in datas[datas.shape[0]-1:]:
            Sw16_FR = data[5]
            Nk16_FR = data[12][0]
            n = 16
            e16_L1_FR2 = (sum(abs(f(x16)-Sw16_FR))*(1/n))
            R16_L1_FR2 = math.log(e8_L1_FR2/e16_L1_FR2,2)

        #datas = np.load('flying/results_BL_Darlan_32_FR_354.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_32_FR_414.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw32_FR = data[5]
            x32 = np.linspace(xD[-1],1,32)
            Nk32_FR = data[12][0]
            n = 32
            e32_L1_FR2 = (sum(abs(f(x32)-Sw32_FR))*(1/n))
            R32_L1_FR2 = math.log(e16_L1_FR2/e32_L1_FR2,2)

        #datas = np.load('flying/results_BL_Darlan_64_FR_713.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_64_FR2_1080.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw64_FR = data[5]
            Nk64_FR = data[12][0]
            n = 64
            e64_L1_FR2 = (sum(abs(f(x64)-Sw64_FR))*(1/n))
            R64_L1_FR2 = math.log(e32_L1_FR2/e64_L1_FR2,2)

        #datas = np.load('flying/results_BL_Darlan_128_FR_1984.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_128_FR2_1679.npy', allow_pickle=True)
        #datas = np.load('flying/results_BL_Darlan_128_FR2_2425.npy', allow_pickle=True)
        for data in datas[datas.shape[0]-1:]:
            Sw128_FR = data[5]
            Nk128_FR = data[12][0]
            n = 128
            e128_L1_FR2 = (sum(abs(f(x128)-Sw128_FR))*(1/n))
            R128_L1_FR2 = math.log(e64_L1_FR2/e128_L1_FR2,2)

        datas = np.load('flying/results_BL_Darlan_256_FR_4003.npy', allow_pickle=True)
        #datas = np.load('flying/results_BL_Darlan_256_FR_5038.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw256_FR = data[5]
            Nk256_FR = data[12][0]
            n = 256
            e256_L1_FR2 = (sum(abs(f(x256)-Sw256_FR))*(1/n))
            R256_L1_FR2 = math.log(e128_L1_FR2/e256_L1_FR2,2)

        #datas = np.load('flying/results_BL_Darlan_512_FR_8019.npy', allow_pickle=True)
        datas = np.load('flying/results_BL_Darlan_512_FR_6921.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR = data[5]
            P512_FR = data[4]
            x512_FR = np.linspace(xD[-1],1,512)
            Nk512_FR = data[12][0]
            n = 512
            e512_L1_FR2 = (sum(abs(f(x512_FR)-Sw512_FR))*(1/n))
            R512_L1_FR2 = math.log(e256_L1_FR2/e512_L1_FR2,2)

        datas = np.load('flying/results_BL_Darlan_512_FR2_u1_8415.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR_u1 = data[5]

        datas = np.load('flying/results_BL_Darlan_512_FR2_vk_2952.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR_vk = data[5]

        datas = np.load('flying/results_BL_Darlan_512_FR2_vk_2954.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR_vk2 = data[5]

        datas = np.load('flying/results_BL_Darlan_512_FR2_vk_2953.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR_vk3 = data[5]

        datas = np.load('flying/results_BL_Darlan_512_FR2_vk8_2953.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR_vk4 = data[5]

        datas = np.load('flying/results_BL_Darlan_1024_FR2_4635.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw1024_FR = data[5]

        datas = np.load('flying/results_BL_Darlan_1024_FR2_vk_5917.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw1024_FR_vk = data[5]



        '''-----------------Flux Reconsctruction (3rd order)-----------------'''
        datas = np.load('flying/results_BL_Darlan_8_FR3_170.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw8_FR3 = data[5]
            Nk8_FR3 = data[12][0]
            n = 8
            e8_L1_FR3 = (sum(abs(f(x8)-Sw8_FR3))*(1/n))

        datas = np.load('flying/results_BL_Darlan_16_FR3_340.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw16_FR3 = data[5]
            Nk16_FR3 = data[12][0]
            n = 16
            e16_L1_FR3 = (sum(abs(f(x16)-Sw16_FR3))*(1/n))
            R16_L1_FR3 = math.log(e8_L1_FR3/e16_L1_FR3,2)

        datas = np.load('flying/results_BL_Darlan_32_FR3_687.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw32_FR3 = data[5]
            Nk32_FR3 = data[12][0]
            n = 32
            e32_L1_FR3 = (sum(abs(f(x32)-Sw32_FR3))*(1/n))
            R32_L1_FR3 = math.log(e16_L1_FR3/e32_L1_FR3,2)

        datas = np.load('flying/results_BL_Darlan_64_FR3_1380.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw64_FR3 = data[5]
            Nk64_FR3 = data[12][0]
            n = 64
            e64_L1_FR3 = (sum(abs(f(x64)-Sw64_FR3))*(1/n))
            R64_L1_FR3 = math.log(e32_L1_FR3/e64_L1_FR3,2)

        datas = np.load('flying/results_BL_Darlan_128_FR3_2772.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw128_FR3 = data[5]
            Nk128_FR3 = data[12][0]
            n = 128
            e128_L1_FR3 = (sum(abs(f(x128)-Sw128_FR3))*(1/n))
            R128_L1_FR3 = math.log(e64_L1_FR3/e128_L1_FR3,2)

        datas = np.load('flying/results_BL_Darlan_256_FR3_5568.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw256_FR3 = data[5]
            Nk256_FR3 = data[12][0]
            n = 256
            e256_L1_FR3 = (sum(abs(f(x256)-Sw256_FR3))*(1/n))
            R256_L1_FR3 = math.log(e128_L1_FR3/e256_L1_FR3,2)

        datas = np.load('flying/results_BL_Darlan_512_FR3_11275.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR3 = data[5]
            Nk512_FR3 = data[12][0]
            n = 512
            e512_L1_FR3 = (sum(abs(f(x512)-Sw512_FR3))*(1/n))
            R512_L1_FR3 = math.log(e256_L1_FR3/e512_L1_FR3,2)

        datas = np.load('flying/results_BL_Darlan_512_FR3_u1_11531.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR3_u1 = data[5]

        datas = np.load('flying/results_BL_Darlan_512_FR3_vk_4145.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR3_vk = data[5]

        datas = np.load('flying/results_BL_Darlan_1024_FR3_7383.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw1024_FR3 = data[5]

        datas = np.load('flying/results_BL_Darlan_1024_FR3_vk_8309.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw1024_FR3_vk = data[5]


        '''-----------------Flux Reconsctruction (3rd order)-----------------'''
        datas = np.load('flying/results_BL_Darlan_8_FR4_237.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw8_FR4 = data[5]
            Nk8_FR4 = data[12][0]
            n = 8
            e8_L1_FR4 = (sum(abs(f(x8)-Sw8_FR4))*(1/n))


        datas = np.load('flying/results_BL_Darlan_16_FR4_477.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw16_FR4 = data[5]
            Nk16_FR4 = data[12][0]
            n = 16
            e16_L1_FR4 = (sum(abs(f(x16)-Sw16_FR4))*(1/n))
            R16_L1_FR4 = math.log(e8_L1_FR4/e16_L1_FR4,2)

        datas = np.load('flying/results_BL_Darlan_32_FR4_959.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw32_FR4 = data[5]
            Nk32_FR4 = data[12][0]
            n = 32
            e32_L1_FR4 = (sum(abs(f(x32)-Sw32_FR4))*(1/n))
            R32_L1_FR4 = math.log(e16_L1_FR4/e32_L1_FR4,2)

        datas = np.load('flying/results_BL_Darlan_64_FR4_1926.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw64_FR4 = data[5]
            Nk64_FR4 = data[12][0]
            n = 64
            e64_L1_FR4 = (sum(abs(f(x64)-Sw64_FR4))*(1/n))
            R64_L1_FR4 = math.log(e32_L1_FR3/e64_L1_FR4,2)

        datas = np.load('flying/results_BL_Darlan_128_FR4_3873.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw128_FR4 = data[5]
            Nk128_FR4 = data[12][0]
            n = 128
            e128_L1_FR4 = (sum(abs(f(x128)-Sw128_FR4))*(1/n))
            R128_L1_FR4 = math.log(e64_L1_FR4/e128_L1_FR4,2)



        datas = np.load('flying/results_BL_Darlan_512_FR4_15677.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR4 = data[5]
            Nk512_FR4 = data[12][0]
            n = 512
            #e512_L1_FR4 = (sum(abs(f(x512)-Sw512_FR4))*(1/n))
            #R512_L1_FR4 = math.log(e256_L1_FR4/e512_L1_FR4,2)

        datas = np.load('flying/results_BL_Darlan_512_FR4_u1_5551.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR4_u1 = data[5]

        datas = np.load('flying/results_BL_Darlan_512_FR4_vk_5466.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw512_FR4_vk = data[5]

        datas = np.load('flying/results_BL_Darlan_1024_FR4_10331.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw1024_FR4 = data[5]

        datas = np.load('flying/results_BL_Darlan_1024_FR4_vk_11035.npy', allow_pickle=True)

        for data in datas[datas.shape[0]-1:]:
            Sw1024_FR4_vk = data[5]




        plt.figure(1)
        plt.plot(x8, Sw8_FR, 'm', x16, Sw16_FR, 'c', x32, Sw32_FR, 'r',
                x64, Sw64_FR, 'b', x128, Sw128_FR, 'y', x256, Sw256_FR, 'g', x512, Sw512_FR,
                'm', xD, SwD, 'k')
        plt.grid()
        plt.legend(('8 elements', '16 elements','32 elements','64 elements',
                    '128 elements', '256 elements', '512 elements', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - FR 2nd order')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_FR2_meshes.png' )

        plt.figure(10)
        plt.plot(x8, Sw8_FR3, 'm', x16, Sw16_FR3, 'c', x32, Sw32_FR3, 'r',
                x64, Sw64_FR3, 'b', x128, Sw128_FR3, 'y', x256, Sw256_FR3, 'g',
                x512, Sw512_FR3, 'm', xD, SwD, 'k')
        plt.grid()
        plt.legend(('8 elements', '16 elements', '32 elements','64 elements',
                    '128 elements', '256 elements', '512 elements', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - FR 3rd order')
        plt.ylabel('Water Saturation')
        plt.xlim(0.60,0.8)
        plt.ylim(0.60,0.71)
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_FR3_meshes.png' )

        plt.figure(2)
        x32 = np.linspace(xD[-1]+1/64,1-1/64,32)
        plt.plot(x32, Sw32, '-r<', x32, Sw32_upw, '-gP', x32, Sw32_FR, '-mo', xD, SwD, 'k', mfc='none')#
                # x32, Sw32_FR4, '-c<', xD, SwD, 'k')
        plt.legend(('MUSCL','FOU', 'FR P1', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - 32 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_32_comparison.png' )

        plt.figure(29)
        x64 = np.linspace(xD[-1]+1/128,1-1/128,64)
        plt.plot(x64, Sw64, '-r<', x64, Sw64_upw, '-gP', x64, Sw64_FR, '-mo', xD, SwD, 'k', mfc='none')#
                # x32, Sw32_FR4, '-c<', xD, SwD, 'k')
        plt.legend(('MUSCL','FOU', 'FR P1', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - 64 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_64_comparison.png' )

        plt.figure(3)
        plt.plot(x32, Sw32, 'r', x64, Sw64, 'b', x128, Sw128, 'y', x256, Sw256, 'g', x512, Sw512,
         'm', xD, SwD, 'k')
        plt.grid()
        plt.legend(('32 elements','64 elements', '128 elements', '256 elements', '512 elements', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - MUSCL')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR_paper/saturation_BL_Darlan_MUSCL_meshes.png' )

        plt.figure(4)
        x1024 = np.linspace(0,1,1024)
        #plt.plot(x512, Sw512_FR4_u1, 'r')
        plt.plot(x1024, Sw1024_FR4_vk, 'mo', xD, SwD, 'k', mfc='none')
        plt.plot(x1024, Sw1024_FR3_vk, 'gv', mfc='none')
        plt.plot(x1024, Sw1024_FR_vk, 'ys', mfc='none')
        plt.xlim(0.65,0.78)
        plt.ylim(0.65,0.7)
        plt.legend(('FR P3 - MLP vk', 'Analytical Solution', 'FR P2 - MLP vk', 'FR P1 - MLP vk'))
        plt.title('Buckley-Leverett Solution Example - 1024 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_1024_vk_comparison.png', format='png')

        plt.figure(16)
        x1024 = np.linspace(0,1,1024)
        plt.plot(x512, Sw512_FR_vk3, 'mo', xD, SwD, '-k', mfc='none')
        plt.plot(x512, Sw512_FR_u1, 'gv', mfc='none')
        plt.plot(x512, Sw512_FR, 'ys', mfc='none')
        plt.xlim(0.65,0.78)
        plt.ylim(0.65,0.7)
        plt.legend(('FR P1-MLP vk', 'Analytical Solution', 'FR P1-MLPu1 mod ', 'FR P1-MLPu1 '))
        plt.title('Buckley-Leverett Solution Example - 512 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_512_vk_P1_comparison.png', format='png')

        plt.figure(22)
        x1024 = np.linspace(0,1,1024)
        #plt.plot(x512, Sw512_FR4_u1, 'r')
        plt.plot(x512, Sw512_FR3_vk, 'mo', xD, SwD, '-k', mfc='none')
        plt.plot(x512, Sw512_FR3_u1, 'gv', mfc='none')
        plt.plot(x512, Sw512_FR3, 'ys', mfc='none')
        plt.xlim(0.65,0.78)
        plt.ylim(0.65,0.7)
        plt.legend(('FR P1-MLP vk', 'Analytical Solution', 'FR P1-MLPu1 mod ', 'FR P1-MLPu1 '))
        plt.title('Buckley-Leverett Solution Example - 512 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_512_vk_P3_comparison.png', format='png')

        plt.figure(20)
        plt.plot(x512, Sw512_FR4_vk, 'mo', xD, SwD, '-k', mfc='none')
        plt.plot(x512, Sw512_FR4_u1, 'gv', mfc='none')
        plt.plot(x512, Sw512_FR4, 'ys', mfc='none')
        plt.xlim(0.65,0.78)
        plt.ylim(0.65,0.7)
        plt.legend(('FR P1-MLP vk', 'Analytical Solution', 'FR P1-MLPu1 mod ', 'FR P1-MLPu1 '))
        plt.title('Buckley-Leverett Solution Example - 512 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_512_vk_P4_comparison.png', format='png')

        plt.figure(30)
        plt.plot(x128, Sw128, '-r<', x128, Sw128_upw, '-gP', x128, Sw128_FR, '-mo', xD, SwD, 'k', mfc='none')#, '-mo',  x32, Sw32_FR3, '-bs',
                # x32, Sw32_FR4, '-c<', xD, SwD, 'k')
        plt.legend(('MUSCL','FOU', 'FR P1', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - 128 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_128_comparison.png' )

        plt.figure(14)
        plt.plot(x128, Sw128, '-r<', x128, Sw128_upw, '-gP', x128, Sw128_FR, '-mo', xD, SwD, 'k')#, '-mo',  x32, Sw32_FR3, '-bs',
                # x32, Sw32_FR4, '-c<', xD, SwD, 'k')
        plt.legend(('MUSCL','FOU', 'FR P1', 'Analytical Solution'))
        plt.title('Buckley-Leverett Solution Example - 128 elements')
        plt.ylabel('Water Saturation')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/saturation_BL_Darlan_128_comparison_all.png' )

        plt.figure(15)
        plt.plot(x512, P512_MUSCL, '-r', x512, P512_upw, 'b', x512, P512_FR, 'm')
        plt.legend(('MUSCL','FOUM', 'FR P1'))
        plt.title('Buckley-Leverett Pressure')
        plt.ylabel('Pressure Solver')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/FR/pressure_BL_Darlan_FR.png' )


        plt.figure(26)
        x = np.log10(np.array([8,16,32,64,128,256, 512]))
        y_FR = np.log10(np.array([e8_L1_FR2, e16_L1_FR2, e32_L1_FR2, e64_L1_FR2, e128_L1_FR2,
            e256_L1_FR2, e512_L1_FR2]))
        y_MUSCL = np.log10(np.array([e8_L1_MUSCL, e16_L1_MUSCL, e32_L1_MUSCL, e64_L1_MUSCL,
            e128_L1_MUSCL, e256_L1_MUSCL, e512_L1_MUSCL]))
        y_upw = np.log10(np.array([e8_L1_upw, e16_L1_upw, e32_L1_upw, e64_L1_upw, e128_L1_upw,
            e256_L1_upw, e512_L1_upw]))

        y_ref = -x-0.8
        plt.plot(x, y_MUSCL, '-g^', x, y_ref, 'b', x, y_upw, '-ys', x, y_FR, '-ro', mfc='none')
        plt.title('Convergence rate - L1 norm')
        plt.ylabel('$log_{10}({E}_{L_1})$')
        plt.xlabel('$log_{10}(N)$')
        plt.legend(('MUSCL-2nd order', 'Reference Line', 'FOU', 'FR P1'))
        plt.grid()
        plt.savefig('results/compositional/FR/BL_Darlan_L1_convergence_order2.png')

        import pdb; pdb.set_trace()
