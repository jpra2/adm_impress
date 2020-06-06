import numpy as np
import matplotlib.pyplot as plt
import os
import math

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
i=1

krw_end = 0.2
kro_end = 1
nw = 2
no = 2
Swr = 0.2
Sor = 0.35
miu_w = 1
miu_o = 20
k = 0.5
W = 0.1*30.48
T = 0.1*30.48
A = W*T
qw = 0.1*0.32774128
q = qw
alpha = math.pi/2
delta_rho = 0.315
g = 980.665
g_gradient = delta_rho*g/1013216
Sw = np.linspace(Swr, 1-Sor, 500)
Sw1 = np.zeros([1,1])
Sw3 = np.zeros([1,1])
xD3 = np.zeros([1,1])
td = 0.2

S = (Sw - Swr)/(1 - Swr - Sor)
krw = krw_end*S**nw
kro = kro_end*(1-S)**no
fw = (1 + (k*kro*A)*(g_gradient*np.sin(alpha))/(q*miu_o))/((krw*miu_o) + (kro*miu_w))*(krw*miu_o)

p1 = np.polyfit(Sw, fw, 12)
p2 = np.polyder(p1)
diff_fw = np.polyval(p2, Sw)

## Find the saturation shock:
for n in range(1, 500):
    if abs((fw[n])/(Sw[n]-Swr) - diff_fw[n]) < 0.009:
        Swf = Sw[n]
        slope_Swf = (fw[n])/(Sw[n]-Swr)
        a = n

##Water saturation profile of water front:
Sw2 = np.linspace(Swr, Swf, 5)
xD2 = np.linspace(slope_Swf*td,slope_Swf*td,5)

## Water saturation profile before water front:
Sw1 = np.zeros(a-1)
for i in range(0,a-1):
    Sw1[i] = Swr
    xD1 = np.linspace(1, slope_Swf*td,a-1)

## Water saturation profile behind water front:
Sw3 = np.zeros(500-a+1)
xD3 = np.zeros(500-a+1)
for j in range(a-1,500):
    Sw3[j-(a-1)] = Sw[j]
    xD3[j-(a-1)] = diff_fw[j]*td

xD = np.append(xD1,xD2)
xD = np.append(xD,xD3)
SwD = np.append(Sw1,Sw2)
SwD = np.append(SwD,Sw3)


for  arq in arquivos:
    if  arq.startswith(name):
        datas = np.load('flying/results_Buckley_Leverett_vertical_caset_506.npy', allow_pickle=True)
        import pdb; pdb.set_trace()

        for data in datas[8:]:
            Sw = data[5]
            x = np.linspace(0,1,500)
            plt.figure(1)
            plt.plot(x, Sw, 'r', xD, SwD, 'y')
            plt.grid()
            loop = data[0]
            plt.legend(('PADMEC', 'Analytical Solution'))
            plt.title('Buckley-Leverett Solution Example')
            plt.ylabel('Water Saturation')
            plt.xlabel('Dimensionless distance')
            plt.savefig('results/compositional/saturation_w_vert_BL_comparison11.png' )


        '''datas = np.load('flying/results_vertical_case_7855.npy', allow_pickle=True)
        import pdb; pdb.set_trace()
        for data in datas[1:]:
            Sw = data[5]
            x = np.linspace(0,1,500)
            plt.figure(1)
            plt.plot(x, Sw, 'r')
            plt.grid()
            loop = data[0]
            plt.title('Teste')
            plt.ylabel('Water Saturation')
            plt.xlabel('Dimensionless distance')
            plt.savefig('results/compositional/saturation_w_vertical_tcase.png' )'''
