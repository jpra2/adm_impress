import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import InterpolatedUnivariateSpline

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
mi_w = 1
mi_o = 20
Sw = np.linspace(Swr, 1-Sor, 100000)
td = 0.2

'CALCULATION'
S = (Sw - Swr)/(1 - Swr - Sor)
krw = krw_end*S**nw
kro = kro_end*(1-S)**no
fw = (krw*mi_o) / ((krw*mi_o)+(kro*mi_w))

#p1 = np.polyfit(Sw, fw, 12)
#p1 = interp1d(Sw,fw)
p1 = InterpolatedUnivariateSpline(Sw,fw)
p2 = InterpolatedUnivariateSpline.derivative(p1)
#p2 = np.polyder(p1)
#diff_fw = np.polyval(p2, Sw)
diff_fw = p2(Sw)

## Find the saturation shock:
for n in range(1, 100000):
    if abs((fw[n] - 0)/(Sw[n]-Swr) - diff_fw[n]) < 0.0009:
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
Sw3 = np.zeros(100000-a)
xD3 = np.zeros(100000-a)
for j in range(a,100000):
    Sw3[j-a] = Sw[j]
    xD3[j-a] = diff_fw[j]*td

#Sw3[xD3<0] = Sw3[xD3<0][0]
#xD3[xD3<0] = 0
xD = np.append(xD1,xD2)
xD = np.append(xD,xD3)
SwD = np.append(Sw1,Sw2)
SwD = np.append(SwD,Sw3)

for  arq in arquivos:
    if  arq.startswith(name):
        #datas = np.load('flying/results_Buckley_Leverett_caset_602.npy', allow_pickle=True)
        datas = np.load('flying/results_Buckley_Leverett_case_IMPEC_523.npy', allow_pickle=True)
        for data in datas[7:]:
            Sw_IMPEC = data[5]
            x = np.linspace(0,1,500)

        datas = np.load('flying/results_Buckley_Leverett_case_IMPSAT_1148.npy', allow_pickle=True)
        for data in datas[-1:]:
            Sw_IMPSAT = data[5]
            x = np.linspace(0,1,500)
            plt.figure(1)
            plt.plot(x, Sw_IMPEC, 'rs', xD, SwD, 'yv', x, Sw_IMPSAT, 'bo', mfc = 'none')
            plt.grid()
            loop = data[0]
            plt.legend(('IMPEC', 'Analytical Solution', 'IMPSAT'))
            plt.title('Buckley-Leverett Solution Example')
            plt.ylabel('Water Saturation')
            plt.xlabel('Dimensionless distance')
            plt.savefig('results/compositional/saturation_w_BL_comparison_IMPEC_IMPSAT.png')
        import pdb; pdb.set_trace()
