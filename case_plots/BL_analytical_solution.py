'Function created to plot analytical solution for the water saturation profile  \
for the buckley leverett problem'
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

'DATA ENTRY'
krw_end = 1
kro_end = 1
nw = 2
no = 2
Swr = 0.
Sor = 0.
mi_w = 1e-3
mi_o = 1e-3
Sw = np.linspace(Swr, 1-Sor, 100000)
td = 0.648

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


with open('Sw_BL_FR_analytical.txt', 'w') as f:
    for item in SwD:
        f.write("%s\n" % item)

with open('x_BL_FR_analytical.txt', 'w') as f:
    for item in xD:
        f.write("%s\n" % item)
import pdb; pdb.set_trace()
