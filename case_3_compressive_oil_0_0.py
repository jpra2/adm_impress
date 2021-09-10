import numpy as np
import matplotlib.pyplot as plt
import os
import math

flying = 'flying'
name = 'results'
arquivos = os.listdir(flying)
i=1
t = np.array([274752, 432000., 500000, 589248.0, 745632., 902880, 1059264.0, 1216512., 1372896.0, 1530144.0])

por = 0.2
k = 500
miu = 0.2489
L = 2000
Pi = 2000
Pe = 1900
Sw = 0.2
So = 0.8
C = 1.04e-5
Cf = 5e-4
Ct = C+Cf
m = 100
x = np.zeros(m)
xd = np.zeros(m)
Pd = np.zeros([len(t),m])
P = np.zeros([len(t),m])
data = np.zeros((m,2))
tds = t/( 0.2 * 2.498e-4 * 7.252e-8 * 609.6**2/4.934617e-13)
j=0
for td in tds:
    for i in range(0,m):
        x[0] = 0
        if i<m-1:
            x[i+1] = x[i] + L/(m-1)
        xd[i] = x[i]/L
        b = 0
        for n in range(1,100):
            lam = 1/2*(2*n-1)*math.pi
            a = 2/lam*np.exp(-(lam**2)*td)*(np.sin(lam*xd[i]))
            b = b+a
        Pd[j,i] = b
        P[j,i] = Pd[j,i]*(Pi-Pe) + Pe
    j+=1

p_ans = P *0.006894736842105263#* 6894.75729/1e6
x_ans = xd

for arq in arquivos:
    if  arq.startswith(name):
        '''
        x_ans = np.array([0.997821, 0.955701, 0.981845, 0.969015, 0.936093, 0.915033, 0.897846, 0.878238, 0.858388, 0.834907, 0.801743,
                          0.777052, 0.748971, 0.724038, 0.696926, 0.668119, 0.635924, 0.601549, 0.572259, 0.54902, 0.520697, 0.47543,
                          0.398935, 0.29436, 0.196078, 0.136771, 0.0723796, 0, 0.033406])
        p_ans = np.array([1986.39, 1986.04, 1986.39, 1986.39, 1985.86, 1985.63, 1985.51, 1984.93, 1984.41, 1983.54, 1982.55, 1981.62,
                          1980.22, 1979.12, 1977.66, 1975.97, 1973.94, 1971.55, 1969.17, 1967.25, 1965.04, 1960.79, 1952.94, 1940.9,
                          1928.27, 1920.54, 1911.46, 1900, 1905.88])'''


        #datas = np.load('flying/results_compressive_oil_case_FR2_148.npy', allow_pickle=True)
        datas = np.load('flying/results_compressive_oil_case_FR_127.npy', allow_pickle=True)
        for data in datas[1:]:
            pressure_FR = data[4]/1e6

        datas = np.load('flying/results_compressive_oil_case_32_FR4_36.npy', allow_pickle=True)
        for data in datas[1:]:
            pressure_FR_32 = data[4]/1e6
            Ft = data[13]
            x32 = np.linspace(0, 1, 32)

        datas = np.load('flying/results_compressive_oil_case_upw_2749.npy', allow_pickle=True)
        for data in datas[1:]:
            pressure = data[4]/1e6
            Ft_upw = data[13]
            x = np.linspace(0, 1, 100)

        datas = np.load('flying/results_compressive_oil_case_MUSCL_2749.npy', allow_pickle=True)
        for data in datas[1:]:
            pressure_MUSCL = data[4]/1e6
            Ft_MUSCL = data[13]


        fig2 = plt.figure(3)
        plt.grid(fig2)
        plt.plot(x, pressure_FR, 'ro', mfc= 'none')
        plt.plot(x, pressure_MUSCL, 'bo',  mfc= 'none')
        plt.plot(x_ans, p_ans[0], 'k')
        plt.legend(('CPR-3$^a$ ordem', 'MUSCL', 'Solução de Referência'))
        plt.ylabel('Nk [mol]')
        plt.xlabel('Distância Adimensional')
        plt.title('Resultados com malha 100x1x1 - t = 5 dias')
        plt.savefig('results/compositional/FR/pressure_comp_oil_5days_FR3' + '.png')

        fig3 = plt.figure(4)
        plt.grid(fig3)
        plt.plot(x32, pressure_FR_32, 'ro', mfc= 'none')
        #plt.plot(x, pressure, 'bP', mfc= 'none')
        plt.plot(x_ans, p_ans[0], 'k')
        plt.legend(('CPR-3$^a$ ordem - 32x1x1', 'Solução de Referência'))
        plt.ylabel('Nk [mol]')
        plt.xlabel('Distância Adimensional')
        plt.title('Resultados para t = 5 dias')
        plt.savefig('results/compositional/FR/pressure_comp_oil_5days_FR3_32' + '.png')
        import pdb; pdb.set_trace()



'''
        p1 = (datas[1][4] - 13.78951458E6 * np.ones(100))/6894.76
        p2 = (datas[2][4] - 13.78951458E6 * np.ones(100))/6894.76
        p3 = (datas[3][4] - 13.78951458E6 * np.ones(100))/6894.76
        p4 = (datas[8][4] - 13.78951458E6 * np.ones(100))/6894.76

        t1 = datas[1][3]
        t2 = datas[2][3]
        t3 = datas[3][3]
        t4 = datas[8][3]
        x = np.linspace(0,1,100)
        plt.figure(1)
        plt.plot(x, p1, 'r', x, p2, 'b', x, p3, 'g', x, p4, 'y')
        plt.grid()
        plt.legend(('%f seg' %t1, '%f seg' %t2, '%f seg' %t3, 'solução'))
        plt.ylabel('Pressure Drop (psi)')
        plt.xlabel('Dimensionless distance')
        plt.savefig('results/compositional/pressure_mcN_more3' + '.png')
        import pdb; pdb.set_trace()'''
