import numpy as np
import math
import matplotlib.pyplot as plt

''' Tentativa de implementar a solução analitica de Lauwerier '''

# Dados de entrada
nx = 200 #blocos da malha
delta_x = 1.524
delta_y = 3.048
delta_z = 3.048
x = np.arange(nx)*delta_x + delta_x

temp_inicial = 288.150
pressao_inicial = 6894.757 #unidade
porosidade = 0.35
Permeability = 9.8692326671601e-12
Cap_rock_thermal_conductivity = 218.00 #unidade
Rock_heat_capacity = 2350.000 #unidade
Heat_capacity_ratio = 1.0

temp_injecao = 366.556
taxa_injecao = 5.465 #unidade
Producer_BHP = 6894.757 #unidade



''' Calculo solucao analitica de Lauwerier de forma adimensional - parametros dados '''
L = 855.0
n = 600
delta = L/n
x = np.arange(n)*delta + delta
xd = np.zeros_like(x)
xd = 0.016*x
U = np.zeros_like(x)
teta = 1.14
td = np.zeros_like(x)
T1_T0 = np.zeros_like(U)

#time_step = np.array([1, 2, 4, 8, 16, 32])
time_step = np.array([5, 9, 18, 32])

for k in time_step:
    
    td[:] = k * teta
    U[(td - xd)>0] = 1

    aux = xd[U==1]/(2*(teta*(td[U==1] - xd[U==1]))**0.5)

    for i in range(aux.size):
        T1_T0[i] = math.erfc(aux[i])

    #print(T1_T0)
    plt.plot(xd/teta, T1_T0)

plt.grid()
plt.axis([0, 12, 0, 1.1])
plt.xlabel("Distancia adimensional")
plt.ylabel("Temperatura adimensional")
plt.title("Problema de Lauwerier")
plt.xticks(np.arange(0, max(xd/teta)+1, 1.0))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.show()

''' Calculo solucao analitica de Lauwerier de forma adimensional - Varavei analitica ok'''
L = 885.0 # Duvida em relação com comprimento L
n = 600
delta = L/n
x = np.arange(n)*delta + delta
xd = np.zeros_like(x)

Cap_rock_thermal_conductivity = 2.52314815 #SI λs ou kc
ht = 3.048 
Vf = 18e-6 
ro_f = 997.0
Cf = 4184.0
teste = 4*Cap_rock_thermal_conductivity/(ht*ht*Vf*ro_f*Cf) #0.01446816583950218
xd = teste*x

U = np.zeros_like(x)
teta = 1.0
td = np.zeros_like(x)
T1_T0 = np.zeros_like(U)

#time_step = np.array([1, 2, 4, 8, 16, 32])
time_step = np.array([5, 9, 18, 32])

for k in time_step:
    
    td[:] = k * teta
    U[(td - xd)>0] = 1

    aux = xd[U==1]/(2*(teta*(td[U==1] - xd[U==1]))**0.5)

    for i in range(aux.size):
        T1_T0[i] = math.erfc(aux[i])

    #print(T1_T0)
    plt.plot(xd/teta, T1_T0)

plt.grid()
plt.axis([0, 12, 0, 1.1])
plt.xlabel("Distancia adimensional")
plt.ylabel("Temperatura adimensional")
plt.title("Problema de Lauwerier")
plt.xticks(np.arange(0, max(xd/teta)+1, 1.0))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.show()

