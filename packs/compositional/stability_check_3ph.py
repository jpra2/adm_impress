"""Check stability of a thermodynamic equilibrium."""
import numpy as np
import math
#import thermo
from scipy.misc import derivative
from .equation_of_state import PengRobinson
import matplotlib.pyplot as plt
from ..utils import constants as ctes
from .stability_check import StabilityCheck as StabilityCheck_2ph
import time
## Encontrar os pontos estacionarios. Estes correspondem aos pontos nos quais a derivada de g com respeito a Y é 0
## Todas as equações foram confirmadas pelo livro do Dandekar e a biblioteca do thermo

class StabilityCheck(StabilityCheck_2ph):
    """Check for stability of a thermodynamic equilibrium and returns the
    equilibrium phase compositions (perform the flash calculation)."""

    def __init__(self, P, T):
        super().__init__(P,T)
        self.A = np.empty(len(P))     # 3 phases
        #self.x = np.empty([ctes.n_components, len(P)])
        #self.y = np.empty_like(self.x)
        self.a = np.empty_like(self.x)         # 3 phases


        #StabilityCheck.TPD(self)

    def run_init(self, P, z, pflash = True, ponteiro_flash = [], ksi_W=[], rho_W=[]):
        self.z = np.empty_like(z)
        self.z[0,:] = z[-1,:]
        self.z[1:,:] = z[0:-1,:]

        self.P = np.copy(P)
        if np.sum(pflash,dtype=bool)==True:
            ponteiro_flash = np.ones(len(P), dtype = bool)

        self.vapor_pressure_pure_substancies()
        self.Kwilson = self.K.copy()
        self.pressure_water_saturation()
        self.equilibrium_ratio_aqueous(z)


        ponteiro_flash = np.zeros(len(self.P), dtype = bool)
        dir_flash = np.argwhere(z <= 0)
        ponteiro_flash[dir_flash[:,1]] = True

        #sp = np.zeros(5) # In the case of 5 trial phases in stability analysis

        if not pflash and any(~ponteiro_flash) and ctes.Nc>1:
            sp, Kvalue = self.Stability_2phase(self.z, np.copy(~ponteiro_flash))
            sp = np.round(sp, 14)
            #import pdb; pdb.set_trace()
            ponteiro_aux = ponteiro_flash[~ponteiro_flash]
            ponteiro_aux[(sp>1).sum(axis=0,dtype=bool)] = True #os que devem passar para o calculo de flash
            ponteiro_flash[~ponteiro_flash] = ponteiro_aux
            self.K[:,sp[4]>1] = self.Kw[:,sp[4]>1]

        self.ksi_L, self.ksi_V, self.rho_L, self.rho_V = \
            self.molar_properties(np.ones_like(ponteiro_flash, dtype=bool)) # cálculo do flash bifásico

        self.y[:,(self.V<0) + (self.V>1)] = z[:,(self.L>1) + (self.V>1)]
        self.x[:,(self.V<0) + (self.V>1)] = z[:,(self.L>1) + (self.V>1)]
        self.L[self.L>1] = 1
        self.L[self.L<0] = 0
        self.V = 1 - self.L


        ponteiro_flash_3phase = np.zeros(len(self.P), dtype = bool)
        ponteiro_flash_3phase[(self.L != 1) & (self.V != 1)] = True
        ponteiro_flash_3phase2 = ponteiro_flash_3phase.copy()

        sp2, Kvalue2 = self.Stability_3phase(self.x, np.copy(ponteiro_flash_3phase))
        sp2 = np.round(sp2, 8)
        ponteiro_aux = ~ponteiro_flash_3phase[ponteiro_flash_3phase]
        #ponteiro_aux = np.zeros(len(self.P), dtype = bool)
        ponteiro_aux[(sp2>1).sum(axis=0,dtype=bool)] = True
        ponteiro_flash_3phase[ponteiro_flash_3phase] = ponteiro_aux
        #ponteiro_flash_3phase = ponteiro_aux

        #ponteiro_flash_3phase2 = ponteiro_flash_3phase[ponteiro_flash_3phase]
        #ponteiro_flash_3phase2 = np.ones(len(self.P), dtype = bool)

        try:
            ponteiro_flash_3phase2[(sp2>1).sum(axis=0,dtype=bool)] = False
        except: import pdb; pdb.set_trace()
        #ponteiro_flash_3phase[ponteiro_flash_3phase] = ponteiro_flash_3phase2

        sp3, Kvalue3 = self.Stability_3phase(self.y, np.copy(ponteiro_flash_3phase2))
        sp3 = np.round(sp3, 8)

        ponteiro_aux2 = ~ponteiro_flash_3phase2[ponteiro_flash_3phase2]
        #ponteiro_aux2 = np.zeros(len(self.P), dtype = bool)
        ponteiro_aux2[(sp3>1).sum(axis=0,dtype=bool)] = True
        ponteiro_flash_3phase2[ponteiro_flash_3phase2] = ponteiro_aux2
        #ponteiro_flash_3phase2 = ponteiro_aux2

        ponteiro_flash_3phase = ponteiro_flash_3phase + ponteiro_flash_3phase2
        #ponteiro_flash_3phase[ponteiro_flash_3phase] = ponteiro_flash_3phase2

        '''
        ponteiro_flash_3phase2 = ~ponteiro_flash_3phase[ponteiro_flash_3phase]
        ponteiro_flash_3phase2[(sp3>1).sum(axis=0,dtype=bool)] = True

        #ponteiro_flash_3phase[ponteiro_flash_3phase] = ponteiro_flash_3phase2
        ponteiro_flash_3phase[ponteiro_flash_3phase2] = ponteiro_flash_3phase2
        '''

        self.A[~ponteiro_flash_3phase] = 0.0
        self.a[:, ~ponteiro_flash_3phase] = 0.0
        self.K_A = self.K[:, ponteiro_flash_3phase]
        self.K_V = self.Kwilson[:, ponteiro_flash_3phase]

        self.ksi_L, self.ksi_V, self.ksi_A, self.rho_L, self.rho_V, self.rho_A = \
            self.molar_properties_3phase(z, ponteiro_flash_3phase)

        self.y[:,(self.V<0) + (self.V>1)] = z[:,(self.L>1) + (self.V>1)]
        self.x[:,(self.V<0) + (self.V>1)] = z[:,(self.L>1) + (self.V>1)]
        self.a[:,(self.A<0) + (self.A>1)] = z[:,(self.A>1) + (self.A>1)]
        self.L[self.L>1] = 1
        self.L[self.L<0] = 0
        self.A[self.A>1] = 1
        self.A[self.A<0] = 0
        self.V = 1 - self.L - self.A



        x_out = self.water_last(self.x)
        y_out = self.water_last(self.y)
        a_out = self.water_last(self.a)

        self.organize_outputs(x_out, y_out, a_out)
        return self.L, self.V, self.A, self.xkj, self.Csi_j, self.rho_j

    def organize_outputs(self, x_out, y_out, a_out):
        self.xkj[:,0,:] = x_out
        self.xkj[:,1,:] = y_out
        self.xkj[:,2,:] = a_out

        self.Csi_j[:,0,:]  = self.ksi_L
        self.Csi_j[:,1,:]  = self.ksi_V
        self.Csi_j[:,2,:]  = self.ksi_A

        self.rho_j[:,0,:]  = self.rho_L
        self.rho_j[:,1,:]  = self.rho_V
        self.rho_j[:,2,:]  = self.rho_A

    def water_last(self, arr):
        arr_out = np.empty_like(arr)
        arr_out[-1,:] = arr[0,:]
        arr_out[0:-1,:] = arr[1:,:]
        return arr_out

    def equilibrium_ratio_aqueous(self, z):

        self.Kw = np.zeros_like(z)
        self.Kw[0] = 0.999 / z[0][0]
        self.Kw[1:] = 0.001 / (len(z) - 1) / z[1:,[0]]

        #self.Kw = 0.001 / (len(z) - 1) / z
        #self.Kw[3] = 0.999 / z[-1]

    def equilibrium_ratio_2flash(self, K_2flash):
        self.K = K_2flash.copy()

    def pressure_water_saturation(self):
        ' Correlation of Wagner and Saul 1987 '
        a1 = -7.85823
        a2 = 1.83991
        a3 = -11.7811
        a4 = 22.6705
        a5 = -15.9393
        a6 = 1.77516

        Tc = 647.14 # Kelvin
        Pc = 22.064e6 # Pascal
        tal = 1 - (self.T/Tc)

        ln_P_Pc = (Tc/self.T)*(a1*tal + a2*(tal**1.5) + a3*(tal**3) + a4*(tal**3.5) + a5*(tal**4) + a6*(tal**7.5))
        self.Pw_sat = np.exp(ln_P_Pc) * Pc

    """-------------Below starts phase stability test calculation -----------------"""

    def Stability_2phase(self, z, ponteiro_stab_check):
        ''' In the lnphi function: 0 stands for vapor phase and 1 for liquid '''

    #*****************************Test one**********************************#
        stationary_points = np.empty((5, len(ponteiro_stab_check[ponteiro_stab_check])))
        # Trial phase is liquid
        Y = np.empty(z.shape)
        lnphiz = np.empty(z.shape)
        ponteiro = np.copy(ponteiro_stab_check)
        Y[:,ponteiro] = z[:,ponteiro] / self.K[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphiz[:,ponteiro] = self.EOS.lnphi(z[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
        K_value1 = 1 / self.K[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = self.EOS.lnphi(y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        #lny = np.log(y)
        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphiz[:] - np.log(z)))
        #print(f'TPD do teste 1: {TPD}')
        #import pdb; pdb.set_trace()
        stationary_point1 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[0] = stationary_point1

    #*****************************Test two**********************************#
        # Trial phase is vapour
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = self.K[:,ponteiro] * z[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphiz[:,ponteiro] = self.EOS.lnphi(z[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
        K_value2 = self.K[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = self.EOS.lnphi(y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        #lny = np.log(y)
        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphiz[:] - np.log(z)))
        #print(f'TPD do teste 2: {TPD}')

        stationary_point2 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[1] = stationary_point2


    #*****************************Test three**********************************#
        # Trial phase is liquid
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = z[:,ponteiro] / ((self.K[:,ponteiro])**(1/3))
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphiz[:,ponteiro] = self.EOS.lnphi(z[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
        K_value3 = 1 / ((self.K[:,ponteiro])**(1/3))

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = self.EOS.lnphi(y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        #lny = np.log(y)
        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphiz[:] - np.log(z)))
        #print(f'TPD do teste 3: {TPD}')

        stationary_point3 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[2] = stationary_point3



    #*****************************Test four**********************************#
        # Trial phase is vapour
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = z[:,ponteiro] * ((self.K[:,ponteiro])**(1/3))
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphiz[:,ponteiro] = self.EOS.lnphi(z[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
        K_value4 = (self.K[:,ponteiro])**(1/3)

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = self.EOS.lnphi(y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        #lny = np.log(y)
        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphiz[:] - np.log(z)))
        #print(f'TPD do teste 4: {TPD}')

        stationary_point4 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[3] = stationary_point4


    #*****************************Test five**********************************#
        # Trial phase is aqueous
        ponteiro = np.copy(ponteiro_stab_check)
        Y[:,ponteiro] = z[:,ponteiro] * self.Kw[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]

        lnphiz[:,ponteiro] = self.EOS.lnphi( z[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro]) # both ph_L ?
        K_value5 = self.Kw[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = self.EOS.lnphi( y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro]) # both ph_L ?
            Y[:,ponteiro] = np.exp(np.log(z[:,ponteiro]) + lnphiz[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux


        #lny = np.log(y)
        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphiz[:] - np.log(z)))
        #print(f'TPD do teste 5: {TPD}')

        stationary_point5 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[4] = stationary_point5

        t1 = time.time()
        print('stability test time:', t1-t0)

        #return stationary_point1, stationary_point2, stationary_point3, stationary_point4, stationary_point5
        return stationary_points, [K_value1, K_value2, K_value3, K_value4, K_value5]




    """-------------Below starts 3 phase stability test calculation -----------------"""

    def Stability_3phase(self, x, ponteiro_stab_check):
        'Testing stability in system'
    #*****************************Test one**********************************#
        t0 = time.time()

        stationary_points = np.empty((5, len(ponteiro_stab_check[ponteiro_stab_check])))
        print(f'stationary_points: {stationary_points}')

        # Trial phase is liquid
        Y = np.empty(x.shape)
        lnphix = np.empty(x.shape)
        ponteiro = np.copy(ponteiro_stab_check)
        Y[:,ponteiro] = x[:,ponteiro] / self.Kwilson[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphix[:,ponteiro] = self.EOS.lnphi( x[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
        K_value1 = 1 / self.Kwilson[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = self.EOS.lnphi( y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(x[:,ponteiro]) + lnphix[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        #lny = np.log(y[:, ponteiro_stab_check])
        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphix[:] - np.log(x)))
        #print(f'TPD do teste 1: {TPD}')

        stationary_point1 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[0] = stationary_point1

    #*****************************Test two**********************************#
        # Trial phase is vapour
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = self.Kwilson[:,ponteiro] * x[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphix[:,ponteiro] = self.EOS.lnphi( x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
        K_value2 = self.Kwilson[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = self.EOS.lnphi( y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(x[:,ponteiro]) + lnphix[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux


        lny = np.log(y)

        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphix[:] - np.log(x)))
        #print(f'TPD do teste 2: {TPD}')

        stationary_point2 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[1] = stationary_point2


    #*****************************Test three**********************************#
        # Trial phase is liquid
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = x[:,ponteiro] / ((self.Kwilson[:,ponteiro])**(1/3))
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphix[:,ponteiro] = self.EOS.lnphi( x[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
        K_value3 = 1 / ((self.Kwilson[:,ponteiro])**(1/3))

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = self.EOS.lnphi( y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(x[:,ponteiro]) + lnphix[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        lny = np.log(y)

        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphix[:] - np.log(x)))
        #print(f'TPD do teste 3: {TPD}')

        stationary_point3 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[2] = stationary_point3



    #*****************************Test four**********************************#
        # Trial phase is vapour
        ponteiro = np.copy(ponteiro_stab_check)

        Y[:,ponteiro] = x[:,ponteiro] * ((self.Kwilson[:,ponteiro])**(1/3))
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]
        lnphix[:,ponteiro] = self.EOS.lnphi( x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
        K_value4 = (self.Kwilson[:,ponteiro])**(1/3)

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = self.EOS.lnphi( y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            Y[:,ponteiro] = np.exp(np.log(x[:,ponteiro]) + lnphix[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

        lny = np.log(y)

        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphix[:] - np.log(x)))
        #print(f'TPD do teste 4: {TPD}')

        stationary_point4 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[3] = stationary_point4


    #*****************************Test five**********************************#
        # Trial phase is aqueous
        ponteiro = np.copy(ponteiro_stab_check)
        Y[:,ponteiro] = x[:,ponteiro] * self.Kw[:,ponteiro]
        y = Y / np.sum(Y, axis = 0)[np.newaxis,:]

        lnphix[:,ponteiro] = self.EOS.lnphi( x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro]) # both ph_L ?
        K_value5 = self.Kw[:,ponteiro]

        while any(ponteiro):
            Y_old = np.copy(Y[:,ponteiro])
            lnphiy = self.EOS.lnphi( y[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro]) # both ph_L ?
            #import pdb; pdb.set_trace()
            Y[:,ponteiro] = np.exp(np.log(x[:,ponteiro]) + lnphix[:,ponteiro] - lnphiy)
            y[:,ponteiro] = Y[:,ponteiro] / np.sum(Y[:,ponteiro], axis = 0)[np.newaxis,:]
            stop_criteria = np.max(abs(Y[:,ponteiro] / Y_old - 1), axis = 0)
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux
            #import pdb; pdb.set_trace()

        #lny = np.log(y)
        #TPD = np.sum(y[:]*(lnphiy[:] + lny[:] - lnphix[:] - np.log(x)))
        #print(f'TPD do teste 5: {TPD}')

        stationary_point5 = np.sum(Y[:,ponteiro_stab_check], axis = 0)
        stationary_points[4] = stationary_point5

        #print(f'pontos estacionarios: {stationary_points}')

        t1 = time.time()
        print('stability time:', t1-t0)

        #return stationary_point1, stationary_point2, stationary_point3, stationary_point4, stationary_point5
        return stationary_points, [K_value1, K_value2, K_value3, K_value4, K_value5]





    """-------------Below starts biphasic flash calculations-----------------"""

    def molar_properties_3phase(self, z, ponteiro):
        ksi_L, ksi_V, ksi_A, rho_L, rho_V, rho_A = self.molar_properties_Whitson_3phase(z, ponteiro)
        #ponteiro = self.molar_properties_Lapene_3phase( z, ponteiro)
        return ksi_L, ksi_V, ksi_A, rho_L, rho_V, rho_A


    ' Try of triphasic flash '
    def molar_properties_Whitson_3phase(self, z, ponteiro):

        self.K_A = self.K.copy()
        self.K_V = self.Kwilson.copy()

        self.V[ponteiro] = 0.5
        self.A[ponteiro] = 0.5

        #self.K_A = self.Kw[:,ponteiro]
        self.K_A[:, ponteiro] = self.Kw[:, ponteiro]

        #t0 = time.time()

        razao_vl = np.ones(z.shape)/2
        razao_al = np.ones(z.shape)/2
        ponteiro_save = np.copy(ponteiro)
        #import pdb; pdb.set_trace()

        Zl = np.empty_like(self.V)
        Zv = np.empty_like(Zl)
        Za = np.empty_like(Zl)

        while any(ponteiro):
            self.solve_objective_function_Whitson_for_V_and_A(z, self.V, self.A, np.copy(ponteiro))

            lnphil, Zl[ponteiro] = self.lnphi_Z_based_on_deltaG(self.x[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])
            lnphiv, Zv[ponteiro] = self.lnphi_Z_based_on_deltaG(self.y[:,ponteiro], self.P[ponteiro], self.ph_V[ponteiro])
            lnphia, Za[ponteiro] = self.lnphi_Z_based_on_deltaG(self.a[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])

            fv = np.exp(lnphiv) * (self.y[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            fl = np.exp(lnphil) * (self.x[:,ponteiro] * self.P[ponteiro][np.newaxis,:])
            fa = np.exp(lnphia) * (self.a[:,ponteiro] * self.P[ponteiro][np.newaxis,:])

            fv[(abs(fl)<1e-300) + (abs(fv)<1e-300)] = fl[(abs(fl)<1e-300) + (abs(fv)<1e-300)]
            fv[fv == 0] = 1e-30
            razao_vl[:,ponteiro] = fl/fv
            fa[(abs(fl)<1e-300) + (abs(fa)<1e-300)] = fl[(abs(fl)<1e-300) + (abs(fa)<1e-300)]
            fa[fa == 0] = 1e-30
            razao_al[:,ponteiro] = fl/fa

            self.K_V[:,ponteiro] *= razao_vl[:,ponteiro]
            self.K_A[:,ponteiro] *= razao_al[:,ponteiro]

            stop_criteria_vl = np.max(abs((fv/(fl + 1e-15)) - 1), axis = 0)
            stop_criteria_al = np.max(abs((fa/(fl + 1e-15)) - 1), axis = 0)
            stop_criteria = stop_criteria_vl.copy()
            stop_criteria[stop_criteria_al > stop_criteria_vl] = stop_criteria_al[stop_criteria_al > stop_criteria_vl]
            #stop_criteria = max(stop_criteria_vl, stop_criteria_al)

            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux
            #print(ponteiro)
            #import pdb; pdb.set_trace()
            #ponteiro[((self.V)<0) + ((self.V)>1)] = False


        #t1 = time.time()
        #print('Whitson-Michelsen time for 3phase flash:', t1-t0)

        ksi_L, ksi_V, ksi_A, rho_L, rho_V, rho_A = \
        self.update_EOS_dependent_properties_3ph(Zl, Zv, Za)

        return ksi_L, ksi_V, ksi_A, rho_L, rho_V, rho_A

    def solve_objective_function_Whitson_for_V_and_A(self, z, V, A, ponteiro):

        ponteiro_save = np.copy(ponteiro)
        i = 0
        V_and_A = np.array([V, A])
        #import pdb; pdb.set_trace()
        #self.K_A = self.Kw
        while any(ponteiro):
            Vold = V[ponteiro].copy()
            Aold = A[ponteiro].copy()
            #V_and_A = np.array([Vold, Aold])

            f = np.zeros([2, Vold.size])

            f[0] = np.sum((self.K_V[:,ponteiro] - 1) * z[:,ponteiro] / (1 + \
                V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - 1) + \
                A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)), axis = 0)

            f[1] = np.sum((self.K_A[:,ponteiro] - 1) * z[:,ponteiro] / (1 + \
                    V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - 1) + \
                    A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)), axis = 0)

            df_VdV = - np.sum(((self.K_V[:,ponteiro] - 1) ** 2) * z[:,ponteiro] / ((1 + \
                    V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - 1) + \
                    A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)) ** 2), axis = 0)

            df_VdA = - np.sum((self.K_V[:,ponteiro] - 1) * (self.K_A[:,ponteiro] - 1) *\
                                    z[:,ponteiro] / ((1 + V[ponteiro][np.newaxis,:] * \
                                    (self.K_V[:,ponteiro] - 1) + A[ponteiro][np.newaxis,:] * \
                                    (self.K_A[:,ponteiro] - 1)) ** 2), axis = 0)

            df_AdV = - np.sum((self.K_V[:,ponteiro] - 1) * (self.K_A[:,ponteiro] - 1) \
                                 * z[:,ponteiro] / ((1 + V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - \
                                   1) + A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)) ** 2), axis = 0)

            df_AdA = - np.sum(((self.K_A[:,ponteiro] - 1) ** 2) * z[:,ponteiro] / ((1 + \
                    V[ponteiro][np.newaxis,:] * (self.K_V[:,ponteiro] - 1) + \
                    A[ponteiro][np.newaxis,:] * (self.K_A[:,ponteiro] - 1)) ** 2), axis = 0)


            Jacobian = np.empty([ponteiro[ponteiro].size, 2, 2])
            Jacobian[:,0,0] = df_VdV
            Jacobian[:,0,1] = df_VdA
            Jacobian[:,1,0] = df_AdV
            Jacobian[:,1,1] = df_AdA

            try:
                Jacobian_inv = np.linalg.inv(Jacobian)
            except: import pdb; pdb.set_trace()

            '''
            solucao = np.zeros_like(f.T)
            for i in range(ponteiro[ponteiro].size):
                solucao[i] = f.T[i].dot(Jacobian_inv[i])
            '''
            solucao_aux = f.T.dot(Jacobian_inv)
            solucao = np.diagonal(solucao_aux).T

            V_and_A[:, ponteiro] = V_and_A[:, ponteiro] - solucao.T
            #import pdb; pdb.set_trace()

            #print(V_and_A)
            V[ponteiro] = V_and_A[0, ponteiro].copy()
            A[ponteiro] = V_and_A[1, ponteiro].copy()

            V_aux = V[ponteiro]
            V_aux[V_aux > 1.0] = 1.0
            V_aux[V_aux < 0.0] = 0.0
            #V_aux[V_aux > 1.0] = 0.5 * (1.0 + Vold[V_aux > 1.0])
            #V_aux[V_aux < 0.0] = 0.5 * (0.0 + Vold[V_aux < 0.0])
            V[ponteiro] = V_aux

            A_aux = A[ponteiro]
            A_aux[A_aux > 1.0] = 1.0
            A_aux[A_aux < 0.0] = 0.0
            #A_aux[A_aux > 1.0] = 0.5 * (1.0 + Aold[A_aux > 1.0])
            #A_aux[A_aux < 0.0] = 0.5 * (0.0 + Aold[A_aux < 0.0])
            A[ponteiro] = A_aux

            #import pdb; pdb.set_trace()
            stop_criteria_V = abs(V[ponteiro] - Vold)
            stop_criteria_A = abs(A[ponteiro] - Aold)
            stop_criteria = stop_criteria_V.copy()
            stop_criteria[stop_criteria_A > stop_criteria_V] = stop_criteria_A[stop_criteria_A > stop_criteria_V]
            #stop_criteria = max(abs(V[ponteiro] - Vold), abs(A[ponteiro] - Aold))

            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

            i+=1
            if i>=3000:
                print('maxit in triphasic flash')
                import pdb; pdb.set_trace()

        self.V[ponteiro_save] = V[ponteiro_save]
        self.A[ponteiro_save] = A[ponteiro_save]
        if any((self.V[ponteiro_save] + self.A[ponteiro_save])> 1.0):
            import pdb; pdb.set_trace()
        self.L[ponteiro_save] = (1 - self.V[ponteiro_save] - self.A[ponteiro_save])

        self.x[:,ponteiro_save] = z[:,ponteiro_save] / (1 + self.V[ponteiro_save][np.newaxis,:] * (self.K_V[:,ponteiro_save] - 1) + \
                                    self.A[ponteiro_save][np.newaxis,:] * (self.K_A[:,ponteiro_save] - 1))
        self.y[:,ponteiro_save] = self.K_V[:,ponteiro_save] * self.x[:,ponteiro_save]
        self.a[:,ponteiro_save] = self.K_A[:,ponteiro_save] * self.x[:,ponteiro_save]

        return True

    def molar_properties_Lapene_3phase(self, z, ponteiro):

        V = self.V[ponteiro]
        x = self.x[:, ponteiro]
        y = self.y[:, ponteiro]

        self.K_V[0] = (ctes.Pc[0]/self.P[ponteiro])*(self.T/ctes.Tc[0]) # Adicionar ponteiro no T para o termico

        y[0] = self.Pw_sat / self.P[ponteiro]
        self.K2w = y[0].copy()
        x[0] = y[0] / self.K_V[0]
        self.a[:, ponteiro] = 0
        self.a[0, ponteiro] = 1

        #t0 = time.time()

        razao_vl = np.ones(x.shape)/2
        razao_av = np.ones(self.P[ponteiro].shape)/2
        ponteiro_save = np.copy(ponteiro)

        Zl = np.empty_like(self.V)
        Zv = np.empty_like(Zl)
        Za = np.empty_like(Zl)

        while any(ponteiro):
            y[0] = self.K2w.copy()
            x[0] = y[0] / self.K_V[0]

            K_V_max = np.amax(self.K_V[1:,:], axis=0)
            K_V_min = np.amin(self.K_V[1:,:], axis=0)
            Kw_ast = (1 - y[0]) / (1 - x[0])
            try:
                Kwz = (1 - z[0, ponteiro_save]) / (1 - x[0])
                #Kwz = (1 - z[0, ponteiro]) / (1 - x[0])
            except: import pdb; pdb.set_trace()

            Vmax = Kwz / (Kw_ast - K_V_min)
            Vmin = Kwz / (Kw_ast - K_V_max)

            try:
                ponteiro[ponteiro_save] = ~((K_V_max < Kw_ast) + (K_V_min > Kw_ast))
            except: import pdb; pdb.set_trace()
            '''
            if (K_V_max < Kw_ast or K_V_min > Kw_ast):
                print('Solve for V a classical RR equation')
                # colocar False no ponteiro
                import pdb; pdb.set_trace()
            '''

            V_dash = (z[0, ponteiro] - x[0]) / (y[0] - x[0])
            Objetive_function_V_dash = np.sum((self.K_V[1:] - Kw_ast) * z[1:,ponteiro] / (Kwz + \
                                        V_dash * (self.K_V[1:] - Kw_ast)), axis = 0)


            ' determine position of V with respect V_dash '
            V_bigger_then_Vdash = np.ones(len(V), dtype = bool)
            ponteiro_aux = np.ones(len(V), dtype = bool)

            try:
                ponteiro_aux[ponteiro_aux] = (V_dash > Vmin) & (V_dash < Vmax)
            except: import pdb; pdb.set_trace()
            V_bigger_then_Vdash[ponteiro_aux] = Objetive_function_V_dash[ponteiro_aux] > 0
            '''
            if (V_dash > Vmin and V_dash < Vmax):
                if (Objetive_function_V_dash > 0):
                    V_bigger_then_Vdash = True
                else:
                    V_bigger_then_Vdash = False
            '''

            #ponteiro_aux[~ponteiro_aux] = (V_dash < Vmin)[~ponteiro_aux]
            ponteiro_aux = V_dash < Vmin
            V_bigger_then_Vdash[ponteiro_aux] = True
            '''
            elif (V_dash < Vmin):
                V_bigger_then_Vdash = True
            '''

            #ponteiro_aux[~ponteiro_aux] = V_dash > Vmax
            ponteiro_aux = V_dash > Vmax
            V_bigger_then_Vdash[ponteiro_aux] = False
            '''
            elif (V_dash > Vmax):
                V_bigger_then_Vdash = False
            '''

            ' determine presence or absence of water '
            cond_ext = ~V_bigger_then_Vdash
            cond_int = y[0] < x[0]

            cond1 = cond_ext & cond_int
            cond2 = ~cond_ext & ~cond_int

            aux = self.A[ponteiro]
            aux[cond1 + cond2] = 0.0
            self.A[ponteiro] = aux

            aux2 = ponteiro[ponteiro]
            aux2[cond1 + cond2] = False
            ponteiro[ponteiro] = aux2

            '''
            if not V_bigger_then_Vdash:
                if (y[0] < x[0]):
                    self.A = 0
                    print('Não tem água no sistema')
                    # ponteiro falso
                else:
                    print('Tem água no sistema')
            if V_bigger_then_Vdash:
                if (y[0] > x[0]):
                    self.A = 0
                    print('Não tem água no sistema')
                    # ponteiro falso
                else:
                    print('Tem água no sistema')
            '''

            if any(ponteiro) == False:
                return ponteiro_save

            V, x, y = self.solve_objective_function_Lapene_for_V(z, x, y, V, Kw_ast, Kwz, np.copy(ponteiro))

            lnphil, Zl[ponteiro] = self.lnphi_Z_based_on_deltaG(x, self.P[ponteiro], self.ph_L[ponteiro])
            lnphiv, Zv[ponteiro] = self.lnphi_Z_based_on_deltaG(y, self.P[ponteiro], self.ph_V[ponteiro])
            lnphia, Za[ponteiro] = self.lnphi_Z_based_on_deltaG(self.a[:,ponteiro], self.P[ponteiro], self.ph_L[ponteiro])


            fv = np.exp(lnphiv) * (y * self.P[ponteiro][np.newaxis,:])
            fl = np.exp(lnphil) * (x * self.P[ponteiro][np.newaxis,:])
            fa = np.exp(lnphia) * (self.a[:,ponteiro] * self.P[ponteiro][np.newaxis,:])

            fl[fl==0] = 1e-30
            razao_vl = fl/fv
            fv[fv==0] = 1e-30
            razao_av = fa[0] / fv[0] * (1 + 1e-10) #?????
            #razao_vl[:,ponteiro] = np.divide(self.fl, self.fv, out = razao_vl[:,ponteiro] / razao_vl[:,ponteiro] * (1 + 1e-10),
                              #where = self.fv != 0)
            #razao_av[ponteiro] = self.fa[0,ponteiro] / self.fv[0,ponteiro] * (1 + 1e-10)

            self.K_V *= razao_vl
            self.K2w *= razao_av
            #self.K_V[ponteiro] = razao_vl[ponteiro] * self.K_V[ponteiro]
            #self.K2w[ponteiro] = self.K2w[ponteiro] * razao_av

            #stop_criteria = np.max(abs((self.fl/(self.fv + 1e-15)) - 1), axis = 0)
            stop_criteria = (np.sum(((fl/fv) - 1)**2, axis=0))**0.5
            ponteiro_aux = ponteiro[ponteiro]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro[ponteiro] = ponteiro_aux

            #ponteiro[((self.V)<0) + ((self.V)>1)] = False

        #t1 = time.time()
        #print('Lapene time for 3phase flash:', t1-t0)

        V_ast = (1 - z[0, ponteiro_save]) / (1 - y[0])
        V[V>V_ast] = V_ast[V>V_ast]
        #if V > V_ast:
            #V = V_ast
            #print('V maior que o V asterisco')
        V[V<0] = 0.0
        #if V < 0:
            #V = 0

        # Aqui atualizar tudo
        #self.V = V[ponteiro_save].copy()
        #self.x[:,ponteiro_save] = x.copy()
        #self.y[:,ponteiro_save] = y.copy()
        self.V[ponteiro_save] = V
        self.x[:, ponteiro_save] = x
        self.y[:, ponteiro_save] = y
        self.A[ponteiro_save] = (z[0, ponteiro_save] + self.V[ponteiro_save]*(self.x[0, ponteiro_save] - self.y[0, ponteiro_save]) - self.x[0, ponteiro_save]) / (1 - self.x[0, ponteiro_save])
        self.L[ponteiro_save] = 1 - self.A[ponteiro_save] - self.V[ponteiro_save]
        self.L[self.L<0] = 0.0


        ksi_L, ksi_V, ksi_A, rho_L, rho_V, rho_A = self.update_EOS_dependent_properties_3ph(Zl, Zv, Za)

        return ksi_L, ksi_V, ksi_A, rho_L, rho_V, rho_A

    def solve_objective_function_Lapene_for_V(self, z, x, y, V, Kw_ast, Kwz, ponteiro):
        ' Leibovici–Neoschil window '

        '''
        # Alternativa para calcular o Vmin
        aux = self.K_V > Kw_ast
        aux_vols = aux.astype(int).sum(axis=0)
        vols = len(aux_vols[aux_vols>0])
        Vmin = np.empty(vols)
        z_aux = z[:, ponteiro]

        for i in range(self.Nc):
            Vmin[aux_vols[aux_vols>0]==(self.Nc-i)] = np.max((self.K_V[:,aux_vols==self.Nc-i] * z_aux[:,aux_vols==self.Nc-i] - \
            Kwz[aux_vols==self.Nc-i])/(self.K_V[:,aux_vols==self.Nc-i] - Kw_ast[aux_vols==self.Nc-i]), axis=0, initial=0)
        '''

        z_aux = z[:, ponteiro]
        equation1 = (self.K_V * z_aux - Kwz)/(self.K_V - Kw_ast)
        try:
            equation1[(self.K_V <= Kw_ast)] = -np.inf
        except: import pdb; pdb.set_trace()

        Vmin = np.max(equation1, axis=0)

        equation2 = (Kwz - z_aux)/(Kw_ast - self.K_V)
        try:
            equation2[(self.K_V >= Kw_ast)] = np.inf
        except: import pdb; pdb.set_trace()

        Vmax = np.min(equation2, axis=0)

        #V[ponteiro] = 0.5*(Vmax + Vmin) # Chute inicial
        V = 0.5*(Vmax + Vmin) # Chute inicial

        ponteiro_save = np.copy(ponteiro)
        i = 0

        while any(ponteiro_save):
            #Vold = V[ponteiro]
            Vold = V

            G = np.sum((self.K_V[1:] - Kw_ast) * z[1:,ponteiro] / (Kwz + \
                    V * (self.K_V[1:] - Kw_ast)), axis = 0)

            dG = - np.sum(((self.K_V[1:] - Kw_ast)**2) * z[1:,ponteiro] / ((Kwz + \
                    V * (self.K_V[1:] - Kw_ast))**2), axis = 0)

            #G = np.sum((self.K_V[1:] - Kw_ast) * z[1:,ponteiro_save] / (Kwz + \
                    #V * (self.K_V[1:] - Kw_ast)), axis = 0)

            #dG = - np.sum(((self.K_V[1:] - Kw_ast)**2) * z[1:,ponteiro_save] / ((Kwz + \
                    #V * (self.K_V[1:] - Kw_ast))**2), axis = 0)

            #V[ponteiro] = V[ponteiro] - G / dG # Newton-Raphson iterative method
            V = V - G / dG # Newton-Raphson iterative method

            #V_aux = V[ponteiro]
            V_aux = V
            #V_aux[V_aux > Vmax[ponteiro]] = 0.5 * (Vmax[ponteiro] + Vold)[V_aux > Vmax[ponteiro]] #(Vmax + Vold)/2
            #V_aux[V_aux < Vmin[ponteiro]] = 0.5 * (Vmin[ponteiro] + Vold)[V_aux < Vmin[ponteiro]]#(Vmin + Vold)/2
            V_aux[V_aux > Vmax] = 0.5 * (Vmax + Vold)[V_aux > Vmax] #(Vmax + Vold)/2
            V_aux[V_aux < Vmin] = 0.5 * (Vmin + Vold)[V_aux < Vmin] #(Vmin + Vold)/2
            #V[ponteiro] = V_aux
            V = V_aux

            #stop_criteria = abs((V[ponteiro] / Vold) - 1)
            #import pdb; pdb.set_trace()
            #stop_criteria = np.empty(ponteiro_save[ponteiro_save].shape)
            #stop_criteria[ponteiro_save[ponteiro_save]] = abs((V / Vold) - 1)
            stop_criteria = abs((V / Vold) - 1)
            ponteiro_aux = ponteiro_save[ponteiro]
            #ponteiro_aux = ponteiro_save[ponteiro_save]
            ponteiro_aux[stop_criteria < 1e-9] = False
            ponteiro_save[ponteiro] = ponteiro_aux
            #ponteiro_save[ponteiro_save] = ponteiro_aux
            i+=1
            if i>=3000:
                print('maxit in triphasic flash')
                import pdb; pdb.set_trace()


        x[1:] = z[1:, ponteiro] / (1 + V*(self.K_V[1:] - 1 + (y[0] - x[0])/(1 - x[0])) + (x[0] - z[0, ponteiro])/(1 - x[0]))
        y[1:] = x[1:] * self.K_V[1:]
        y[0] = self.K2w.copy()
        x[0] = y[0] / self.K_V[0]

        #self.V[ponteiro_save] = V[ponteiro_save]
        #self.x[:,ponteiro_save] = x
        #self.y[:,ponteiro_save] = y
        return V, x, y


    def update_EOS_dependent_properties_3ph(self, Zl, Zv, Za):
        #self.EOS = ctes.EOS_class(self.P, self.T)
        ksi_L, rho_L = self.get_EOS_dependent_properties(self.T, self.x, self.P, Zl)
        ksi_V, rho_V = self.get_EOS_dependent_properties(self.T, self.y, self.P, Zv)
        ksi_A, rho_A = self.get_EOS_dependent_properties(self.T, self.a, self.P, Za)

        return ksi_L, ksi_V, ksi_A, rho_L, rho_V, rho_A
