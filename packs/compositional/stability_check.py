"""Check stability of a system and find its composition at a thermodynamic equilibrium."""
import numpy as np
from ..directories import data_loaded
from scipy.misc import derivative
from . import equation_of_state
import math
# import matplotlib.pyplot as plt
## Encontrar os pontos estacionarios. Estes correspondem aos pontos nos quais a derivada de g com respeito a Y é 0
## Todas as equações foram confirmadas pelo livro do Dandekar e a biblioteca do thermo
# ph=0 - vapor. ph=1 - liquid.

class StabilityCheck:
    """Check for stability of a thermodynamic equilibrium and returns the
    equilibrium phase compositions (perform the flash calculation)."""

    def __init__(self, P, T, kprop):
        self.T = T
        self.P = P
        EOS_class = getattr(equation_of_state, data_loaded['compositional_data']['equation_of_state'])
        self.EOS = EOS_class(P, T, kprop)

    def run(self, z, kprop):
        self.equilibrium_ratio_Wilson(kprop)

        if any(z <= 0):
            self.molar_properties(kprop, z)
        else:
            sp1,sp2 = self.StabilityTest(kprop, z)
            if np.round(sp1,8) > 1 or np.round(sp2,8) > 1: self.molar_properties(kprop, z)
            else: #tiny manipulation
                self.x = z; self.y = z
                self.bubble_point_pressure()
                if self.P > self.Pbubble: self.L = 1; self.V = 0; ph = 1
                else: self.L = 0; self.V = 1; ph = 0

    def equilibrium_ratio_Wilson(self, kprop):
        self.K = np.exp(5.37 * (1 + kprop.w) * (1 - kprop.Tc / self.T)) * \
                (kprop.Pc / self.P)


    """------------------- Stability test calculation -----------------------"""

    def StabilityTest(self, kprop, z):
        ''' In the lnphi function: 0 stands for vapor phase and 1 for liquid '''
    #****************************INITIAL GUESS******************************#
    ## Both approaches bellow should be used in case the phase is in the critical region

    #*****************************Test one**********************************#
        #Used alone when the phase investigated (z) is clearly vapor like (ph = 0)

        Y = z / self.K
        Yold = 0.9 * Y
        y = Y / sum(Y)
        lnphiz = self.EOS.lnphi(kprop, z, 0)
        while max(abs(Y / Yold - 1)) > 1e-9: #convergência
            Yold = np.copy(Y)
            lnphiy = self.EOS.lnphi(kprop, y, 1)
            Y = np.exp(np.log(z) + lnphiz - lnphiy)
            y = Y / sum(Y)
        stationary_point1 = sum(Y)

    #*****************************Test two**********************************#
        #Used alone when the phase investigated (z) is clearly liquid like (ph == 1)

        Y = self.K * z
        Y_old = 0.9 * Y
        y = Y / sum(Y)
        lnphiz = self.EOS.lnphi(kprop, z, 1)
        while max(abs(Y / Y_old - 1)) > 1e-9:
            Y_old = np.copy(Y)
            lnphiy = self.EOS.lnphi(kprop, y, 0)
            Y = np.exp(np.log(z) + lnphiz - lnphiy)
            y = Y / sum(Y)
        stationary_point2 = sum(Y)

        return stationary_point1,stationary_point2


    """-------------------- Biphasic flash calculations ---------------------"""

    def molar_properties(self, kprop, z):
        self.fv = 2 * np.ones(len(z)); self.fl = np.ones(len(z)) #entrar no primeiro loop
        if kprop.Nc<= 2: self.molar_properties_Whitson(kprop, z)
        else: self.molar_properties_Yinghui(kprop, z)

    def deltaG_molar(self, kprop, l, ph):
        lnphi = [self.EOS.lnphi(kprop, l, 1 - ph), self.EOS.lnphi(kprop, l, ph)]
        deltaG_molar = sum(l * (lnphi[1 - ph] - lnphi[ph]))
        if deltaG_molar < 0: ph = 1 - ph
        return ph

    def deltaG_molar_vectorized(EOS, kprop, l, ph):
        lnphi = np.empty([2, len(ph)])
        lnphi[0,:] = EOS.lnphi(kprop, l, 1 - ph)
        lnphi[1,:] = EOS.lnphi(kprop, l, ph)

        deltaG_molar = np.sum(l * (lnphi[1 - ph] - lnphi[ph]), axis=0)
        ph[deltaG_molar<0] = 1 - ph[deltaG_molar<0]
        return ph

    def lnphi_based_on_deltaG(self,kprop, l, ph):
        ph = np.array(ph)[:,np.newaxis]
        ph = StabilityCheck.deltaG_molar_vectorized(self.EOS, kprop, l, ph)
        return self.EOS.lnphi(kprop, l, ph)

    def solve_objective_function_Yinghui(self, z1, zi, z, K1, KNc, Ki):

        x1_min = z1 * (1 - KNc) / (K1 - KNc)
        x1_max = (1 - KNc) / (K1 - KNc)

        if any(z < 0):
            theta = np.ones(len(z))
            theta[self.K > 1] = (1 - KNc) / (self.K[self.K > 1] - KNc)
            aux_eq = (self.K - 1) * z1 / (z * (K1 - 1) / theta - (K1 - self.K))
            if all((self.K[z != 0] - 1) * z1 / z[z != 0] > 0):
                aux_eq = aux_eq[aux_eq >= 0] #se arr<0 é descartado
                x1_max = min(aux_eq)
            else: x1_min = max(np.append(arr, 0))

        if x1_min > x1_max: raise ValueError('There is no physical root')

        x1 = (x1_min + x1_max) / 2
        f = 1
        while abs(f) > 1e-10:
            f = 1 + ((K1 - KNc) / (KNc - 1)) * x1 + sum(((Ki - KNc) / (KNc - 1))
               * zi * (K1 - 1) * x1 / ((Ki - 1) * z1 + (K1 - Ki) * x1))
            df = ((K1 - KNc) / (KNc - 1)) + sum(((Ki - KNc) / (KNc - 1)) * zi *
                z1 * (K1 - 1)* (Ki - 1) / ((Ki - 1) * z1 + (K1 - Ki) * x1) ** 2)
            x1 = x1 - f/df #Newton-Raphson iterative method
            if (x1 > x1_max) or (x1 < x1_min): x1 = (x1_min + x1_max)/2 #recommended
            if f * df > 0: x1_max = x1
            if f * df < 0: x1_min = x1

        xi = (K1 - 1) * zi * x1 / ((Ki - 1) * z1 + (K1 - Ki) * x1)
        self.x[self.K == K1] = x1
        self.x[self.K == KNc] = 1 - sum(xi) - x1
        return xi

    def Yinghui_method(self, z):

        """ Shaping K to Nc-2 components by removing K1 and KNc """
        K1 = max(self.K); KNc = min(self.K)
        Ki = self.K[(self.K != K1) & (self.K != KNc)]

        """ Shaping z to Nc-2 components by removing z1 and zNc """
        z1 = z[self.K == K1]
        index1 = np.argwhere(self.K == K1)
        zi = np.delete(z, index1)
        indexNc = np.argwhere(self.K == KNc)
        zi = np.delete(zi, indexNc)

        #starting x
        # self.x = np.zeros(kprop.Nc)

        """ Solution """
        if z1 != 0:
            xi = self.solve_objective_function_Yinghui(z1, zi, z, K1, KNc, Ki)
        else:
            '''Explicit Calculation of xi'''
            xi = (K1 - 1) * zi / (K1 - Ki)
            self.x[self.K == KNc] = (K1 - 1) * z[self.K == KNc] / (K1 -
                                    self.K[self.K == KNc])
            self.x[self.K == K1] = 1 - sum(xi) - sum(self.x)

        #ainda não sei como tirar esse for
        for j in range(0, len(xi)):
            self.x[self.K == Ki[j]] = xi[j]

        self.y = self.K * self.x

    def molar_properties_Yinghui(self, kprop, z):
        #razao = fl/fv -> an arbitrary vector to enter in the iterative mode
        razao = np.ones(kprop.Nc)/2
        while max(abs(razao - 1)) > 1e-9:
            self.Yinghui_method(z)
            lnphil = self.lnphi_based_on_deltaG(kprop, self.x, 1)
            lnphiv = self.lnphi_based_on_deltaG(kprop, self.y, 0)
            self.fl = np.exp(lnphil) * (self.x * self.P)
            self.fv = np.exp(lnphiv) * (self.y * self.P)
            razao = np.divide(self.fl, self.fv, out = razao / razao * (1 + 1e-10),
                              where = self.fv != 0)
            self.K = razao * self.K
        self.V = (z[self.x != 0] - self.x[self.x != 0]) / (self.y[self.x != 0] -
                self.x[self.x != 0])
        self.V = self.V[0]

    def solve_objective_function_Whitson(self, z):
        """ Solving for V """
        Vmax = 1 / (1 - min(self.K))
        Vmin = 1 / (1 - max(self.K))
        #Vmin = ((K1-KNc)*z[self.K==K1]-(1-KNc))/((1-KNc)*(K1-1))
        #proposed by Li et al for Whitson method
        self.V = (Vmin + Vmax) / 2
        Vold = self.V / 2 #just to get into the loop

        while abs(self.V / Vold - 1) > 1e-8:
            Vold = self.V
            f = sum((self.K - 1) * z / (1 + self.V * (self.K - 1)))
            df = -sum((self.K - 1) ** 2 * z / (1 + self.V * (self.K - 1)) ** 2)
            self.V = self.V - f / df #Newton-Raphson iterative method

            if self.V > Vmax: self.V = Vmax #(Vmax + Vold)/2
            elif self.V < Vmin: self.V = Vmin #(Vmin + Vold)/2

        self.x = z / (1 + self.V * (self.K - 1))
        self.y = self.K * self.x

    def molar_properties_Whitson(self, kprop, z):
        razao = np.ones(kprop.Nc)/2
        while max(abs(self.fv / self.fl - 1)) > 1e-9:
            self.solve_objective_function_Whitson(z)
            lnphil = self.lnphi_based_on_deltaG(kprop, self.x, 1)
            lnphiv = self.lnphi_based_on_deltaG(kprop, self.y, 0)
            self.fv = np.exp(lnphiv) * (self.y * self.P)
            self.fl = np.exp(lnphil) * (self.x * self.P)
            razao = np.divide(self.fl, self.fv, out = razao / razao * (1 + 1e-10),
                              where = self.fv != 0)
            self.K = razao * self.K

    def bubble_point_pressure(self):
        #Isso vem de uma junção da Lei de Dalton com a Lei de Raoult
        Pv = self.K * self.P
        self.Pbubble = sum(self.x * Pv)

    '''def TPD(self, z): #ainda não sei onde usar isso
        x = np.zeros(self.Nc)

        #**********************Tangent Plane distance plot*********************#
        t = np.linspace(0.01, 0.99, 0.9 / 0.002) #vetor auxiliar
        TPD = np.zeros(len(t)) ##F

        for i in range(0, len(t)):
            aux = 0;
            lnphiz = self.lnphi(z, 1) #original phase

            #x = np.array([1-t[i],t[i]]) #new phase composition (1-t e t) - apenas válido para Nc=2 acredito eu.
            for k in range(0, kprop.Nc- 1):
                x[k] = (1 - t[i]) / (kprop.Nc- 1)
                x[kprop.Nc- 1] = t[i]

            ''''''O modo que x varia implica no formato de TPD. No presente exemplo,
            a fração molar do segundo componente de x varia direto com t, que é a
            variável de plotagem. Logo, a distancia dos planos tangentes será
            zero em z[Nc-1]. O contrário ocorreria''''''
            lnphix = self.lnphi(x, 0); #new phase (vapor- ph=2)
            for j in range(0,self.Nc):
                fix = math.exp(lnphix[j]) * x[j] * self.P
                fiz = math.exp(lnphiz[j]) * z[j] * self.P
                aux = aux + x[j] * ctes.R* self.T * (math.log(fix / fiz))
                TPD[i] = aux

        plt.figure(0)
        plt.plot(t, TPD)
        plt.xlabel('x')
        plt.ylabel('TPD')
        plt.show()
        return TPD'''
