"""Check stability of a thermodynamic equilibrium."""
import numpy as np
from scipy.misc import derivative
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
        # self.run(z, kprop)
        #StabilityCheck.TPD(self)

    def run(self, z, kprop):
        #if all(z==0):
        #    self.x = np.zeros(kprop.Nc); self.y = np.zeros(kprop.Nc); self.L = 0.; self.V = 0.; self.K = 1.;
        #    self.Mw_L = 0; self.Mw_V = 0.; self.rho_L = 0.; self.rho_V = 0; self.ksi_L = 0; self.ksi_V = 0
        #    self.z = np.zeros(kprop.Nc)
        self.R = kprop.R
        self.equilibrium_ratio_Wilson(kprop)

        #if sum(z < 0):
        #    self.x = np.zeros(kprop.Nc); self.y = np.zeros(kprop.Nc); self.L = -1; self.V = -1
        #else:
        if any(z <= 0):
            self.molar_properties(kprop, z)
        else:
            sp1,sp2 = self.StabilityTest(kprop, z)
            if sp1 > 1 or sp2 > 1: self.molar_properties(kprop, z)
            else: #tiny manipulation
                self.x = z; self.y = z
                self.bubble_point_pressure()
                if self.P > self.Pbubble: self.L = 1; self.V = 0
                else: self.L = 0; self.V = 1
                lnphil = self.lnphi_based_on_deltaG(kprop, self.x, 1)
                lnphiv = self.lnphi_based_on_deltaG(kprop, self.y, 0)
                self.fv = np.exp(lnphiv) * (self.y * self.P)
                self.fl = np.exp(lnphil) * (self.x * self.P)
                self.K = self.y/self.x


        self.z = self.x * self.L + self.y * self.V
        self.Mw_L, self.ksi_L, self.rho_L = self.other_properties(kprop, self.x, 1)
        self.Mw_V, self.ksi_V, self.rho_V = self.other_properties(kprop, self.y, 0)

        #if self.V == 0: self.ksi_V = self.V * self.ksi_V; self.rho_V = self.V * self.ksi_V

    def equilibrium_ratio_Wilson(self, kprop):
        self.K = np.exp(5.37 * (1 + kprop.w) * (1 - kprop.Tc / self.T)) * \
                (kprop.Pc / self.P)

    def coefficientsPR(self, kprop, l):
        #l - any phase molar composition
        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.54226, 0.26992])

        k = (PR_kC7[0] + PR_kC7[1] * kprop.w - PR_kC7[2] * kprop.w ** 2 + \
            PR_kC7[3] * kprop.w ** 3) * kprop.C7 + (PR_k[0] + PR_k[1] * kprop.w - \
            PR_k[2] * kprop.w ** 2) * (1 - kprop.C7)
        alpha = (1 + k * (1 - (self.T / kprop.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (kprop.R* kprop.Tc) ** 2 / kprop.Pc * alpha
        self.b = 0.07780 * kprop.R* kprop.Tc / kprop.Pc
        aalpha_i_reshape = np.ones((kprop.Nc,kprop.Nc)) * aalpha_i[:,np.newaxis]
        aalpha_ij = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - kprop.Bin)
        self.bm = sum(l * self.b)
        B = self.bm * self.P / (kprop.R* self.T)
        l_reshape = np.ones((aalpha_ij).shape) * l[:, np.newaxis]
        self.aalpha = (l_reshape.T * l[:,np.newaxis] * aalpha_ij).sum()
        A = self.aalpha * self.P / (kprop.R* self.T) ** 2
        self.psi = (l_reshape * aalpha_ij).sum(axis = 0)
        return A, B


    def Z_PR(B, A, ph):
        # PR cubic EOS: Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = [1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)]
        Z = np.roots(coef)
        roots = np.isreal(Z)
        Z_reais = np.real(Z[roots]) #Saving the real roots
        Z_ans = min(Z_reais) * ph + max(Z_reais) * (1 - ph)

        ''' This last line, considers that the phase is composed by a pure
         component, so the EOS model can return more than one real root.
            If liquid, Zl = min(Z) and gas, Zv = max(Z).
            You can notice that, if there's only one real root,
        it works as well.'''
        return Z_ans

    def lnphi(self, kprop, l, ph):
        #l - any phase molar composition
        A, B = self.coefficientsPR(kprop, l)
        Z = StabilityCheck.Z_PR(B, A, ph)
        lnphi = self.b / self.bm * (Z - 1) - np.log(Z - B) - A / (2 * (2 ** (1/2))
                * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * np.log((Z + (1 +
                2 ** (1/2)) * B) / (Z + (1 - 2 ** (1/2)) * B))
        return lnphi

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
        lnphiz = self.lnphi(kprop, z, 0)
        while max(abs(Y / Yold - 1)) > 1e-9: #convergência
            Yold = np.copy(Y)
            lnphiy = self.lnphi(kprop, y, 1)
            Y = np.exp(np.log(z) + lnphiz - lnphiy)
            y = Y / sum(Y)
        stationary_point1 = sum(Y)

    #*****************************Test two**********************************#
        #Used alone when the phase investigated (z) is clearly liquid like (ph == 1)

        Y = self.K * z
        Y_old = 0.9 * Y
        y = Y / sum(Y)
        lnphiz = self.lnphi(kprop, z, 1)
        while max(abs(Y / Y_old - 1)) > 1e-9:
            Y_old = np.copy(Y)
            lnphiy = self.lnphi(kprop, y, 0)
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
        lnphi = [self.lnphi(kprop, l, 1 - ph), self.lnphi(kprop, l, ph)]
        deltaG_molar = sum(l * (lnphi[1 - ph] - lnphi[ph]))
        if deltaG_molar >= 0: ph = ph
        else: ph = 1 - ph
        return ph

    def lnphi_based_on_deltaG(self,kprop, l, ph):
        ph = self.deltaG_molar(kprop, l, ph)
        return self.lnphi(kprop, l, ph)

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
            xi = self.solve_objective_function_Yinghui(z1, zi, z, K1, KNc, Ki, i)
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

    def other_properties(self, kprop, l, ph):
        s = 0.085
        #l - any phase molar composition
        A, B = self.coefficientsPR(kprop, l)
        ph = self.deltaG_molar(kprop, l, ph)
        Z = StabilityCheck.Z_PR(B, A, ph)
        #v = Z*kprop.R*self.T/self.P - s*self.b
        ksi_phase = self.P / (Z * kprop.R* self.T)
        Mw_phase = sum(l * kprop.Mw)
        rho_phase = ksi_phase * sum(l * kprop.Mw)
        # se precisar retornar mais coisa, entra aqui
        return Mw_phase, ksi_phase, rho_phase

    def bubble_point_pressure(self):
        #Isso vem de uma junção da Lei de Dalton com a Lei de Raoult
        Pv = self.K * self.P
        self.Pbubble = sum(self.x * Pv)

    def TPD(self, z): #ainda não sei onde usar isso
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

            '''O modo que x varia implica no formato de TPD. No presente exemplo,
            a fração molar do segundo componente de x varia direto com t, que é a
            variável de plotagem. Logo, a distancia dos planos tangentes será
            zero em z[Nc-1]. O contrário ocorreria'''
            lnphix = self.lnphi(x, 0); #new phase (vapor- ph=2)
            for j in range(0,self.Nc):
                fix = math.exp(lnphix[j]) * x[j] * self.P
                fiz = math.exp(lnphiz[j]) * z[j] * self.P
                aux = aux + x[j] * kprop.R* self.T * (math.log(fix / fiz))
                TPD[i] = aux

        plt.figure(0)
        plt.plot(t, TPD)
        plt.xlabel('x')
        plt.ylabel('TPD')
        plt.show()
        return TPD
