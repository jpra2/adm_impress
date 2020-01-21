"""Check stability of a thermodynamic equilibrium."""
import numpy as np
import math
# import matplotlib.pyplot as plt
## Encontrar os pontos estacionarios. Estes correspondem aos pontos nos quais a derivada de g com respeito a Y é 0
## Todas as equações foram confirmadas pelo livro do Dandekar e a biblioteca do thermo

class StabilityCheck:
    """Check for stability of a thermodynamic equilibrium and returns the
    equilibrium phase compositions (perform the flash calculation)."""

    def __init__(self, w, Bin, R, Tc, Pc, Vc, T, P, Mw, C7):
        self.w = np.array(w).astype(float)
        self.Bin = np.array(Bin).astype(float)
        self.Mw = np.array(Mw).astype(float)
        self.R = np.array(R).astype(float)
        self.Tc = np.array(Tc).astype(float)
        self.Pc = np.array(Pc).astype(float)
        self.Vc = np.array(Vc).astype(float)
        self.T = np.array(T).astype(float)
        self.P = np.array(P).astype(float)
        self.Nc = len(w)
        self.C7 = np.array(C7)
        #StabilityCheck.TPD(self)

    def run(self, z):
        if any(z <= 0):
            self.molar_properties(z)
        else:
            sp1,sp2 = self.Stability(z)
            # import pdb; pdb.set_trace()
            if sp1 > 1 or sp2 > 1:
                self.molar_properties(z)
            '''if sp1<1 and sp2<1:
                TPD = obj.TPD(z)
                if TPD.any()<0: #checar se isso iria funcionar
                    obj.molar_properties(z,Mw)'''

    def coefficientsPR(self):
        # The equilibrium ratio (K) is calculated by Wilson's equation

        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.5422, 0.26992])

        k = (PR_kC7[0] + PR_kC7[1] * self.w - PR_kC7[2] * self.w ** 2 + \
            PR_kC7[3] * self.w ** 3) * self.C7 + (PR_k[0] + PR_k[1] * self.w - \
            PR_k[2] * self.w ** 2) * (1 - self.C7)
        alpha = (1 + k * (1 - (self.T / self.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (self.R * self.Tc) ** 2 / self.Pc * alpha
        self.b = 0.07780 * self.R * self.Tc / self.Pc
        self.K = np.exp(5.37 * (1 + self.w) * (1 - self.Tc / self.T)) * \
                (self.Pc / self.P)
        aalpha_i_reshape = np.ones((self.Nc,self.Nc)) * aalpha_i[:,np.newaxis]
        self.aalpha_ij = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - self.Bin)


    def Z_PR(B, A, ph):
        # PR cubic EOS: Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = [1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)]
        Z = np.roots(coef)
        root = np.isreal(Z) # return True for real roots
        #position where the real roots are - crated for organization only
        real_roots_position = np.where(root == True)
        Z_reais = np.real(Z[real_roots_position[:]]) #Saving the real roots
        Z_ans = min(Z_reais) * ph + max(Z_reais) * (1 - ph)

        ''' This last line, considers that the phase is composed by a pure
         component, so the EOS model can return more than one real root.
            If liquid, Zl = min(Z) and gas, Zv = max(Z).
            You can notice that, if there's only one real root,
        it works as well.'''
        return Z_ans


    def lnphi(self, x, ph):
        bm = sum(x * self.b)
        B = bm * self.P / (self.R * self.T)
        x_reshape = np.ones((self.aalpha_ij).shape) * x[:, np.newaxis]
        aalpha = (x_reshape.T * x[:,np.newaxis] * self.aalpha_ij).sum()
        A = aalpha * self.P / (self.R * self.T) ** 2
        Z = StabilityCheck.Z_PR(B, A, ph)
        psi = (x_reshape * self.aalpha_ij).sum(axis = 0)
        lnphi = self.b / bm * (Z - 1) - np.log(Z - B) - A / (2 * (2 ** (1/2))
                * B) * (2 * psi / aalpha - self.b / bm) * np.log((Z + (1 +
                2 ** (1/2)) * B) / (Z + (1 - 2 ** (1/2)) * B))

        return lnphi

    def Stability(self,z):
        self.coefficientsPR()
    #****************************INITIAL GUESS******************************#
    ## Both approaches bellow should be used in case the phase is in the critical region

        #identifying the initial phase 'z' mole fraction
        ''' In the lnphi function: 0 stands for vapor phase and 1 for liquid'''

    #*****************************Test one**********************************#
        #Used alone when the phase investigated (z) is clearly vapor like (ph = 0)

        Y = z / self.K
        Yold = 0.9 * Y
        y = Y / sum(Y)
        lnphiz = self.lnphi(z, 0)
        while max(abs(Y / Yold - 1)) > 1e-9: #convergência
            Yold = np.copy(Y)
            lnphiy = self.lnphi(y, 1)
            Y = np.exp(np.log(z) + lnphiz - lnphiy)
            y = Y / sum(Y);

        stationary_point1 = sum(Y)
        # print(sum(Y))
        # if sum(Y) <= 1: print('estavel')
        # else: print('instavel') #estavel2 = 0

        ''' If this first test returns stable: There might be a chance that the
        Gibbs free energy is at its global minimum. However, as will be
        explaned later, this alone doesn't guarantee the phase stability.
         If it returns unstable: There is another composition that makes the
        Gibbs free energy at its global minimum indicated by sum(Y)>1'''
    #*****************************Test two**********************************#
        #Used alone when the phase investigated (z) is clearly liquid like (ph == 1)

        Y = self.K * z
        Y_old = 0.9 * Y
        y = Y / sum(Y)
        lnphiz = self.lnphi(z, 1)
        while max(abs(Y / Y_old - 1)) > 1e-9:
            Y_old = np.copy(Y)
            lnphiy = self.lnphi(y, 0)
            Y = np.exp(np.log(z) + lnphiz - lnphiy)
            y = Y / sum(Y)

        stationary_point2 = sum(Y)
        # print(sum(Y))
        # if sum(Y) <= 1: print('estavel')#estavel2 = 1
        # else: print('instavel') #estavel2 = 0

        """
        The same thing happens here. The difference is that, the original
        phase is gas, and then the "new" phase is supposed to be liquid.
        In cases that the compressibility equation returns only one root,
        both tests work like two different initial guess for the same problem,
        being more likely to find a stationary point that makes the phase
        unstable.
        """

        ''' If one of these approaches returns unstable the system is unstable.
        The stability of the phase is something a little bit more complex
        to guarantee. '''

        return stationary_point1,stationary_point2

    """-------------Below starts biphasic flash calculations-----------------"""

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
        zi = np.delete(zi, indexNc - 1)

        #starting x
        self.x = np.zeros(self.Nc)

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
        for i in range(0, len(xi)):
            self.x[self.K == Ki[i]] = xi[i]

        self.y = self.K * self.x

    def molar_properties_Yinghui(self,z):
        #razao = fl/fv -> an arbitrary vector to enter in the iterative mode
        razao = np.ones(self.Nc)/2
        while max(abs(razao - 1)) > 1e-9:
            self.Yinghui_method(z)

            lnphix = [self.lnphi(self.x, 0), self.lnphi(self.x, 1)]
            lnphiy = [self.lnphi(self.y, 0), self.lnphi(self.y, 1)]
            deltaGx_molar = sum(self.x * (lnphix[0] - lnphix[1]))
            deltaGy_molar = sum(self.y * (lnphiy[1] - lnphiy[0]))
            lnphil = -lnphix[0] * np.sign(1 - np.sign(deltaGx_molar)) * \
                    np.sign(deltaGx_molar) + lnphix[1] * np.sign(1 +
                    np.sign(deltaGx_molar))
            lnphiv = -lnphiy[1] * np.sign(1 - np.sign(deltaGy_molar)) * \
                    np.sign(deltaGy_molar) + lnphiy[0] * np.sign(1 +
                    np.sign(deltaGy_molar))

            """
            This part above means that:
            if deltaGx_molar>=0: lnphil = lnphix[1]
                                (lnphi for x in the liquid phase)
            if deltaGx_molar<0: lnphil = lnphix[0]
                                (lnphi for x in the vapour phase)
                                because it minimizes the Gibbs molar energy
            The same analysis is made for the y phase composition.
            """
            self.fl = np.exp(lnphil) * (self.x * self.P)
            self.fv = np.exp(lnphiv) * (self.y * self.P)

            razao = np.divide(self.fl, self.fv, out = razao / razao * (1 + 1e-10),
                              where = self.fv != 0)
            self.K = razao * self.K


    def solve_objective_function_Whitson(self, z):
        """ Solving for V """
        #K1 = max(self.K); KNc = min(self.K)
        Vmax = max(1, 1 / (1 - min(self.K)))#1 / (1 - min(self.K))#
        Vmin = min(0, 1 / (1 - max(self.K))) # 1 / (1 - max(self.K))

        #Vmin = ((K1-KNc)*z[self.K==K1]-(1-KNc))/((1-KNc)*(K1-1))
        #proposed by Li et al for Whitson method
        V = (Vmin + Vmax) / 2
        Vold = V / 2 #just to get into the loop

        while abs(V / Vold - 1) > 1e-8:
            Vold = V
            f = sum((self.K - 1) * z / (1 + V * (self.K - 1)))
            df = -sum((self.K - 1) ** 2 * z / (1 + V * (self.K - 1)) ** 2)
            # import pdb; pdb.set_trace()
            V = V - f / df #Newton-Raphson iterative method
            if V > Vmax: V = Vmax #(Vold + Vmax)/2
            elif V < Vmin: V = Vmin #(Vold + Vmin)/2

        self.x = z / (1 + V * (self.K - 1))
        self.y = self.K * self.x


    def molar_properties_Whitson(self, z):

        while max(abs(self.fv / self.fl - 1)) > 1e-9:
            self.solve_objective_function_Whitson(z)
            # In order to verify the state of each new phase, in consequence,
            #the new phases stabilities, we have to identify the composition
            #that makes the Gibbs free energy smaller
            lnphix = [self.lnphi(self.x, 0), self.lnphi(self.x, 1)]
            lnphiy = [self.lnphi(self.y, 0), self.lnphi(self.y, 1)]
            deltaGx_molar = sum(self.x * (lnphix[0] - lnphix[1]))
            deltaGy_molar = sum(self.y * (lnphiy[1] - lnphiy[0]))
            lnphil = -lnphix[0] * np.sign(1 - np.sign(deltaGx_molar)) * \
                    np.sign(deltaGx_molar) + lnphix[1] * np.sign(1 +
                    np.sign(deltaGx_molar))
            lnphiv = -lnphiy[1] * np.sign(1 - np.sign(deltaGy_molar)) * \
                    np.sign(deltaGy_molar) + lnphiy[0] * np.sign(1 +
                    np.sign(deltaGy_molar))
            self.fv = np.exp(lnphiv) * (self.y * self.P)
            self.fl = np.exp(lnphil) * (self.x * self.P)
            self.K = (self.fl / self.fv) * self.K

    def molar_properties(self,z):
        self.coefficientsPR()
        # Aqui a correlação de wilson é utilizada apenas para achar o K inicial
        self.fv = 2 * np.ones(self.Nc); self.fl = np.ones(self.Nc) #entrar no primeiro loop

        if self.Nc <= 2: self.molar_properties_Whitson(z)
        else: self.molar_properties_Yinghui(z)

        ''' Molar Volume Fractions '''
        self.V = (z[self.x != 0] - self.x[self.x != 0]) / (self.y[self.x != 0] -
                self.x[self.x != 0])
        self.L = 1 - self.V

        ''' Phase Molecular Weight '''
        self.Mw_L = sum(self.x * self.Mw)
        self.Mw_V = sum(self.y * self.Mw)

        ''' Phase Mass Densities '''
        if self.V!=0:
            self.rho_V = self.Mw_V / self.V
        else: self.rho_V = 0
        if self.L!=0:
            self.rho_L = self.Mw_L / self.L
        else: self.rho_L = 0
        
        ''' Phase molar densities '''
        self.eta_L = self.rho_L/self.Mw_L
        self.eta_V = self.rho_V/self.Mw_V

        # se precisar retornar mais coisa, entra aqui



    def TPD(self, z): #ainda não sei onde usar isso
        x = np.zeros(self.Nc)

        #**********************Tangent Plane distance plot*********************#
        t = np.linspace(0.01, 0.99, 0.9 / 0.002) #vetor auxiliar
        TPD = np.zeros(len(t)) ##F

        for i in range(0, len(t)):
            aux = 0;
            lnphiz = self.lnphi(z, 1) #original phase

            #x = np.array([1-t[i],t[i]]) #new phase composition (1-t e t) - apenas válido para Nc=2 acredito eu.
            for k in range(0, self.Nc - 1):
                x[k] = (1 - t[i]) / (self.Nc - 1)
                x[self.Nc - 1] = t[i]

            '''O modo que x varia implica no formato de TPD. No presente exemplo,
            a fração molar do segundo componente de x varia direto com t, que é a
            variável de plotagem. Logo, a distancia dos planos tangentes será
            zero em z[Nc-1]. O contrário ocorreria'''
            lnphix = self.lnphi(x, 0); #new phase (vapor- ph=2)
            for j in range(0,self.Nc):
                fix = math.exp(lnphix[j]) * x[j] * self.P
                fiz = math.exp(lnphiz[j]) * z[j] * self.P
                aux = aux + x[j] * self.R * self.T * (math.log(fix / fiz))
                TPD[i] = aux

        plt.figure(0)
        plt.plot(t, TPD)
        plt.xlabel('x')
        plt.ylabel('TPD')
        plt.show()
        return TPD
