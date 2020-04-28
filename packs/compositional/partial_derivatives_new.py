import numpy as np
from .stability_check import StabilityCheck
import math
from scipy.misc import derivative
import sympy
from sympy.utilities import lambdify
from .properties_calculation import PropertiesCalc

# FALTA CHECAR TUDO COM A FORMA ANALÍTICA DA DERIVADA
class PartialDerivatives:

    def __init__(self, fprop):
        self.n_phases = 2
        self.T = fprop.T
        self.P = fprop.P

    def coefficientsPR(self, kprop, l):
        #l - any phase molar composition
        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.54226, 0.26992])

        k = (PR_kC7[0] + PR_kC7[1] * kprop.w - PR_kC7[2] * kprop.w ** 2 + \
            PR_kC7[3] * kprop.w ** 3) * kprop.C7 + (PR_k[0] + PR_k[1] * kprop.w - \
            PR_k[2] * kprop.w ** 2) * (1 - kprop.C7)
        alpha = (1 + k * (1 - (self.T / kprop.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (kprop.R* kprop.Tc) ** 2 / kprop.Pc * alpha
        b = 0.07780 * kprop.R* kprop.Tc / kprop.Pc
        aalpha_i_reshape = np.ones((kprop.Nc,kprop.Nc)) * aalpha_i[:,np.newaxis]
        aalpha_ij = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - kprop.Bin)
        bm = sum(l * b)
        l_reshape = np.ones((aalpha_ij).shape) * l[:, np.newaxis]
        aalpha = (l_reshape.T * l[:,np.newaxis] * aalpha_ij).sum()
        return bm, aalpha

    def PR_EOS(self, P, kprop, Sj, xkj, Nk, ph, bm ,a):
        C = np.array([P*(Sj*xkj/Nk)**3, (bm - kprop.R*self.T)*(Sj*xkj/Nk)**2, (a - 3*bm**2 - 2*kprop.R*self.T*bm)*(Sj*xkj/Nk), P*bm**3 + kprop.R*self.T*bm**2 - a*bm])
        Vt = np.roots(C)
        roots = np.isreal(Vt)
        Vt_reais = np.real(Vt[roots]) #Saving the real roots
        Vt_ans = min(Vt_reais) * ph + max(Vt_reais) * (1 - ph)
        ind = np.argwhere(Vt_ans==Vt)
        return Vt_ans

    def get_dVt_dNk(self, P, kprop, Sj, xkj, Nk, ph):
        delta = 0.0001
        bm, am = self.coefficientsPR(kprop, xkj)
        dVt_dNk = np.zeros(kprop.Nc)
        Nj = Nk/xkj
        for nc in range(kprop.Nc):
            Nk_plus = Nk[nc] + delta/2
            Nk_minus = Nk[nc] - delta/2
            Vt_plus = self.PR_EOS(P, kprop, Sj, xkj[nc], Nk_plus, 1, bm, am)
            Vt_minus = self.PR_EOS(P, kprop, Sj, xkj[nc], Nk_minus, 1, bm, am)
            dVt_dNk[nc] = (Vt_plus - Vt_minus)/delta
        return dVt_dNk

    def get_dVt_dP(self, P, kprop, Sj, xkj, Nk, ph):
        delta = 0.0001
        P_plus = P + delta/2
        P_minus = P - delta/2

        bm, am = self.coefficientsPR(kprop, xkj)

        Vt_plus = self.PR_EOS(P_plus, kprop, Sj, xkj[0], Nk[0], 1, bm, am)
        Vt_minus = self.PR_EOS(P_minus, kprop, Sj, xkj[0], Nk[0], 1, bm, am)
        dVt_dP = (Vt_plus - Vt_minus)/delta
        return dVt_dP

    def PR_EOS_sym(self, P, kprop, Sj, xkj, Nk, bm, a):
        C = np.array([P*(Sj*xkj/Nk)**3, (bm - kprop.R*self.T)*(Sj*xkj/Nk)**2, (a - 3*bm**2 - 2*kprop.R*self.T*bm)*(Sj*xkj/Nk), P*bm**3 + kprop.R*self.T*bm**2 - a*bm])
        Vt = sympy.symbols('Vt')
        a = C[0]*Vt**3 + C[1]*Vt**2 + C[2]*Vt + C[3]
        Vt_a = sympy.solve(a,Vt)
        return Vt_a

    def get_all_derivatives(self, kprop, fprop):
        n_blocks = len(fprop.P)
        dVt_dP = np.zeros(n_blocks)
        dVt_dNk = np.zeros([kprop.Nc, n_blocks])
        for b in range(n_blocks):
            P = fprop.P[b]
            So = fprop.So[b]
            xkj = fprop.component_molar_fractions[0:kprop.Nc,0,b]
            Nk = fprop.component_mole_numbers[0:kprop.Nc,b]
            dVt_dNk[:,b] = self.get_dVt_dNk(P, kprop, So, xkj, Nk, 1)
            dVt_dP[b] = self.get_dVt_dP(P, kprop, So, xkj, Nk, 1)
        return dVt_dNk, dVt_dP

    def get_dVt_dNk_symbolic(self, P, kprop, Sj, xkj, Nk, bm, am):

        Vt = self.PR_EOS_sym(P, kprop, Sj, xkj, Nk, bm, am)
        dVt_dNk = sympy.diff(Vt[2], Nk)
        dVt_dP = sympy.diff(Vt[2], P)
        dVt_dNk = lambdify([Nk, P, Sj, xkj, bm, am], dVt_dNk,'numpy')
        import pdb; pdb.set_trace()
        return dVt_dNk, dVt_dP


    def get_all_symbolic(self, kprop, fprop):
        bm = sympy.symbols('bm')
        am = sympy.symbols('am')
        P = sympy.symbols('P')
        Sj = sympy.symbols('Sj')
        xkj = sympy.symbols('xkj')
        Nk = sympy.symbols('Nk')
        n_blocks = len(fprop.P)
        dVt_dNk = np.zeros([kprop.Nc, n_blocks])
        dVt_dp = np.zeros(n_blocks)
        dVt_dNk_f, dVt_dP_f = self.get_dVt_dNk_symbolic(P, kprop, Sj, xkj, Nk, bm, am)
        for b in range(n_blocks):
            P_value = fprop.P[b]
            So = fprop.So[b]
            xkj_value = fprop.component_molar_fractions[0:kprop.Nc,0,b]
            Nk_value = fprop.component_mole_numbers[0:kprop.Nc,b]
            bm_value, am_value = self.coefficientsPR(kprop, xkj_value)
            for nc in range(kprop.Nc):
                import pdb; pdb.set_trace()
                dVt_dNk[nc,b] = dVt_dNk_f([(xkj, xkj_value[nc]), (P, P_value), (Sj, So), (Nk, Nk_value[nc]), (am, am_value), (bm, bm_value)])
            dVt_dP[b] = dVtk_dP_f([(xkj, xkj_value), (P, P_value), (Sj, So), (Nk, Nk_value), (am, am_value), (bm, bm_value)])
