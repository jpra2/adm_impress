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
        self.T = fprop.T

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
        aalpha_ik = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - kprop.Bin)
        bm = np.sum(l * b[:,np.newaxis], axis=0)
        l_reshape = np.ones((kprop.Nc, kprop.Nc, len(l[0,:]))) * l[:, np.newaxis, :]
        aalpha = (l_reshape.T * l.T[:, :, np.newaxis] * aalpha_ik[np.newaxis,:,:]).sum(axis=2).sum(axis=1)
        return bm, aalpha

    def coefficientsPR2(self, kprop, l):
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
        aalpha_ik = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - kprop.Bin)
        bm = sum(l * b)
        l_reshape = np.ones((aalpha_ik).shape) * l[:, np.newaxis]
        aalpha = (l_reshape.T * l[:,np.newaxis] * aalpha_ik).sum()

        return bm, aalpha

    def get_dVt_dNk_analytically(self, P, kprop, Vt, Sj, xkj, Nk):
        bm, am = self.coefficientsPR(kprop, xkj)
        C = np.array([P*(Sj*xkj)**3, (bm*P - kprop.R*self.T)*(Sj*xkj)**2, (am - 3*P*bm**2 - 2*kprop.R*self.T*bm)*(Sj*xkj)])

        Num = 3*Vt**3/(Nk**4)*C[0] + 2*Vt**2/(Nk**3)*C[1] + Vt/(Nk**2)*C[2]
        Den = 3*Vt**2/(Nk**3)*C[0] + 2*Vt/(Nk**2)*C[1] + 1/Nk*C[2]

        dVt_dNk = Num/Den
        return dVt_dNk

    def get_dVt_dP_analytically(self, P, kprop, Vt, Nj, xkj):
        bm, am = self.coefficientsPR(kprop, xkj)
        C = np.array([(1/Nj)**3, (1/Nj)**2, (1/Nj)])

        Num = - Vt**3*C[0] - bm**3 - Vt**2*bm*C[1] + 3*bm**2*Vt*C[2]
        Den = 3*Vt**2*P*C[0] + 2*Vt*(bm*P - kprop.R*self.T)*C[1] + (am - 3*bm**2*P - 2*kprop.R*self.T*bm)*C[2]

        dVt_dP = Num/Den
        return dVt_dP

    def get_all_derivatives(self, kprop, fprop):
        n_blocks = len(fprop.P)
        #dVt_dP = np.zeros(n_blocks)
        #dVt_dNk = np.zeros([kprop.Nc, n_blocks])

        P = fprop.P_old
        So = fprop.So
        Vt = fprop.Vt
        xkj = fprop.component_molar_fractions[0:kprop.Nc,0,:]
        Nk = fprop.component_mole_numbers[0:kprop.Nc,:]
        dVt_dNk = self.get_dVt_dNk_analytically(P, kprop, Vt, So, xkj, Nk)
        No = fprop.phase_mole_numbers[0,0,:]
        dVt_dP = self.get_dVt_dP_analytically(P, kprop, Vt, No, xkj)

        '''for b in range(n_blocks):
            P = fprop.P_old[b]
            So = fprop.So[b]
            Vt = fprop.Vt[b]
            xkj = fprop.component_molar_fractions[0:kprop.Nc,0,b]
            Nk = fprop.component_mole_numbers[0:kprop.Nc,b]
            dVt_dNk[:,b] = self.get_dVt_dNk(P, kprop, So, xkj, Nk, 1)
            dVt_dP[b] = self.get_dVt_dP(P, kprop, So, xkj, Nk, 1)'''

        return dVt_dNk, dVt_dP

    def PR_EOS(self, P, kprop, Sj, xkj, Nk, ph, bm ,a):
        C = np.array([P*(Sj*xkj/Nk)**3, (bm*P - kprop.R*self.T)*(Sj*xkj/Nk)**2, (a - 3*P*bm**2 - 2*kprop.R*self.T*bm)*(Sj*xkj/Nk), P*bm**3 + kprop.R*self.T*bm**2 - a*bm])
        Vt = np.roots(C)
        roots = np.isreal(Vt)
        Vt_reais = np.real(Vt[roots]) #Saving the real roots
        Vt_ans = min(Vt_reais) * ph + max(Vt_reais) * (1 - ph)
        ind = np.argwhere(Vt_ans==Vt)
        return Vt_ans

    def get_dVt_dNk(self, P, kprop, Sj, xkj, Nk, ph):
        delta = 0.0001
        bm, am = self.coefficientsPR2(kprop, xkj)
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

        bm, am = self.coefficientsPR2(kprop, xkj)

        Vt_plus = self.PR_EOS(P_plus, kprop, Sj, xkj[0], Nk[0], 1, bm, am)
        Vt_minus = self.PR_EOS(P_minus, kprop, Sj, xkj[0], Nk[0], 1, bm, am)
        dVt_dP = (Vt_plus - Vt_minus)/delta
        return dVt_dP

    '''def get_all_derivatives(self, kprop, fprop):
        n_blocks = len(fprop.P)
        dVt_dP = np.zeros(n_blocks)
        dVt_dNk = np.zeros([kprop.Nc, n_blocks])
        for b in range(n_blocks):
            P = fprop.P[b]
            So = fprop.So[b]
            Vt = fprop.Vt[b]
            xkj = fprop.component_molar_fractions[0:kprop.Nc,0,b]
            Nk = fprop.component_mole_numbers[0:kprop.Nc,b]
            dVt_dNk[:,b] = self.get_dVt_dNk_analytically(P, kprop, Vt, So, xkj, Nk)
            #dVt_dNkk = self.get_dVt_dNk(P, kprop, So, xkj, Nk, 1)
            dVt_dP[b] = self.get_dVt_dP_analytically(P, kprop, Vt, So, xkj, Nk)'''
            #dVt_dPp = self.get_dVt_dP(P, kprop, So, xkj, Nk, 1)
