import numpy as np
from cmath import acos
from ..utils import constants as ctes
from ..solvers.EOS_solver.solver import CubicRoots

class PengRobinson:
    def __init__(self, P, T, kprop):
        self.P = P
        self.T = T
        self.coefficientsPR(kprop)

    def coefficientsPR(self, kprop):
        #l - any phase molar composition
        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.54226, 0.26992])

        k = (PR_kC7[0] + PR_kC7[1] * kprop.w - PR_kC7[2] * kprop.w ** 2 + \
            PR_kC7[3] * kprop.w ** 3) * (1*(kprop.w >= 0.49))  + (PR_k[0] + PR_k[1] * kprop.w - \
            PR_k[2] * kprop.w ** 2) * (1*(kprop.w < 0.49))
        alpha = (1 + k * (1 - (self.T/ kprop.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (ctes.R * kprop.Tc) ** 2 / kprop.Pc * alpha
        self.b = 0.07780 * ctes.R * kprop.Tc / kprop.Pc
        aalpha_i_reshape = np.ones((kprop.Nc,kprop.Nc)) * aalpha_i[:,np.newaxis]
        self.aalpha_ij = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - kprop.Bin)

    '''def coefficients_cubic_EOS(self, kprop, l):
        self.bm = sum(l * self.b)
        l_reshape = np.ones((self.aalpha_ij).shape) * l[:, np.newaxis]
        self.aalpha = (l_reshape.T * l[:,np.newaxis] * self.aalpha_ij).sum()
        B = self.bm * self.P/ (ctes.R* self.T)
        A = self.aalpha * self.P/ (ctes.R* self.T) ** 2
        self.psi = (l_reshape * self.aalpha_ij).sum(axis = 0)
        return A, B

    def Z(B, A, ph):
        # PR cubic EOS: Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = [1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)]
        Z = np.roots(coef)
        root = np.isreal(Z) # return True for real roots
        #position where the real roots are - crated for organization only
        real_roots_position = np.where(root == True)
        Z_reais = np.real(Z[real_roots_position[:]]) #Saving the real roots
        Z = min(Z_reais) * ph + max(Z_reais) * (1 - ph)
        return Z_reais'''
    ''' This last line, considers that the phase is composed by a pure
         component, so the EOS model can return more than one real root.
            If liquid, Zl = min(Z) and gas, Zv = max(Z).
            You can notice that, if there's only one real root,
    it works as well.'''


    def lnphi(self, kprop, l, ph):
        #l - any phase molar composition
        l = l[:,np.newaxis]
        A, B = self.coefficients_cubic_EOS_vectorized(kprop,l)
        Z = PengRobinson.Z_vectorized(A, B)
        Z = min(Z) * ph + max(Z) * (1 - ph)
        lnphi = self.b / self.bm * (Z - 1) - np.log(Z - B) - A / (2 * (2 ** (1/2))
                * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * np.log((Z + (1 +
                2 ** (1/2)) * B) / (Z + (1 - 2 ** (1/2)) * B))

        return lnphi

    """ Bellow is showed two functions of the PR EOS equation, that were constructed in a vectorized manner.
     I modified the code so this two functions are used in the vectorized and non vectorized part. For now they
     are working in both of them. But, just because I didn't tested everything possible, I don't completelly
     trust them. Thats why they are kind of separeted from the original ones(that I do trust completelly but
     you can't calculate in a vectorized manner with them."""

    def coefficients_cubic_EOS_vectorized(self, kprop, l):
        self.bm = np.sum(l * self.b[:,np.newaxis], axis=0)
        l_reshape = np.ones((self.aalpha_ij).shape)[:,:,np.newaxis] * l[:,np.newaxis,:]
        self.aalpha = (l_reshape * l[np.newaxis,:,:] * self.aalpha_ij[:,:,np.newaxis]).sum(axis=0).sum(axis=0)
        B = self.bm * self.P / (ctes.R* self.T)
        A = self.aalpha * self.P / (ctes.R* self.T) ** 2
        self.psi = (l_reshape * self.aalpha_ij).sum(axis = 0)
        return A, B

    def Z_vectorized(A, B):
        coef = np.empty([4,len(B.ravel())])
        coef[0,:] = np.ones(len(B))
        coef[1,:] = -(1 - B)
        coef[2,:] = (A - 2*B - 3*B**2)
        coef[3,:] = -(A*B - B**2 - B**3)
        Z = CubicRoots().run(coef)
        return Z

    """ Derivatives - Still need to organize this"""

    def get_dVt_dNk_analytically(self, P, kprop, Vt, Sj, xkj, Nk):
        self.coefficients_cubic_EOS_vectorized(kprop, xkj)
        bm = self.bm; am = self.aalpha
        C = np.array([P*(Sj*xkj)**3, (bm*P - ctes.R*self.T)*(Sj*xkj)**2, (am - 3*P*bm**2 - 2*ctes.R*self.T*bm)*(Sj*xkj)])

        Num = 3*Vt**3/(Nk**4)*C[0] + 2*Vt**2/(Nk**3)*C[1] + Vt/(Nk**2)*C[2]
        Den = 3*Vt**2/(Nk**3)*C[0] + 2*Vt/(Nk**2)*C[1] + 1/Nk*C[2]

        dVt_dNk = Num/Den
        return dVt_dNk

    def get_dVt_dP_analytically(self, P, kprop, Vt, Nj, xkj):
        self.coefficients_cubic_EOS_vectorized(kprop, xkj)
        bm = self.bm; am = self.aalpha
        C = np.array([(1/Nj)**3, (1/Nj)**2, (1/Nj)])

        Num = - Vt**3*C[0] - bm**3 - Vt**2*bm*C[1] + 3*bm**2*Vt*C[2]
        Den = 3*Vt**2*P*C[0] + 2*Vt*(bm*P - ctes.R*self.T)*C[1] + (am - 3*bm**2*P - 2*ctes.R*self.T*bm)*C[2]

        dVt_dP = Num/Den
        return dVt_dP

    def get_all_derivatives(self, kprop, fprop):
        n_blocks = len(fprop.P)

        P = fprop.P
        So = fprop.So
        Vt = fprop.Vt
        xkj = fprop.component_molar_fractions[0:kprop.Nc,0,:]
        Nk = fprop.component_mole_numbers[0:kprop.Nc,:]
        dVt_dNk = self.get_dVt_dNk_analytically(P, kprop, Vt, So, xkj, Nk)
        No = fprop.phase_mole_numbers[0,0,:]
        dVt_dP = self.get_dVt_dP_analytically(P, kprop, Vt, No, xkj)
        return dVt_dNk, dVt_dP
