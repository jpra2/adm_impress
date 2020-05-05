import numpy as np


class PengRobinson:
    def __init__(self, P, T, kprop):
        self.P = P
        self.T = T
        self.coefficientsPR(kprop)

    def coefficientsPR(self, kprop):
        #l - any phase molar composition
        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.5422, 0.26992])

        k = (PR_kC7[0] + PR_kC7[1] * kprop.w - PR_kC7[2] * kprop.w ** 2 + \
            PR_kC7[3] * kprop.w ** 3) * kprop.C7 + (PR_k[0] + PR_k[1] * kprop.w - \
            PR_k[2] * kprop.w ** 2) * (1 - kprop.C7)
        alpha = (1 + k * (1 - (self.T/ kprop.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (kprop.R * kprop.Tc) ** 2 / kprop.Pc * alpha
        self.b = 0.07780 * kprop.R * kprop.Tc / kprop.Pc
        aalpha_i_reshape = np.ones((kprop.Nc,kprop.Nc)) * aalpha_i[:,np.newaxis]
        self.aalpha_ij = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - kprop.Bin)

    def coefficients_cubic_EOS(self, kprop, l):
        self.bm = sum(l * self.b)
        l_reshape = np.ones((self.aalpha_ij).shape) * l[:, np.newaxis]
        self.aalpha = (l_reshape.T * l[:,np.newaxis] * self.aalpha_ij).sum()
        B = self.bm * self.P/ (kprop.R* self.T)
        A = self.aalpha * self.P/ (kprop.R* self.T) ** 2
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
        Z_ans = min(Z_reais) * ph + max(Z_reais) * (1 - ph)

        ''' This last line, considers that the phase is composed by a pure
         component, so the EOS model can return more than one real root.
            If liquid, Zl = min(Z) and gas, Zv = max(Z).
            You can notice that, if there's only one real root,
        it works as well.'''
        return Z_ans

    def lnphi(self, kprop, l, ph):
        #l - any phase molar composition
        A, B = self.coefficients_cubic_EOS(kprop,l)
        Z = PengRobinson.Z(B, A, ph)
        lnphi = self.b / self.bm * (Z - 1) - np.log(Z - B) - A / (2 * (2 ** (1/2))
                * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * np.log((Z + (1 +
                2 ** (1/2)) * B) / (Z + (1 - 2 ** (1/2)) * B))

        return lnphi
