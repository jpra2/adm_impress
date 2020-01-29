import numpy as np
import sympy
from sympy import *
from mpmath import *
from .stability_check import StabilityCheck

class FugacityParcialDerivative:

    def get_dlnphi_dnk(self, fluid_properties, z, Nphase, ph):
        nkphase = sympy.symarray('nkphase', len(z))
        # Nphase = self.Nphase[ph]
        lnphi_all = self.lnphi_sym(fluid_properties, nkphase, Nphase, ph)

        A, B = self.coefficientsPR(fluid_properties, z)
        Z = np.array(FugacityParcialDerivative.Z_PR_sym(B, A, ph))
        reais = np.full(len(Z), False, dtype=bool)
        for r in range(0,len(Z)):
            reais[r] = sympify(Z[r]).is_real
        Z_reais = Z[reais]
        Z_ans = min(Z_reais) * ph + max(Z_reais) * (1 - ph)

        a = np.argwhere(Z == Z_ans)
        lnphi = lnphi_all[a,:]
        dlnphi_dnk = np.zeros([fluid_properties.Nc,fluid_properties.Nc])

        for i in range(0,fluid_properties.Nc):
            dlnphi_dn = sympy.diff(lnphi,nkphase[i])
            import pdb; pdb.set_trace()
            func = lambdify(nkphase, dlnphi_dn,'numpy')
            dlnphi_dnk[:,i] = np.array(func(*z)[0][0])
        return dlnphi_dnk

    def get_dlnphi_dnk_numerically(self, fluid_properties, l, Nphase, ph):
        nkphase = l * Nphase
        delta = 0.0001
        dlnphi_dnk = np.zeros([fluid_properties.Nc,fluid_properties.Nc])

        for i in range(0,fluid_properties.Nc):
            nkphased = np.copy(nkphase)
            nkphased[i] = nkphase[i] + delta
            dlnphi_dnk[:,i] = (self.lnphi_sym2(fluid_properties, nkphased, Nphase, ph) - self.lnphi_sym2(fluid_properties, nkphase, Nphase, ph))/delta
        return dlnphi_dnk

    def dlnphi_dn_matrix(self, fluid_properties, Nphase):
        dlnphil_dn_all = self.get_dlnphi_dnk(fluid_properties, fluid_properties.x, Nphase[1], 1)
        dlnphiv_dn_all = self.get_dlnphi_dnk(fluid_properties, fluid_properties.y, Nphase[0], 0)
        matrix = dlnphil_dn_all + dlnphiv_dn_all
        return matrix

    def dlnphi_dn_matrix_numerically(self, fluid_properties, Nphase):
        dlnphil_dn_all = self.get_dlnphi_dnk_numerically(fluid_properties, fluid_properties.x, Nphase[1], 1)
        dlnphiv_dn_all = self.get_dlnphi_dnk_numerically(fluid_properties, fluid_properties.y, Nphase[0], 0)
        matrix = dlnphil_dn_all + dlnphiv_dn_all
        return matrix

    def dlnphi_dn_matrix(self, fluid_properties, Nphase):
        dlnphil_dn_all = self.get_dlnphi_dnk(fluid_properties, fluid_properties.x, Nphase[1], 1)
        dlnphiv_dn_all = self.get_dlnphi_dnk(fluid_properties, fluid_properties.y, Nphase[0], 0)
        matrix = dlnphil_dn_all + dlnphiv_dn_all
        return matrix

    def Z_PR_sym(B, A, ph):
        Z = symbols('Z')
        a = Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = [1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)]
        Za = solve(a,Z)
        return Za

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

    def coefficientsPR(self, fluid_properties, l):
        #l - any phase molar composition
        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.5422, 0.26992])

        k = (PR_kC7[0] + PR_kC7[1] * fluid_properties.w - PR_kC7[2] * fluid_properties.w ** 2 + \
            PR_kC7[3] * fluid_properties.w ** 3) * fluid_properties.C7 + (PR_k[0] + PR_k[1] * fluid_properties.w - \
            PR_k[2] * fluid_properties.w ** 2) * (1 - fluid_properties.C7)
        alpha = (1 + k * (1 - (fluid_properties.T / fluid_properties.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (fluid_properties.R * fluid_properties.Tc) ** 2 / fluid_properties.Pc * alpha
        self.b = 0.07780 * fluid_properties.R * fluid_properties.Tc / fluid_properties.Pc
        aalpha_i_reshape = np.ones((fluid_properties.Nc,fluid_properties.Nc)) * aalpha_i[:,np.newaxis]
        aalpha_ij = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - fluid_properties.Bin)
        self.bm = sum(l * self.b)
        B = self.bm * fluid_properties.P / (fluid_properties.R * fluid_properties.T)
        l_reshape = np.ones((aalpha_ij).shape) * l[:, np.newaxis]
        self.aalpha = (l_reshape.T * l[:,np.newaxis] * aalpha_ij).sum()
        A = self.aalpha * fluid_properties.P / (fluid_properties.R * fluid_properties.T) ** 2
        self.psi = (l_reshape * aalpha_ij).sum(axis = 0)
        return A, B

    def lnphi_sym(self, fluid_properties, nkphase, Nphase, ph):
        if Nphase == 0:
            lnphi = np.zeros(fluid_properties.Nc)
        else:
            l = nkphase/Nphase
            A, B = self.coefficientsPR(fluid_properties, l)
            Z = FugacityParcialDerivative.Z_PR_sym(B, A, ph)

            lnphi1 =  (self.b / self.bm * (Z[0] - 1) - sympy.log(Z[0] - B) - A / (2 * (2 ** (1/2))
                        * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[0] + (1 +
                        2 ** (1/2)) * B) / (Z[0] + (1 - 2 ** (1/2)) * B)))
            lnphi2 =  self.b / self.bm * (Z[1] - 1) - sympy.log(Z[1] - B) - A / (2 * (2 ** (1/2))
                        * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[1] + (1 +
                        2 ** (1/2)) * B) / (Z[1] + (1 - 2 ** (1/2)) * B))
            lnphi3 =  self.b / self.bm * (Z[2] - 1) - sympy.log(Z[2] - B) - A / (2 * (2 ** (1/2))
                        * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[2] + (1 +
                        2 ** (1/2)) * B) / (Z[2] + (1 - 2 ** (1/2)) * B))

            lnphi = np.array([lnphi1,lnphi2,lnphi3])
        return lnphi

    def lnphi_sym2(self, fluid_properties, nkphase, Nphase, ph):
        if Nphase == 0: l = np.zeros(len(nkphase))
        else: l = nkphase/Nphase
        A, B = self.coefficientsPR(fluid_properties, l)
        Z = FugacityParcialDerivative.Z_PR(B, A, ph)

        lnphi = self.b / self.bm * (Z - 1) - np.log(Z - B) - A / (2 * (2 ** (1/2))
                * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * np.log((Z + (1 +
                2 ** (1/2)) * B) / (Z + (1 - 2 ** (1/2)) * B))
        return lnphi
