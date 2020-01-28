import numpy as np
import sympy
from sympy import *
from mpmath import *


class Fugacity_parcial_derivative:
    def __init__(self, w, Bin, R, Tc, Pc, T, P, C7):
        self.w = w
        self.Bin = Bin
        self.R = R
        self.Tc = Tc
        self.Pc = Pc
        self.T = T
        self.P = P
        self.Nc = len(w)
        self.C7 = C7

    def get_dlnphidn(self, z, ph):
        l = sympy.symarray('l', len(z))
        lnphil_all = self.lnphis(l, ph)

        ''' Get Z value '''
        A, B = self.coefficientsPR(z)
        Z = StabilityCheck._Z_PR(B, A, ph)
        root = np.isreal(Z)
        real_roots_position = np.where(root == True)
        Z_reais = np.real(Z[real_roots_position[:]]) #Saving the real roots
        Z_ans = min(Z_reais) * ph + max(Z_reais) * (1 - ph)
        index_Z = np.argwhere(Z == Z_ans)

        lnphi = lnphil_all[index_Z,:]
        dlnphi_dn = np.zeros([self.Nc,self.Nc])
        for i in range(0,self.Nc):
            dlndl = sympy.diff(lnphi,l[i])
            func = lambdify(l, dlndl,'numpy')
            dlnphi_dn[i,:] = func(*z)
        print(dlnphi_dn)

    def dlnphi_dn(self,z):
        dlnphil_dn = self.get_dlnphidn(x, 1)
        dlnphiv_dn = self.get_dlnphidn(y, 0)

    def coefficientsPR(self, l):
        #l - any phase molar composition
        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.5422, 0.26992])

        k = (PR_kC7[0] + PR_kC7[1] * self.w - PR_kC7[2] * self.w ** 2 + \
            PR_kC7[3] * self.w ** 3) * self.C7 + (PR_k[0] + PR_k[1] * self.w - \
            PR_k[2] * self.w ** 2) * (1 - self.C7)
        alpha = (1 + k * (1 - (self.T / self.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (self.R * self.Tc) ** 2 / self.Pc * alpha
        self.b = 0.07780 * self.R * self.Tc / self.Pc
        aalpha_i_reshape = np.ones((self.Nc,self.Nc)) * aalpha_i[:,np.newaxis]
        aalpha_ij = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - self.Bin)
        self.bm = sum(l * self.b)
        B = self.bm * self.P / (self.R * self.T)
        l_reshape = np.ones((aalpha_ij).shape) * l[:, np.newaxis]
        self.aalpha = (l_reshape.T * l[:,np.newaxis] * aalpha_ij).sum()
        A = self.aalpha * self.P / (self.R * self.T) ** 2
        self.psi = (l_reshape * aalpha_ij).sum(axis = 0)

        return A, B

    def Zs_PR(B, A):
        # PR cubic EOS: Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = [1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)]
        a,b,c,d = coef
        x1 = (3*a*c - b**2)/3/a**2
        x2 = (2*b**3 - 9*a*b*c + 27*a**2*d)/(27*a**3)
        x3 = 18*a*b*c*d - 4*b**3*d + b**2*c**2 - 4*a*c**3 - 27*a**2*d**2
        x = np.array([x1,x2,x3])

        return x

    def lnphis(self, l):
        A, B = self.coefficientsPR(l)
        Z = StabilityCheck.Zs_PR(B, A)

        lnphi1 =  (self.b / self.bm * (Z[0] - 1) - sympy.log(Z[0] - B) - A / (2 * (2 ** (1/2))
                    * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[0] + (1 +
                    2 ** (1/2)) * B) / (Z[0] + (1 - 2 ** (1/2)) * B)))
        lnphi2 =  self.b / self.bm * (Z[1] - 1) - sympy.log(Z[1] - B) - A / (2 * (2 ** (1/2))
                    * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[1] + (1 +
                    2 ** (1/2)) * B) / (Z[1] + (1 - 2 ** (1/2)) * B))
        lnphi3 =  self.b / self.bm * (Z[2] - 1) - sympy.log(Z[2] - B) - A / (2 * (2 ** (1/2))
                    * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[2] + (1 +
                    2 ** (1/2)) * B) / (Z[2] + (1 - 2 ** (1/2)) * B))

        lnphis = np.array([lnphi1,lnphi2,lnphi3])

        return lnphis

        # dlnphi_dl = [derivative(self._lnphi,self.x[i],n = 1,args = [1,i]) for i in self.Nc]
