import numpy as np
from ..utils import constants as ctes
from scipy import linalg
from ..solvers.EOS_solver.solver import CubicRoots

class PengRobinson:
    def __init__(self, T):
        self.T = T
        self.coefficientsPR()

    def coefficientsPR(self):
        #l - any phase molar composition
        PR_kC7 = np.array([0.379642, 1.48503, 0.1644, 0.016667])
        PR_k = np.array([0.37464, 1.54226, 0.26992])
        k = (PR_kC7[0] + PR_kC7[1] * ctes.w - PR_kC7[2] * ctes.w ** 2 + \
            PR_kC7[3] * ctes.w ** 3) * (1*(ctes.w >= 0.4884))  + (PR_k[0] + PR_k[1] * ctes.w - \
            PR_k[2] * ctes.w ** 2) * (1*(ctes.w < 0.4884))
        alpha = (1 + k * (1 - (self.T/ ctes.Tc) ** (1 / 2))) ** 2
        aalpha_i = 0.45724 * (ctes.R * ctes.Tc) ** 2 / ctes.Pc * alpha
        self.b = 0.07780 * ctes.R * ctes.Tc / ctes.Pc
        aalpha_i_reshape = np.ones((ctes.Nc,ctes.Nc)) * aalpha_i[:,np.newaxis]
        self.aalpha_ik = np.sqrt(aalpha_i_reshape.T * aalpha_i[:,np.newaxis]) \
                        * (1 - ctes.Bin)

    def coefficients_cubic_EOS_vectorized(self, l, P):
        self.bm = np.sum(l * self.b[:,np.newaxis], axis = 0)
        l_reshape = np.ones((self.aalpha_ik).shape)[:,:,np.newaxis] * l[:,np.newaxis,:]
        self.aalpha = (l_reshape * l[np.newaxis,:,:] * self.aalpha_ik[:,:,np.newaxis]).sum(axis=0).sum(axis=0)
        B = self.bm * P / (ctes.R* self.T)
        A = self.aalpha * P / (ctes.R* self.T) ** 2
        self.psi = (l_reshape * self.aalpha_ik[:,:,np.newaxis]).sum(axis = 0)
        return A, B

    def Z_vectorized(self, A, B, ph):
        coef = np.empty([4,len(B.ravel())])
        coef[0,:] = 1
        coef[1,:] = -(1 - B)
        coef[2,:] = (A - 2*B - 3*B**2)
        coef[3,:] = -(A*B - B**2 - B**3)
        Z = CubicRoots().run(coef)
        root = np.isreal(Z)
        n_reais = np.sum(root*1, axis = 1)

        aux_reais = (n_reais<3) & (n_reais>1)
        if any(aux_reais): Z[~root[aux_reais]] = Z[root[aux_reais]][0]

        Z[~root[n_reais==1]] = np.repeat(Z[root[n_reais == 1]], 2)

        aux_neg = np.zeros(Z.shape,dtype=bool)
        aux_neg[Z<0] = True
        Z[aux_neg] = Z[~aux_neg][0]
        Zsave = Z
        Z = np.min(Z, axis = 1) * ph + np.max(Z, axis = 1) * (1 - ph)
        Z = np.real(Z)
        return Z

    def lnphi(self, l, P, ph):
        xkj = l
        A, B = self.coefficients_cubic_EOS_vectorized(l, P)
        Z = self.Z_vectorized(A, B, ph)
        dd = (Z[np.newaxis,:] + (1 - 2 ** (1/2)) * B[np.newaxis,:])
        lnphi = self.lnphi_calculation(A, B, Z)
        return lnphi

    def lnphi_calculation(self, A, B, Z):
        lnphi = self.b[:,np.newaxis] / self.bm[np.newaxis,:] * (Z[np.newaxis,:] - 1) - \
        np.log((Z[np.newaxis,:] - B[np.newaxis,:])) - A[np.newaxis,:] / (2 * (2
        ** (1/2)) * B[np.newaxis,:]) * (2 * self.psi / self.aalpha[np.newaxis,:] - \
        self.b[:,np.newaxis] / self.bm[np.newaxis,:]) * np.log((Z[np.newaxis,:] + (1 +
        2 ** (1/2)) * B[np.newaxis,:]) / (Z[np.newaxis,:] + (1 - 2 ** (1/2)) * B[np.newaxis,:]))

        return lnphi


    """ Derivatives - Still need to organize this"""

    def get_all_derivatives(self, fprop):

        x = fprop.xkj[0:ctes.Nc,0,:]
        y = fprop.xkj[0:ctes.Nc,1,:]
        Nl = fprop.Nj[0,0,:]
        Nv = fprop.Nj[0,1,:]
        P = fprop.P
        T = fprop.T

        dlnfildP, dlnfildnij, dZldP_parcial, dZldnij_parcial, Zl = \
                self.get_phase_derivatives(P, T, x, Nl, np.ones(ctes.n_volumes))
        dlnfivdP, dlnfivdnij, dZvdP_parcial, dZvdnij_parcial, Zv = \
                self.get_phase_derivatives(P, T, y, Nv, np.zeros(ctes.n_volumes))
        dnldP, dnvdP, dnldNk, dnvdNk, dnildP, dnivdP, dnildNk, dnivdNk = \
                self.dnij_dNk_dP(dlnfildP, dlnfivdP, dlnfildnij, dlnfivdnij, Nl, Nv)
        dZldP, dZvdP, dZldNk, dZvdNk = self.dZ_dP_dNk(dZldP_parcial,
                dZvdP_parcial, dZldnij_parcial, dZvdnij_parcial, dnildP,
                dnivdP, dnildNk, dnivdNk)
        dVtdP, dVtdNk = self.dVt_dP_dNk(dnldP, dnvdP, dnldNk, dnvdNk, dZldP,
                        dZvdP, dZldNk, dZvdNk, P, T, Zl, Zv, Nl, Nv)
        return dVtdNk, dVtdP

    def get_phase_derivatives(self, P, T, xij, Nj, ph):
        A, B = self.coefficients_cubic_EOS_vectorized(xij, P)
        Z = self.Z_vectorized(A, B, ph)
        dAdP = self.dA_dP()
        dBdP = self.dB_dP()
        dbdnij = self.db_dnij(Nj)
        dadnij = self.da_dnij(Nj, xij)
        dAdnij = self.dA_dnij(P, T, Nj, dadnij)
        dBdnij = self.dB_dnij(P, T, Nj, dbdnij)
        dZdP_parcial = self.dZ_dP_parcial(dAdP, dBdP, Z, A, B)
        dZdnij_parcial = self.dZ_dnij_parcial(dAdnij, dBdnij, Z, A, B)
        dlnphidP = self.dlnphi_dP(dAdP, dBdP, dZdP_parcial, Z, A, B)
        dlnphidnij = self.dlnphi_dnij(dAdnij, dBdnij, dZdnij_parcial, Z, A, B, dbdnij, dadnij, Nj, P)
        dlnfijdP = self.dlnfij_dP(P, dlnphidP)
        dlnfijdnij = self.dlnfij_dnij(dlnphidnij, xij, Nj)
        return dlnfijdP, dlnfijdnij, dZdP_parcial, dZdnij_parcial, Z

    def da_dnij(self, Nj, xkj):
        dadnij = np.empty((ctes.Nc,ctes.n_volumes))
        dadnij[:,Nj!=0] = 2 / Nj[Nj!=0] * ((xkj.T[Nj!=0] @ self.aalpha_ik).T - self.aalpha[Nj!=0])
        dadnij[:,Nj==0] = 0.
        dadnij = dadnij[np.newaxis,:,:]
        return dadnij

    def db_dnij(self, Nj):
        dbdnij = np.empty((ctes.Nc,ctes.n_volumes))
        dbdnij[:,Nj!=0] = (self.b[:,np.newaxis] - self.bm[Nj!=0]) / (Nj[Nj!=0])
        dbdnij[:,Nj==0] = 0.
        dbdnij = dbdnij[np.newaxis,:,:]
        return dbdnij

    def dA_dnij(self, P, T, Nj, dadnij):
        dAdnij = P[np.newaxis,np.newaxis,:] / (ctes.R * T)**2 * dadnij
        return dAdnij

    def dB_dnij(self, P, T, Nj, dbdnij):
        dBdnij = P[np.newaxis,np.newaxis,:] / (ctes.R * T) * dbdnij
        return dBdnij

    def dA_dP(self):
        dAdP = self.aalpha / (ctes.R * self.T) ** 2
        return dAdP

    def dB_dP(self):
        dBdP = self.bm / (ctes.R * self.T)
        return dBdP

    def dZ_dP_parcial(self, dAdP, dBdP, Z, A, B): #OK
        coef0 = -(Z - B)
        coef1 = -(Z ** 2 - (6 * B + 2) * Z - A + 2 * B + 3 * B **2)
        coef2 = 3 * Z ** 2 - 2 * Z * (1 - B) + (A - 3 * B ** 2 - 2 * B)
        dZdP =  (dAdP * coef0 + dBdP * coef1) / coef2
        return dZdP

    def dZ_dnij_parcial(self,dAdnij, dBdnij, Z, A, B):
        coef0 = -(Z - B)
        coef1 = -(Z ** 2 - (6 * B + 2) * Z - A + 2 * B + 3 * B **2)
        coef2 = 3 * Z ** 2 - 2 * Z * (1 - B) + (A - 3 * B ** 2 - 2 * B)
        dZdnij = (dAdnij * coef0 + dBdnij * coef1) / coef2
        return dZdnij

    def dlnphi_dP(self, dAdP, dBdP, dZdP, Z, A, B): #OK i guess
        #ai = np.sum(self.aalpha_ik, axis=1)
        coef0 = self.b[:,np.newaxis] / self.bm #coef8
        coef1 = 1 / (Z - B) #coef5
        coef2 = 2 * self.psi / self.aalpha - coef0 #coef0
        coef3 = np.log((Z + (1 + 2**(1/2))*B)/(Z + (1 - 2**(1/2))*B)) #coef1
        coef4 = 2 * (2 **(1/2)) / (Z**2 + 2*Z*B - B**2) #coef2

        dlnphiijdP = coef0 * dZdP - coef1 * (dZdP - dBdP) - coef2 / (2*(2 **(1/2))) * ((B * dAdP - A * dBdP) /  (B **2) * coef3 +
                    A / B *  coef4 * (Z * dBdP - B * dZdP))
        return dlnphiijdP

    def dlnphi_dnij(self, dAdnij, dBdnij, dZdnij, Z, A, B, dbdnij, dadnij, Nj, P):
        coef0 = self.b[:,np.newaxis] / self.bm #coef8
        coef1 = 1 / (Z - B) #coef5
        coef2 = 2 * self.psi / self.aalpha - coef0 #coef0
        coef3 = np.log((Z + (1 + 2**(1/2))*B)/(Z + (1 - 2**(1/2))*B)) #coef1
        coef4 = 2 * (2 **(1/2)) / (Z**2 + 2*Z*B - B**2) #coef2

        dlnphiijdnij = np.empty((ctes.Nc, ctes.Nc, len(P)))
        dlnphiijdnij[:,:,Nj==0] = 0

        aux1 = coef0[:,np.newaxis,:] / self.bm * (self.bm * dZdnij - (Z - 1) * dbdnij) - coef1 * (dZdnij - dBdnij)
        aux2_1 = coef2[:,np.newaxis,:] * ((B * dAdnij - A * dBdnij) / B**2 * coef3 + A / B * coef4 * (Z * dBdnij - B * dZdnij))
        aux2_2_1 =  2 * (self.aalpha_ik[:,:,np.newaxis] - (Nj / self.aalpha * dadnij + 1) * self.psi[:,np.newaxis,:])
        aux2_2_2 = coef0[:,np.newaxis,:] / self.bm * dbdnij
        dlnphiijdnij[:,:,Nj!=0] = aux1[:,:,Nj!=0] - 1 / (2 * (2 **(1/2))) * (aux2_1[:,:,Nj!=0] + A[Nj!=0] / B[Nj!=0] * coef3[Nj!=0] *
                                (aux2_2_1 [:,:,Nj!=0] / (Nj[Nj!=0] * self.aalpha[Nj!=0]) + aux2_2_2[:,:,Nj!=0]))
        return dlnphiijdnij

    def dlnfij_dP(self, P, dlnphidP):
        dlnfijdP = dlnphidP + 1 / P[np.newaxis,:]
        return dlnfijdP

    def dlnfij_dnij(self, dlnphidnij, xij, Nj):
        dlnfijdnij = np.empty(dlnphidnij.shape)
        dlnfijdnij[:,:,Nj!=0] = (dlnphidnij[:,:,Nj!=0] + 1 / xij[:,np.newaxis,Nj!=0] * (np.identity(ctes.Nc)[:,:,np.newaxis] -
                                xij[:,np.newaxis, Nj!=0]) / Nj[Nj!=0])

        dlnfijdnij[:,:,Nj==0] = 0
        return dlnfijdnij

    def dnij_dNk_dP(self, dlnfildP, dlnfivdP, dlnfildnij, dlnfivdnij, Nl, Nv):
        dnivdNk = np.empty(dlnfildnij.shape)
        dnivdP = np.empty(dlnfildP.shape)
        dnildNk = np.empty(dlnfildnij.shape)
        dnildP = np.empty(dlnfildP.shape)

        aux = np.ones(len(Nv),dtype=bool)
        aux[Nv == 0] = False
        aux[Nl == 0] = False

        '''Making dlnf's matrix'''
        dlnfildnij = np.transpose(dlnfildnij[:,:,aux], (2,0,1))
        dlnfivdnij = np.transpose(dlnfivdnij[:,:,aux], (2,0,1))
        matrix = dlnfildnij + dlnfivdnij

        '''Independent tem vector for dnsdP '''
        v = dlnfildP[:,aux] - dlnfivdP[:,aux]

        '''Solving system of equations for dnivdNk and dnivdP'''
        dnivdP_aux = dnivdP[:,aux]
        dnivdNk_aux = dnivdNk[:,:,aux]
        dnivdNk_aux = np.linalg.solve(matrix, dlnfildnij).transpose((1,2,0))
        dnivdP_aux =  np.linalg.solve(matrix, v.T).T


        dnivdNk[:,:,aux] = dnivdNk_aux
        dnivdNk[:,:,Nv==0] = 0
        dnivdNk[:,:,Nl==0] = np.identity(ctes.Nc)[:,:,np.newaxis]

        dnivdP[:,aux] = dnivdP_aux
        dnivdP[:,Nv==0] = 0
        dnivdP[:,Nl==0] = 0

        '''Calculate dnildP and dnildNk'''
        dnildP = - dnivdP
        dnildNk = np.identity(ctes.Nc)[:,:,np.newaxis] - dnivdNk

        dnvdP = np.sum(dnivdP, axis = 0)
        dnldP = - dnvdP
        dnldNk = np.sum(dnildNk, axis = 0)
        dnvdNk = np.sum(dnivdNk, axis = 0)
        return dnldP, dnvdP, dnldNk, dnvdNk, dnildP, dnivdP, dnildNk, dnivdNk

    def dZ_dP_dNk(self, dZldP, dZvdP, dZldnij, dZvdnij, dnildP, dnivdP, dnildNk, dnivdNk):
        dZldP = dZldP + np.sum(dZldnij.sum(axis=0) * dnildP, axis = 0)
        dZvdP = dZvdP + np.sum(dZvdnij.sum(axis=0) * dnivdP, axis = 0)
        dZldNk = np.sum(dZldnij[0][:,np.newaxis,:] * dnildNk, axis = 0)
        dZvdNk = np.sum(dZvdnij[0][:,np.newaxis,:] * dnivdNk, axis = 0)
        return dZldP, dZvdP, dZldNk, dZvdNk

    def dVt_dP_dNk(self, dnldP, dnvdP, dnldNk, dnvdNk, dZldP, dZvdP, dZldNk, dZvdNk, P, T, Zl, Zv, Nl, Nv):
        coef = ctes.R * T / P
        dVtdP = coef * (Nl * dZldP + Nv * dZvdP + Zl * dnldP + Zv * dnvdP - (Zl * Nl / P + Zv * Nv / P))
        dVtdNk = coef * (Nl * dZldNk + Nv * dZvdNk + Zl * dnldNk + Zv * dnvdNk)
        return dVtdP, dVtdNk
