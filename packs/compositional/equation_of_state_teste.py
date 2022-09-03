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
        l_reshape = np.ones_like(self.aalpha_ik)[:,:,np.newaxis] * l[:,np.newaxis,:]
        self.aalpha = (l_reshape * l[np.newaxis,:,:] * self.aalpha_ik[:,:,np.newaxis]).sum(axis=0).sum(axis=0)
        oper = P / (ctes.R* self.T)
        B = self.bm * oper
        A = self.aalpha * oper / (ctes.R* self.T)
        self.psi = (l_reshape * self.aalpha_ik[:,:,np.newaxis]).sum(axis = 0)
        return A, B

    def Z_vectorized(self, A, B):
        coef = np.empty([4,len(B.ravel())])
        B_square = B * B

        coef[0,:] = 1
        coef[1,:] = -(1 - B)
        coef[2,:] = (A - 2*B - 3*B_square)
        coef[3,:] = -(A*B - B_square - B_square * B)
        Z = CubicRoots().run(coef)
        root = np.isreal(Z)
        n_reais = np.sum(root.astype(np.int), axis=1)

        'if n_reais == 2'
        aux_reais = (n_reais==2)
        Z_reais_n2 = np.reshape(Z[aux_reais][root[aux_reais]],(len(aux_reais[aux_reais]),2))
        Z_aux_reais = Z[aux_reais]
        Z_aux_reais[~root[aux_reais]] = np.max(Z_reais_n2,axis=1, initial=0)
        Z[aux_reais] = Z_aux_reais

        'if n_reais==1'

        aux_reais = (n_reais==1)
        Z_aux_reais = Z[aux_reais]
        Z_aux_reais[~root[aux_reais]] = np.repeat(Z[aux_reais][root[aux_reais]], 2)
        Z[aux_reais] = Z_aux_reais

        'if any Z<0'

        aux_neg = np.zeros(Z.shape,dtype=bool)
        aux_neg[Z<B[:,np.newaxis]] = True
        aux_neg_sum = aux_neg.sum(axis=1)
        lines_Zneg1 = (aux_neg_sum == 1).astype(bool)
        lines_Zneg2 = (aux_neg_sum == 2).astype(bool)
        Zneg1 = Z[lines_Zneg1]
        Zneg1[Zneg1<B[lines_Zneg1,np.newaxis]] = np.max(Z[lines_Zneg1],axis=1, initial=0)
        Z[lines_Zneg1] = Zneg1
        Zneg2 = Z[lines_Zneg2]
        Zneg2[Zneg2<B[lines_Zneg2,np.newaxis]] = np.repeat(np.max(Z[lines_Zneg2],axis=1, initial=0),2)
        Z[lines_Zneg2] = Zneg2
        return Z

    def getZ_by_phase(self, Z, ph):
        Z = np.min(Z, axis = 1) * ph + np.max(Z, axis = 1) * (1 - ph)
        Z = np.real(Z)
        return Z

    def lnphi_Z_deltaG(self, xkj, P, ph):
        A, B = self.coefficients_cubic_EOS_vectorized(xkj, P)
        Z_all = self.Z_vectorized(A, B)
        Z_ph1 = self.getZ_by_phase(Z_all, ph)
        Z_ph2 = self.getZ_by_phase(Z_all, 1-ph)
        Z = np.concatenate((Z_ph1[:,np.newaxis],Z_ph2[:,np.newaxis]),axis=-1)
        lnphi = self.lnphi_calculation_deltaG(A, B, Z)
        return lnphi, Z

    def lnphi(self, xkj, P, ph):
        A, B = self.coefficients_cubic_EOS_vectorized(xkj, P)
        Z_all = self.Z_vectorized(A, B)
        Z = self.getZ_by_phase(Z_all, ph)
        lnphi = self.lnphi_calculation(A, B, Z)
        #if any(np.isnan(lnphi).ravel()): import pdb; pdb.set_trace()
        return lnphi

    def lnphi_calculation(self, A, B, Z):
        b_bm = self.b[:,np.newaxis] / self.bm[np.newaxis,:]
        lnphi = b_bm * (Z[np.newaxis,:] - 1) - np.log((Z - B))[np.newaxis,:] - \
                (A / (2 * (2 ** (1/2)) * B))[np.newaxis,:] \
                * (2 * self.psi / self.aalpha[np.newaxis,:] - b_bm) * \
                np.log((Z + (1 + 2 ** (1/2)) * B) /
                (Z + (1 - 2 ** (1/2)) * B))[np.newaxis,:]
        return lnphi

    def lnphi_calculation_deltaG(self, A, B, Z):
        b_bm = self.b[:,np.newaxis] / self.bm[np.newaxis,:]
        lnphi = b_bm[...,np.newaxis] * (Z[np.newaxis,:] - 1) - np.log((Z - B[:,np.newaxis])[np.newaxis,:]) - \
                ((A / (2 * (2 ** (1/2)) * B))[np.newaxis,:] \
                * (2 * self.psi / self.aalpha[np.newaxis,:] - b_bm))[...,np.newaxis] * \
                np.log((Z + (1 + 2 ** (1/2)) * B[:,np.newaxis]) /
                (Z + (1 - 2 ** (1/2)) * B[:,np.newaxis]))[np.newaxis,:]
        return lnphi

    """ Derivatives - Still need to organize this"""

    def get_all_derivatives(self, fprop):

        x = fprop.xkj[:,0,:]
        y = fprop.xkj[:,1,:]
        w_ = fprop.xkj[:,2,:]

        Nl = fprop.Nj[0,0,:]
        Nv = fprop.Nj[0,1,:]
        Nw = fprop.Nj[0,2,:]

        P = fprop.P
        T = fprop.T

        dlnfildP, dlnfildnij, dZldP_parcial, dZldnij_parcial, Zl = \
                self.get_phase_derivatives(P, T, x, Nl, np.ones(ctes.n_volumes))
        dlnfivdP, dlnfivdnij, dZvdP_parcial, dZvdnij_parcial, Zv = \
                self.get_phase_derivatives(P, T, y, Nv, np.zeros(ctes.n_volumes))
        dlnfiwdP, dlnfiwdnij, dZwdP_parcial, dZwdnij_parcial, Zw = \
                self.get_phase_derivatives(P, T, w_, Nw, np.ones(ctes.n_volumes))
        #import pdb; pdb.set_trace()

        dnldP, dnvdP, dnwdP, dnldNk, dnvdNk, dnwdNk, dnildP, dnivdP, dniwdP, dnildNk, dnivdNk, dniwdNk = \
                self.dnij_dNk_dP(dlnfildP, dlnfivdP, dlnfiwdP, dlnfildnij, dlnfivdnij, dlnfiwdnij, Nl, Nv, Nw)
        #import pdb; pdb.set_trace()
        dZldP, dZvdP, dZwdP, dZldNk, dZvdNk, dZwdNk = self.dZ_dP_dNk(dZldP_parcial,
               dZvdP_parcial, dZwdP_parcial, dZldnij_parcial, dZvdnij_parcial, dZwdnij_parcial,
               dnildP, dnivdP, dniwdP, dnildNk, dnivdNk, dnivdNk)
        #import pdb; pdb.set_trace()
        dVldP, dVvdP, dVwdP, dVldNk, dVvdNk, dVwdNk = self.dVj_dP_dNk(dnldP, dnvdP, dnwdP,
                                                      dnldNk, dnvdNk, dnwdNk, dZldP, dZvdP, dZwdP,
                                                      dZldNk, dZvdNk, dZwdNk, P, T, Zl, Zv, Zw, Nl, Nv, Nw)

        return dVldP, dVvdP, dVwdP, dVldNk, dVvdNk, dVwdNk

    def get_phase_derivatives(self, P, T, xij, Nj, ph):
        A, B = self.coefficients_cubic_EOS_vectorized(xij, P)
        Z = self.Z_vectorized(A, B)
        Z = self.getZ_by_phase(Z, ph)
        xij[xij==0] = 1e-300 #OBSERVAR SE ISSO TA FUNCIONANDO
        self.Z_square = Z * Z
        self.B_square = B * B

        dAdP = self.dA_dP(T)
        dBdP = self.dB_dP(T)
        dbdnij = self.db_dnij(Nj)
        dadnij = self.da_dnij(Nj, xij)
        dAdnij = self.dA_dnij(P, T, Nj, dadnij)
        dBdnij = self.dB_dnij(P, T, Nj, dbdnij)
        dZdP_parcial, dZdnij_parcial = self.dZ_dP_dnij_parcial(dAdP, dBdP, dAdnij, dBdnij, Z, A, B)
        dlnphidP, dlnphidnij = self.dlnphi_dP_dnij(dAdP, dBdP, dZdP_parcial,dAdnij, dBdnij, dZdnij_parcial, Z, A, B, dbdnij, dadnij, Nj, P)
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

    def dA_dP(self, T):
        dAdP = self.aalpha / (ctes.R * T) ** 2
        return dAdP

    def dB_dP(self, T):
        dBdP = self.bm / (ctes.R * T)
        return dBdP

    def dZ_dP_dnij_parcial(self, dAdP, dBdP, dAdnij, dBdnij, Z, A, B): #OK

        'coefs'
        coef0 = -(Z - B)
        coef1 = -(self.Z_square - (6 * B + 2) * Z - A + 2 * B + 3 * self.B_square)
        coef2 = 3 * self.Z_square - 2 * Z * (1 - B) + (A - 3 * self.B_square - 2 * B)

        'dZ_dP'
        dZdP =  (dAdP * coef0 + dBdP * coef1) / coef2

        'dZ_dnij'
        dZdnij = (dAdnij * coef0 + dBdnij * coef1) / coef2
        return dZdP, dZdnij

    def dlnphi_dP_dnij(self, dAdP, dBdP, dZdP,dAdnij, dBdnij, dZdnij, Z, A, B, dbdnij, dadnij, Nj, P): #OK i guess
        #ai = np.sum(self.aalpha_ik, axis=1)
        coef0 = self.b[:,np.newaxis] / self.bm #coef8
        coef1 = 1 / (Z - B) #coef5
        coef2 = 2 * self.psi / self.aalpha - coef0 #coef0
        coef3 = np.log((Z + (1 + 2**(1/2))*B)/(Z + (1 - 2**(1/2))*B)) #coef1
        coef4 = 2 * (2 **(1/2)) / (self.Z_square + 2*Z*B - self.B_square) #coef2

        dlnphiijdP = coef0 * dZdP - coef1 * (dZdP - dBdP) - coef2 / (2*(2 **(1/2))) \
                    * ((B * dAdP - A * dBdP) /  (self.B_square) * coef3 +
                    A / B *  coef4 * (Z * dBdP - B * dZdP))

        dlnphiijdnij = np.empty((ctes.Nc, ctes.Nc, len(P)))
        dlnphiijdnij[:,:,Nj==0] = 0

        aux1 = coef0[:,np.newaxis,:] / self.bm * (self.bm * dZdnij - (Z - 1) \
                * dbdnij) - coef1 * (dZdnij - dBdnij)
        aux2_1 = coef2[:,np.newaxis,:] * ((B * dAdnij - A * dBdnij) / self.B_square * coef3 \
                + A / B * coef4 * (Z * dBdnij - B * dZdnij))
        aux2_2_1 =  2 * (self.aalpha_ik[:,:,np.newaxis] - (Nj / self.aalpha * dadnij + 1) \
                * self.psi[:,np.newaxis,:])
        aux2_2_2 = coef0[:,np.newaxis,:] / self.bm * dbdnij
        dlnphiijdnij[:,:,Nj!=0] = aux1[:,:,Nj!=0] - 1 / (2 * (2 **(1/2))) * \
                (aux2_1[:,:,Nj!=0] + A[Nj!=0] / B[Nj!=0] * coef3[Nj!=0] *
                (aux2_2_1 [:,:,Nj!=0] / (Nj[Nj!=0] * self.aalpha[Nj!=0]) + \
                aux2_2_2[:,:,Nj!=0]))
        return dlnphiijdP, dlnphiijdnij

    def dlnfij_dP(self, P, dlnphidP):
        dlnfijdP = dlnphidP + 1 / P[np.newaxis,:]
        return dlnfijdP

    def dlnfij_dnij(self, dlnphidnij, xij, Nj):
        dlnfijdnij = np.empty(dlnphidnij.shape)
        dlnfijdnij[:,:,Nj!=0] = (dlnphidnij[:,:,Nj!=0] + 1 / xij[:,np.newaxis,Nj!=0] * \
                                (np.identity(ctes.Nc)[:,:,np.newaxis] -
                                xij[:,np.newaxis, Nj!=0]) / Nj[Nj!=0])

        #dlnfijdnij[np.isnan(dlnfijdnij)] = 0
        #dlnfijdnij[np.isinf(dlnfijdnij)] = 0
        dlnfijdnij[:,:,Nj==0] = 0
        return dlnfijdnij

    def dnij_dNk_dP(self, dlnfildP, dlnfivdP, dlnfiwdP, dlnfildnij, dlnfivdnij, dlnfiwdnij, Nl, Nv, Nw):
        dnivdNk = np.empty(dlnfildnij.shape)
        dnivdP = np.empty(dlnfildP.shape)
        dnildNk = np.empty(dlnfildnij.shape)
        dnildP = np.empty(dlnfildP.shape)
        dniwdNk = np.empty(dlnfildnij.shape)
        dniwdP = np.empty(dlnfildP.shape)

        aux = np.ones(len(Nv),dtype=bool)
        aux[Nv == 0] = False
        aux[Nl == 0] = False
        aux[Nw == 0] = False

        '''Making dlnf's matrix'''
        dlnfildnij = np.transpose(dlnfildnij[:,:,aux], (2,0,1))
        dlnfivdnij = np.transpose(dlnfivdnij[:,:,aux], (2,0,1))
        dlnfiwdnij = np.transpose(dlnfiwdnij[:,:,aux], (2,0,1))
        A11 = dlnfildnij + dlnfivdnij
        A12 = dlnfivdnij
        A21 = dlnfivdnij
        A22 = dlnfiwdnij + dlnfivdnij

        '''Independent tem vector for dnsdP '''
        B1 = dlnfivdP[:,aux] - dlnfildP[:,aux]
        B2 = dlnfivdP[:,aux] - dlnfiwdP[:,aux]

        ''' Montagem do sistema three-phase'''
        A_3ph = np.empty((1, 2*len(B1), 2*len(B1)))
        B_3ph = np.empty((2*len(B1),1))

        B_3ph[0:len(B1),:] = B1
        #import pdb; pdb.set_trace()

        B_3ph[(len(B1)):(2*len(B1)),:] = B2
        #ver = np.shape(B1)
        #import pdb; pdb.set_trace()
        A_3ph[:, 0:len(B1), 0:len(B1)] = A11
        A_3ph[:,len(B1):2*len(B1),0:len(B1)] = A21
        A_3ph[:,0:len(B1),len(B1):2*len(B1)] = A12
        A_3ph[:,len(B1):2*len(B1),len(B1):2*len(B1)] = A22

        '''Solving system of equations for dnijdP'''
        aux_3ph = np.ones(2*len(Nv),dtype=bool)
        aux_3ph[0:len(Nv),] = aux
        aux_3ph[len(Nv):2*len(Nv),] = aux
        dnilwdP = np.empty((len(B_3ph), len(aux_3ph)))
        dnilwdP_aux = dnilwdP[:,aux_3ph]

        dnilwdP_aux =  np.linalg.solve(A_3ph, B_3ph.T).T

        '''Calculate dnijdP '''
        dnildP[:,aux] = dnilwdP_aux[0:len(B1),:]
        dnildP[:,Nv==0] = 0
        dnildP[:,Nl==0] = 0
        dnildP[:,Nw==0] = 0

        dniwdP[:,aux] = dnilwdP_aux[len(B1):2*len(B1),:]
        dniwdP[:,Nv==0] = 0
        dniwdP[:,Nl==0] = 0
        dniwdP[:,Nw==0] = 0

        dnivdP = - dnildP - dniwdP

        '''Solving system of equations for dnilwdNk'''
        dnilwdNk = np.empty((1,len(B_3ph), len(B_3ph)))
        dlnfivdnij_3ph = np.empty((1, 2*len(B1), 2*len(B1)))

        dlnfivdnij_3ph[:,0:len(B1),0:len(B1)] = dlnfivdnij
        dlnfivdnij_3ph[:,len(B1):2*len(B1),len(B1):2*len(B1)] = dlnfivdnij
        try:
            dnilwdNk_aux = np.linalg.solve(A_3ph, dlnfivdnij_3ph).transpose((1,2,0))
        except: import pdb; pdb.set_trace()
        dnildNk[:,:,aux] = dnilwdNk_aux[0:len(B1),0:len(B1),:]
        dniwdNk[:,:,aux] = dnilwdNk_aux[len(B1):2*len(B1),len(B1):2*len(B1),:]
        dnildNk[:,:,Nl==0] = 0
        #dnildNk[:,:,Nw==0] = 0
        #dniwdNk[:,:,Nl==0] = 0
        dniwdNk[:,:,Nw==0] = 0
        dnildNk[:,:,Nv==0] = np.identity(ctes.Nc)[:,:,np.newaxis]
        dnildNk[:,:,Nw==0] = np.identity(ctes.Nc)[:,:,np.newaxis]
        dniwdNk[:,:,Nv==0] = np.identity(ctes.Nc)[:,:,np.newaxis]
        dniwdNk[:,:,Nl==0] = np.identity(ctes.Nc)[:,:,np.newaxis]

        dnivdNk = np.identity(ctes.Nc)[:,:,np.newaxis] - dnildNk - dniwdNk

        '''Calculate dnildP and dnildNk'''
        dnldP = np.sum(dnildP, axis = 0)
        dnwdP = np.sum(dniwdP, axis = 0)
        dnvdP = - dnldP - dnwdP
        dnldNk = np.sum(dnildNk, axis = 0)
        dnvdNk = np.sum(dnivdNk, axis = 0)
        dnwdNk = np.sum(dniwdNk, axis = 0)
        return dnldP, dnvdP,  dnwdP, dnldNk, dnvdNk, dnwdNk, dnildP, dnivdP, dniwdP, dnildNk, dnivdNk, dniwdNk


    def dZ_dP_dNk(self, dZldP, dZvdP, dZwdP, dZldnij, dZvdnij, dZwdnij,
                        dnildP, dnivdP, dniwdP, dnildNk, dnivdNk, dniwdNk):
        dZldP = dZldP + np.sum(dZldnij.sum(axis=0) * dnildP, axis = 0)
        dZvdP = dZvdP + np.sum(dZvdnij.sum(axis=0) * dnivdP, axis = 0)
        dZwdP = dZwdP + np.sum(dZwdnij.sum(axis=0) * dniwdP, axis = 0)
        dZldNk = np.sum(dZldnij[0][:,np.newaxis,:] * dnildNk, axis = 0)
        dZvdNk = np.sum(dZvdnij[0][:,np.newaxis,:] * dnivdNk, axis = 0)
        dZwdNk = np.sum(dZwdnij[0][:,np.newaxis,:] * dniwdNk, axis = 0)
        return dZldP, dZvdP, dZwdP, dZldNk, dZvdNk, dZwdNk

    def dVj_dP_dNk(self, dnldP, dnvdP, dnwdP, dnldNk, dnvdNk, dnwdNk,
                   dZldP, dZvdP, dZwdP, dZldNk, dZvdNk, dZwdNk, P, T, Zl, Zv, Zw, Nl, Nv, Nw):
        coef = ctes.R * T / P
        dVldP = coef * (Nl * dZldP + Zl * dnldP - (Zl * Nl / P))
        dVvdP = coef * (Nv * dZvdP + Zv * dnvdP - (Zv * Nv / P))
        dVwdP = coef * (Nw * dZwdP + Zw * dnwdP - (Zw * Nw / P))
        #dVvdP = coef * (Nv * dZvdP + Zv * dnvdP - (Zv * Nv / P))
        dVldNk = coef * (Nl * dZldNk + Zl * dnldNk)
        dVvdNk = coef * (Nv * dZvdNk + Zv * dnvdNk)
        dVwdNk = coef * (Nw * dZwdNk + Zw * dnwdNk)
        return dVldP, dVvdP, dVwdP, dVldNk, dVvdNk, dVwdNk
