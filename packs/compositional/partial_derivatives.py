import numpy as np
from .stability_check import StabilityCheck

class PartialDerivatives:

    def __init__(self):
        self.n_phases = 2

    def d_dn_all_derivatives(self, fprop_block, Nphase, nkphase, kprop, l):
        delta = 0.0001
        dlnf_dn = np.zeros([fprop_block.Nc, fprop_block.Nc, self.n_phases])
        dZ_dn = np.zeros([fprop_block.Nc, self.n_phases])
        for ph in range(self.n_phases):
            if Nphase[ph] != 0:
                for i in range(0,fprop_block.Nc):
                    l_plus = np.copy(l[:,ph]); l_minus = np.copy(l[:,ph])
                    l_plus[i] = (nkphase[i,ph] + delta / 2) / Nphase[ph]
                    l_minus[i] = (nkphase[i,ph] - delta / 2) / Nphase[ph]
                    dlnf_dn[:,i,ph] = (PartialDerivatives.lnf(fprop_block, kprop, l_plus, ph)
                     - PartialDerivatives.lnf(fprop_block, kprop, l_minus, ph)) / delta
                    dZ_dn[i,ph] = (self.Z(fprop_block, kprop, l_plus, ph)
                     - self.Z(fprop_block, kprop, l_minus, ph)) / delta
        return dlnf_dn, dZ_dn

    def d_dP_all_derivatives(self, fprop_block, Nphase, kprop, l, b):
        delta = 0.0001
        dZ_dP = np.zeros(self.n_phases)
        dlnf_dP = np.zeros([fprop_block.Nc, self.n_phases])
        for ph in range(0, self.n_phases):
            if Nphase[ph] != 0:
                fprop_block.P = self.Pvolume + delta/2
                Z_plus = self.Z(fprop_block, kprop, l[:,ph], ph)
                lnf_plus = PartialDerivatives.lnf(fprop_block, kprop, l[:,ph], ph)
                fprop_block.P = self.Pvolume - delta / 2
                dZ_dP[ph] = (Z_plus - self.Z(fprop_block, kprop, l[:,ph], ph)) / delta
                dlnf_dP[:,ph] = (lnf_plus - PartialDerivatives.lnf(fprop_block, kprop, l[:,ph], ph)) / delta
        return dlnf_dP, dZ_dP

    def lnf(fprop_block, kprop, l, ph):
        lnf = fprop_block.lnphi(kprop, l, ph) + np.log(fprop_block.P * l)
        return lnf

    def Z(self, fprop_block, kprop, l, ph):
        A, B = fprop_block.coefficientsPR(kprop, l)
        return StabilityCheck.Z_PR(B, A, ph)

    def dVt_derivatives(self, fprop, fprop_block, kprop):
        Nphase_allvolumes = fprop.mole_numbers_o_and_g
        nkphase_allvolumes = fprop.component_phase_mole_numbers
        l_allvolumes = fprop.component_molar_fractions
        ksi_j = fprop.ksi_o_and_g

        n_vols = len(Nphase_allvolumes[0,0,:])

        """ Initializing dN derivatives """
        dnij_dNk = np.zeros([fprop_block.Nc, fprop_block.Nc, self.n_phases, n_vols])
        dZj_dnij = np.zeros([fprop_block.Nc, self.n_phases, n_vols])

        """ Initializing dP derivatives """
        dnij_dP = np.zeros([fprop.Nc, self.n_phases, n_vols])
        dZj_dP = np.zeros([1, self.n_phases, n_vols])
        P = np.copy(fprop.P)
        Zj = np.zeros([1,self.n_phases, n_vols])
        for b in range(n_vols):
            self.Pvolume = P[b]
            fprop_block.P = self.Pvolume
            self.y = l_allvolumes[:,0,b]
            self.x = l_allvolumes[:,1,b]
            l = np.zeros([fprop_block.Nc, self.n_phases])
            l[:,0] = self.y[0:fprop_block.Nc]; l[:,1] = self.x[0:fprop_block.Nc]
            Nphase = Nphase_allvolumes[0,:,b]
            nkphase = nkphase_allvolumes[:,:,b]
            dlnfij_dnij, dZj_dnij[:,:,b] = self.d_dn_all_derivatives(fprop_block, Nphase, nkphase, kprop, l)
            dlnfj_dP, dZj_dP[0,:,b] = self.d_dP_all_derivatives(fprop_block, Nphase, kprop, l, b)
            Zj[0,:,b] = np.array([self.Z(fprop_block, kprop, self.y, 0), self.Z(fprop_block, kprop, self.x, 1)])
            matrix = np.sum(dlnfij_dnij, axis = 2)

            dlnf_dP_vector = dlnfj_dP[:,0] - dlnfj_dP[:,1]
            dnij_dP[:,1,b] =  np.linalg.inv(matrix)@dlnf_dP_vector
            dnij_dP[:,0,b] = - dnij_dP[:,1,b]

            for k in range(0, fprop_block.Nc):
                dlnfl_dnk = dlnfij_dnij[:,k,1]
                dnij_dNk[:,k,1,b] = np.linalg.inv(matrix)@dlnfl_dnk
                dnij_dNk[:,k,0,b] = - dnij_dNk[:,k,1,b]
                dnij_dNk[k,k,0,b] = 1 - dnij_dNk[k,k,1,b]

        dvj_dnij = fprop.R * fprop.T / fprop.P * dZj_dnij
        dVj_dNk = np.sum(dnij_dNk * (1 / ksi_j + Nphase[:,np.newaxis] * dvj_dnij), axis = 0)
        dVt_dNk = np.sum(dVj_dNk, axis = 1)

        dvj_dP = fprop.R * fprop.T / fprop.P * \
                (dZj_dP - Zj / fprop.P)
        dVj_dP = np.sum(Nphase_allvolumes * dvj_dP + np.sum(dnij_dP * (1 / ksi_j +
                Nphase[:,np.newaxis] * dvj_dnij), axis = 0),axis = 0)
        dVt_dP = np.sum(dVj_dP,axis = 0)

        #fprop.P = P #comming back
        return dVt_dNk, dVt_dP
