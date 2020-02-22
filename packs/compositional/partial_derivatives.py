import numpy as np
from .stability_check import StabilityCheck

class PartialDerivatives:

    def __init__(self):
        self.n_phases = 2

    def d_dn_all_derivatives(self, fp, Nphase, nkphase, kprop, l):
        delta = 0.0001
        dlnf_dn = np.zeros([fp.Nc, fp.Nc, self.n_phases])
        dZ_dn = np.zeros([fp.Nc, self.n_phases])
        for ph in range(self.n_phases):
            if Nphase[ph] != 0:
                for i in range(0,fp.Nc):
                    l_plus = np.copy(l[:,ph]); l_minus = np.copy(l[:,ph])
                    l_plus[i] = (nkphase[i,ph] + delta / 2) / Nphase[ph]
                    l_minus[i] = (nkphase[i,ph] - delta / 2) / Nphase[ph]
                    dlnf_dn[:,i,ph] = (PartialDerivatives.lnf(fp, kprop, l_plus, ph)
                     - PartialDerivatives.lnf(fp, kprop, l_minus, ph)) / delta
                    dZ_dn[i,ph] = (self.Z(fp, kprop, l_plus, ph)
                     - self.Z(fp, kprop, l_minus, ph)) / delta
        return dlnf_dn, dZ_dn

    def d_dP_all_derivatives(self, fp, Nphase, kprop, l, b):
        delta = 0.0001
        dZ_dP = np.zeros(self.n_phases)
        dlnf_dP = np.zeros([fp.Nc, self.n_phases])
        for ph in range(0, self.n_phases):
            if Nphase[ph] != 0:
                fp.P = self.Pvolume + delta/2
                Z_plus = self.Z(fp, kprop, l[:,ph], ph)
                lnf_plus = PartialDerivatives.lnf(fp, kprop, l[:,ph], ph)
                fp.P = self.Pvolume - delta / 2
                dZ_dP[ph] = (Z_plus - self.Z(fp, kprop, l[:,ph], ph)) / delta
                dlnf_dP[:,ph] = (lnf_plus - PartialDerivatives.lnf(fp, kprop, l[:,ph], ph)) / delta
        return dlnf_dP, dZ_dP

    def lnf(fp, kprop, l, ph):
        lnf = fp.lnphi(kprop, l, ph) + np.log(fp.P * l)
        return lnf

    def Z(self, fp, kprop, l, ph):
        A, B = fp.coefficientsPR(kprop, l)
        return StabilityCheck.Z_PR(B, A, ph)

    def dVt_derivatives(self, fluid_properties, fp, kprop):
        Nphase_allvolumes = fluid_properties.mole_numbers_o_and_g
        nkphase_allvolumes = fluid_properties.component_phase_mole_numbers
        l_allvolumes = fluid_properties.component_molar_fractions
        ksi_j = fluid_properties.ksi_o_and_g

        n_vols = len(Nphase_allvolumes[0,0,:])

        """ Initializing dN derivatives """
        dnij_dNk = np.zeros([fp.Nc, fp.Nc, self.n_phases, n_vols])
        dZj_dnij = np.zeros([fp.Nc, self.n_phases, n_vols])

        """ Initializing dP derivatives """
        dnij_dP = np.zeros([fluid_properties.Nc, self.n_phases, n_vols])
        dZj_dP = np.zeros([1, self.n_phases, n_vols])
        P = np.copy(fluid_properties.P)
        Zj = np.zeros([1,self.n_phases, n_vols])
        for b in range(n_vols):
            self.Pvolume = P[b]
            fp.P = self.Pvolume
            self.y = l_allvolumes[:,0,b]
            self.x = l_allvolumes[:,1,b]
            l = np.zeros([fp.Nc, self.n_phases])
            l[:,0] = self.y[0:fp.Nc]; l[:,1] = self.x[0:fp.Nc]
            Nphase = Nphase_allvolumes[0,:,b]
            nkphase = nkphase_allvolumes[:,:,b]
            dlnfij_dnij, dZj_dnij[:,:,b] = self.d_dn_all_derivatives(fp, Nphase, nkphase, kprop, l)
            dlnfj_dP, dZj_dP[0,:,b] = self.d_dP_all_derivatives(fp, Nphase, kprop, l, b)
            Zj[0,:,b] = np.array([self.Z(fp, kprop, self.y, 0), self.Z(fp, kprop, self.x, 1)])
            matrix = np.sum(dlnfij_dnij, axis = 2)

            dlnf_dP_vector = dlnfj_dP[:,0] - dlnfj_dP[:,1]
            dnij_dP[:,1,b] =  np.linalg.inv(matrix)@dlnf_dP_vector
            dnij_dP[:,0,b] = - dnij_dP[:,1,b]

            for k in range(0, fp.Nc):
                dlnfl_dnk = dlnfij_dnij[:,k,1]
                dnij_dNk[:,k,1,b] = np.linalg.inv(matrix)@dlnfl_dnk
                dnij_dNk[:,k,0,b] = - dnij_dNk[:,k,1,b]
                dnij_dNk[k,k,0,b] = 1 - dnij_dNk[k,k,1,b]

        dvj_dnij = fluid_properties.R * fluid_properties.T / fluid_properties.P * dZj_dnij
        dVj_dNk = np.sum(dnij_dNk * (1 / ksi_j + Nphase[:,np.newaxis] * dvj_dnij), axis = 0)
        dVt_dNk = np.sum(dVj_dNk, axis = 1)

        dvj_dP = fluid_properties.R * fluid_properties.T / fluid_properties.P * \
                (dZj_dP - Zj / fluid_properties.P)
        dVj_dP = np.sum(Nphase_allvolumes * dvj_dP + np.sum(dnij_dP * (1 / ksi_j +
                Nphase[:,np.newaxis] * dvj_dnij), axis = 0),axis = 0)
        dVt_dP = np.sum(dVj_dP,axis = 0)

        #fluid_properties.P = P #comming back
        return dVt_dNk, dVt_dP
