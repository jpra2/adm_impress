import numpy as np
from .stability_check import StabilityCheck

class PartialDerivatives:

    def __init__(self):
        self.n_phases = 2

    def d_dn_all_derivatives(self, fluid_properties, Nphase, l):
        nkphase = l * Nphase
        delta = 0.0001
        dlnf_dn = np.zeros([fluid_properties.Nc, fluid_properties.Nc, self.n_phases])
        dZ_dn = np.zeros([fluid_properties.Nc, self.n_phases])
        for ph in range(self.n_phases):
            if Nphase[ph] != 0:
                for i in range(0,fluid_properties.Nc):
                    l_plus = np.copy(l[:,ph]); l_minus = np.copy(l[:,ph])
                    l_plus[i] = (nkphase[i,ph] + delta / 2) / Nphase[ph]
                    l_minus[i] = (nkphase[i,ph] - delta / 2) / Nphase[ph]
                    dlnf_dn[:,i,ph] = (PartialDerivatives.lnf(fluid_properties, l_plus, ph)
                     - PartialDerivatives.lnf(fluid_properties, l_minus, ph)) / delta
                    dZ_dn[i,ph] = (self.Z(fluid_properties, l_plus, ph)
                     - self.Z(fluid_properties, l_minus, ph)) / delta
        return dlnf_dn, dZ_dn

    def d_dP_all_derivatives(self, fluid_properties, Nphase, l, b):
        delta = 0.0001
        dZ_dP = np.zeros(self.n_phases)
        dlnf_dP = np.zeros([fluid_properties.Nc, self.n_phases])
        for ph in range(0, self.n_phases):
            if Nphase[ph] != 0:
                fluid_properties.P = self.Pvolume + delta/2
                Z_plus = self.Z(fluid_properties, l[:,ph], ph)
                lnf_plus = PartialDerivatives.lnf(fluid_properties, l[:,ph], ph)
                fluid_properties.P = self.Pvolume - delta / 2
                dZ_dP[ph] = (Z_plus - self.Z(fluid_properties, l[:,ph], ph)) / delta
                dlnf_dP[:,ph] = (lnf_plus - PartialDerivatives.lnf(fluid_properties, l[:,ph], ph)) / delta
        return dlnf_dP, dZ_dP

    def lnf(fluid_properties, l, ph):
        lnf = fluid_properties.lnphi(l, ph) + np.log(fluid_properties.P * l)
        return lnf

    def Z(self, fluid_properties, l, ph):
        A, B = fluid_properties.coefficientsPR(l)
        return StabilityCheck.Z_PR(B, A, ph)

    def dVt_derivatives(self, fluid_properties, Nphase_allvolumes, l_allvolumes, eta_j):
        n_vols = len(Nphase_allvolumes[0,0,:])

        """ Initializing dN derivatives """
        dnij_dNk = np.zeros([fluid_properties.Nc, fluid_properties.Nc, self.n_phases, n_vols])
        dZj_dnij = np.zeros([fluid_properties.Nc, self.n_phases, n_vols])

        """ Initializing dP derivatives """
        dnij_dP = np.zeros([fluid_properties.Nc, self.n_phases, n_vols])
        dZj_dP = np.zeros([1, self.n_phases, n_vols])
        P = np.copy(fluid_properties.P)
        Zj = np.zeros([1,self.n_phases, n_vols])
        for b in range(n_vols):
            self.Pvolume = P[b]
            fluid_properties.P = self.Pvolume 
            self.y = l_allvolumes[:,0,b]
            self.x = l_allvolumes[:,1,b]
            l = np.zeros([fluid_properties.Nc, self.n_phases])
            l[:,0] = self.y[0:fluid_properties.Nc]; l[:,1] = self.x[0:fluid_properties.Nc]
            Nphase = Nphase_allvolumes[0,:,b]
            dlnfij_dnij, dZj_dnij[:,:,b] = self.d_dn_all_derivatives(fluid_properties, Nphase, l)
            dlnfj_dP, dZj_dP[0,:,b] = self.d_dP_all_derivatives(fluid_properties, Nphase, l, b)
            Zj[0,:,b] = np.array([self.Z(fluid_properties,self.y,0), self.Z(fluid_properties, self.x,1)])
            matrix = np.sum(dlnfij_dnij, axis = 2)

            dlnf_dP_vector = dlnfj_dP[:,0] - dlnfj_dP[:,1]
            dnij_dP[:,1,b] =  np.linalg.inv(matrix)@dlnf_dP_vector
            dnij_dP[:,0,b] = - dnij_dP[:,1,b]

            for k in range(0, fluid_properties.Nc):
                dlnfl_dnk = dlnfij_dnij[:,k,1]
                dnij_dNk[:,k,1,b] = np.linalg.inv(matrix)@dlnfl_dnk
                dnij_dNk[:,k,0,b] = - dnij_dNk[:,k,1,b]
                dnij_dNk[k,k,0,b] = 1 - dnij_dNk[k,k,1,b]

        dvj_dnij = fluid_properties.R * fluid_properties.T / fluid_properties.P * dZj_dnij
        dVj_dNk = np.sum(dnij_dNk * (1 / eta_j + Nphase[:,np.newaxis] * dvj_dnij), axis = 0)
        dVt_dNk = np.sum(dVj_dNk, axis = 1)

        dvj_dP = fluid_properties.R * fluid_properties.T / fluid_properties.P * \
                (dZj_dP - Zj / fluid_properties.P)
        dVj_dP = np.sum(Nphase_allvolumes * dvj_dP + np.sum(dnij_dP * (1 / eta_j +
                Nphase[:,np.newaxis] * dvj_dnij), axis = 0),axis = 0)
        dVt_dP = np.sum(dVj_dP,axis = 0)

        fluid_properties.P = P #comming back
        return dVt_dNk, dVt_dP
