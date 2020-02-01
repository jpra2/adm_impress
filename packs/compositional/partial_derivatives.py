import numpy as np
from .stability_check import StabilityCheck

class ParcialDerivatives:

    def __init__(self):
        self.n_phases = 2

    def d_dn_all_derivatives (self, fluid_properties, Nphase):
        """
        Numerical way to get lnphi and Z derivative with respect to n component
        mols.
        Here we use a centered approach to the numerical derivative
        """
        l = np.zeros([fluid_properties.Nc,2])
        l[:,0] = self.y; l[:,1] = self.x
        nkphase = l * Nphase
        delta = 0.0001
        dlnphi_dn = np.zeros([fluid_properties.Nc, fluid_properties.Nc,2])
        dZ_dn = np.zeros([fluid_properties.Nc, 2])
        for ph in range(2):
            if Nphase[ph] != 0:
                for i in range(0,fluid_properties.Nc):
                    l_plus = np.copy(l[:,ph]); l_minus = np.copy(l[:,ph])
                    l_plus[i] = (nkphase[i,ph] + delta/2) / Nphase[ph]
                    l_minus[i] = (nkphase[i,ph] - delta/2) / Nphase[ph]
                    dlnphi_dn[:,i,ph] = (StabilityCheck.lnphi(l_plus, ph)
                     - StabilityCheck.lnphi(l_minus, ph)) / delta
                    dZ_dn[i,ph] = (self.Z(l_plus, ph)
                     - self.Z(l_minus, ph)) / delta
                # dlnphi_dn[:]
        return dlnphi_dn, dZ_dn

    def d_dP_all_derivatives(self, fluid_properties, Nphase):
        delta = 0.0001
        l = np.zeros([fluid_properties.Nc, self.n_phases])
        P = np.copy(fluid_properties.P); P_plus = np.copy(P); P_minus = np.copy(P)
        dZ_dP = np.zeros(self.n_phases)
        dlnphi_dP = np.zeros([fluid_properties.Nc, self.n_phases])
        l[:,0] = self.y[:,b]; l[:,1] = self.x[:,b]
        for ph in range(0, self.n_phases):
            if Nphase[ph] != 0:
                self.P = P + delta/2
                Z_plus = self.Z(l[:,ph], ph)
                lnphi_plus = StabilityCheck.lnphi(l[:,ph], ph)
                self.P = P - delta / 2
                dZ_dP[ph] = (Z_plus - self.Z(l[:,ph], ph)) / delta
                dlnphi_dP[:,ph] = (lnphi_plus - StabilityCheck.lnphi(l[:,ph], ph)) / delta
        self.P = P #voltando ao original
        return dlnphi_dP, dZ_dP

    def Z(self, nkphase, Nphase, ph):
        l = nkphase / Nphase
        A, B = StabilityCheck.coefficientsPR(l)
        return StabilityCheck.Z_PR(B, A, ph)

    def dVt_derivatives(self, fluid_properties, Nphase_all, l_allvolumes, eta_j):
        n_vols = len(Nphase_all[0,0,:])

        """ Initializing dN derivatives """
        dnij_dNk = np.zeros([fluid_properties.Nc, fluid_properties.Nc, self.n_phases, n_vols])
        dZj_dnij = np.zeros([fluid_properties.Nc, self.n_phases, n_vols])

        """ Initializing dP derivatives """
        dnij_dP = np.zeros([fluid_properties.Nc, self.n_phases, n_vols])
        dZj_dP = np.zeros([1, self.n_phases, n_vols])

        Zj = np.zeros([1,self.n_phases, n_vols])

        for b in range(n_vols):
            self.y = l_allvolumes[:,0,b]
            self.x = l_allvolumes[:,1,b]
            Nphase = Nphase_allvolumes[0,:,b]
            dlnphiij_dnij, dZj_dnij[:,:,b] = self.d_dn_all_derivatives(fluid_properties, Nphase)
            dlnphij_dP, dZj_dP[0,:,b] = self.d_dP_all_derivatives(fluid_properties, Nphase)
            Zj[0,:,b] = np.array([self.Z(self.y,0), self.Z(self.x,1)])
            matrix = np.sum(dlnphiij_dnij, axis = 2)

            dlnphi_dP_vector = dlnphij_dP[:,0] - dlnphij_dP[:,1]
            dnij_dP[:,1,b] =  np.linalg.inv(matrix)@dlnphi_dP_vector
            dnij_dP[:,0,b] = - dnij_dP[:,1,b]

            for k in range(0, self.Nc):
                dlnphil_dnk = dlnphi_dn[:,k,1]
                dnij_dNk[:,k,1,b] = np.linalg.inv(matrix)@dlnphil_dnk
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
        return dVt_dNk, dVt_dP
