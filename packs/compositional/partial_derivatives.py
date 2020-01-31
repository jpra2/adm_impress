import numpy as np
from .stability_check import StabilityCheck

class FugacityParcialDerivative:

    def get_dlnphi_dnk (self, fluid_properties, l, Nphase, ph):
        """
        Numerical way to get lnphi derivative with respect to n component mols
        here we use a centered approach to the numerical derivative
        """
        nkphase = l * Nphase
        delta = 0.0001
        dlnphi_dnk = np.zeros([fluid_properties.Nc,fluid_properties.Nc])
        if Nphase != 0:
            for i in range(0,fluid_properties.Nc):
                nkphase_plus = np.copy(nkphase)
                nkphase_minus = np.copy(nkphase)
                nkphase_plus[i] = nkphase[i] + delta/2
                nkphase_minus[i] = nkphase[i] - delta/2
                dlnphi_dnk[:,i] = (self.lnphi(nkphase_plus, self.Nphase[ph], ph)
                 - self.lnphi(nkphase_minus, self.Nphase[ph], ph))/delta
        return dlnphi_dnk

    def lnphi(self, nkphase, Nphase, ph):
        l = nkphase / Nphase
        return StabilityCheck.lnphi(l,ph)

    def dlnphi_dn_matrix(self, fluid_properties, Nphase):
        dlnphil_dn_all = self.get_dlnphi_dnk(fluid_properties, fluid_properties.x, Nphase[1], 1)
        dlnphiv_dn_all = self.get_dlnphi_dnk(fluid_properties, fluid_properties.y, Nphase[0], 0)
        matrix = dlnphil_dn_all + dlnphiv_dn_all
        return matrix
