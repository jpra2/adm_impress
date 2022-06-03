from .stability_check import StabilityCheck as SC
from ..utils import constants as ctes
import numpy as np

class StabilityCheck(SC):
    def run_init(self, P, z):
        self.L, self.V, self.x, self.y, ksi_L, ksi_V, rho_L, rho_V = self.run_constant_K(z)
        return self.L, self.V, self.x, self.y, ksi_L, ksi_V, rho_L, rho_V

    def run(self, P, z, wells):
        self.L, self.V, self.x, self.y, ksi_L, ksi_V, rho_L, rho_V = self.run_constant_K(z)
        return self.L, self.V, self.x, self.y, ksi_L, ksi_V, rho_L, rho_V

    def run_constant_K(self, z):
        self.z = np.copy(z)
        self.K = self.K * np.ones_like(self.z[0,:])
        Lmax = np.max(self.K, axis = 0)/(np.max(self.K, axis = 0) - 1)
        Lmin = np.min(self.K, axis = 0)/(np.min(self.K, axis = 0) - 1)

        Vmax = 1. - Lmin
        Vmin = 1. - Lmax

        self.z[self.z==0] = 1e-300

        #Vmin = ((K1-KNc)*z[self.K==K1]-(1-KNc))/((1-KNc)*(K1-1))
        #proposed by Li et al for Whitson method
        Vmin[Vmin>Vmax] = -Vmax[Vmin>Vmax]
        self.V = (Vmin + Vmax) * 0.5
        ponteiro = np.ones_like(self.V,dtype=bool)

        self.V, self.x, self.y = self.solve_objective_function_Whitson_for_V(np.copy(self.V), Vmax, \
            Vmin, np.copy(self.K), self.z)
        #self.Yinghui_method(ponteiro) #ajeitar!
        #import pdb; pdb.set_trace()
        self.z[z==0] = 0#z[z==0]

        self.x[:,((self.V)<=0) + ((self.V)>=1)] = self.z[:,((self.V)<=0) + ((self.V)>=1)]
        self.y[:,((self.V)<=0) + ((self.V)>=1)] = self.z[:,((self.V)<=0) + ((self.V)>=1)]
        self.V[self.V<0] = 0
        self.V[self.V>1] = 1
        self.L = 1 - self.V

        #ksi = np.array([37342.0019279, 37138.91334958, 13792.42036739,  5248.87665093, 3013.74120719])
        #ksi = np.array([37342.0019279, 37342.0019279, 37342.0019279, 37342.0019279 , 37342.0019279 ])
        ksi = 37342.0019279 * np.ones_like(self.z[:,0])

        #ksi = np.array([36948.96102266, 13766.00368349, 5544.36875657])
        ctes.Mw[:] = 10

        '''lnphil, Zl = self.lnphi_Z_based_on_deltaG(self.x, self.P, self.ph_L)
        lnphiv, Zv = self.lnphi_Z_based_on_deltaG(self.y, self.P, self.ph_V)
        ksi_L, rho_L = self.get_EOS_dependent_properties(self.T, self.x, self.P, Zl)
        ksi_V, rho_V = self.get_EOS_dependent_properties(self.T, self.y, self.P, Zv)
        '''
        ksi_L = np.sum(self.x / ksi[:,np.newaxis], axis=0)
        ksi_V = np.sum(self.y / ksi[:,np.newaxis], axis=0)
        Mw_L = np.sum(self.x * ctes.Mw[:,np.newaxis], axis = 0)
        Mw_V = np.sum(self.y * ctes.Mw[:,np.newaxis], axis = 0)
        rho_L = ksi_L * Mw_L
        rho_V = ksi_V * Mw_V
        return self.L, self.V, self.x, self.y, ksi_L, ksi_V, rho_L, rho_V
