import numpy as np
from . import constants as ctes

## Attention to units: Temperature[K] and Pressure[atm]
# P  = 9.87E-6*P #[Pa] to [atm]
# Still going to add StoneII method so the user can choose
class LorenzBrayClark:

    def __init__(self, fprop, phase_molar_densities):
        self.n_volumes = ctes.n_volumes
        self.phase_molar_densities = phase_molar_densities * 10**(-6) # mole/m³ to mole/cm³
        self.Mw = ctes.Mw * 10**(3)  # Kg/mole to grams/mole
        self.Pc = ctes.Pc / 101325 # Pa to atm
        self.vc = ctes.vc * 10 ** 6 # m³/mole to cm³/mole
        #self.vc[self.w>=0.489] = 21.573 + 0.015122 * self.Mw[self.w>=0.489] - \
        #                        27.656 * SG + 0.070615 * self.Mw[self.w>=0.489] * SG

    def component_viscosity(self, fprop):
        mi_components = np.zeros(ctes.Nc)
        Trs = fprop.T / ctes.Tc
        ind_Tr_lower = np.argwhere(Trs <= 1.5)
        ind_Tr_higher = np.argwhere(Trs > 1.5)

        # in this formula, the component molecular weight is in grams/mole. and
        # the pressure is in atm.
        # ctes.Mw and ctes.Pc are in kg/mole and Pa, respectively.
        pure_components_viscosity_parameters = ctes.Tc ** (1/6) \
                 / ((self.Mw) ** (1/2) * (self.Pc) ** (2/3))

        mi_components[ind_Tr_lower] = 3.4e-4 * Trs[ind_Tr_lower] ** 0.94 / \
                        pure_components_viscosity_parameters[ind_Tr_lower]
        mi_components[ind_Tr_higher] = (1.778e-4 * (4.58 * Trs[ind_Tr_higher] -
        1.67) ** (5/8)) / pure_components_viscosity_parameters[ind_Tr_higher]

        self.mi_components = mi_components[:, np.newaxis, np.newaxis]

    def phase_viscosity_atm(self, component_molar_fractions):
        '''Hernig and Zipperer equation'''
        self.component_molar_fractions = component_molar_fractions[0:ctes.Nc,0:2,:]

        self.mi_atm = np.sum(self.component_molar_fractions * self.mi_components *
                                self.Mw[:, np.newaxis, np.newaxis] ** (1/2), axis = 0) \
                        /np.sum(self.component_molar_fractions *
                                self.Mw[:, np.newaxis, np.newaxis] ** (1/2) , axis = 0)

        self.mi_atm = self.mi_atm[np.newaxis,:,:]

    def phase_viscosity(self):
        #include in the entry parameters the vc: component critical molar volume
        # mi_phase = np.zeros([1,2,self.n_volumes])
        a = np.array([0.1023, 0.023364, 0.058533, -0.040758, 0.0093324])

        phase_reduced_molar_density = (self.phase_molar_densities[:,0:2,:]) * \
            np.sum(self.component_molar_fractions * self.vc[:,np.newaxis,np.newaxis], axis = 0)

        # ind_lower = np.argwhere(phase_reduced_molar_density <= 0.18)
        # ind_higher = np.argwhere(phase_reduced_molar_density > 0.18)

        Xs = a[0]+ a[1] * phase_reduced_molar_density + a[2] * \
            phase_reduced_molar_density ** 2 + a[3] * \
            phase_reduced_molar_density ** 3 + a[4] * \
            phase_reduced_molar_density ** 4

        # parametro de viscosidade para mistura
        neta = (np.sum(self.component_molar_fractions *
                ctes.Tc[:,np.newaxis,np.newaxis] , axis = 0)) ** (1/6) \
                / (np.sum(self.component_molar_fractions *
                self.Mw[:,np.newaxis,np.newaxis], axis = 0) ** (1/2)
                * np.sum(self.component_molar_fractions *
                self.Pc[:,np.newaxis,np.newaxis], axis = 0) ** (2/3))

        # mi_phase[ind_lower] = self.mi_mix[ind_lower] + 2.05e-4 * \
        #                             self.zetas_r[ind_lower] / neta[ind_lower]
        # mi_phase[ind_higher] = self.mi_mix[ind_higher] + (Xs[ind_higher] ** 4 - 1) \
        #                         /(1e4*neta[ind_higher])

        mi_phase = (self.mi_atm + (Xs ** 4 - 1e-4) / neta) * 1e-3 #return in Pa.s
        return mi_phase

    def __call__(self, fprop, component_molar_fractions):
        self.component_viscosity(fprop)
        self.phase_viscosity_atm(component_molar_fractions)
        mi_phase = self.phase_viscosity()
        return mi_phase
