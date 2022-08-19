import numpy as np
from . import constants as ctes
from packs.directories import data_loaded

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
        self.component_molar_fractions = component_molar_fractions[0:ctes.Nc,0:ctes.n_phases,:]

        self.mi_atm = np.sum(self.component_molar_fractions * self.mi_components *
                                self.Mw[:, np.newaxis, np.newaxis] ** (1/2), axis = 0) \
                        /np.sum(self.component_molar_fractions *
                                self.Mw[:, np.newaxis, np.newaxis] ** (1/2) , axis = 0)

        self.mi_atm = self.mi_atm[np.newaxis,:,:]

    def phase_viscosity(self, fprop):
        #include in the entry parameters the vc: component critical molar volume
        # mi_phase = np.zeros([1,2,self.n_volumes])
        a = np.array([0.1023, 0.023364, 0.058533, -0.040758, 0.0093324])

        phase_reduced_molar_density = (self.phase_molar_densities[:,0:ctes.n_phases,:]) * \
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
        # Viscosidade da água
        Csw = data_loaded['compositional_data']['water_data']['Csw']
        a11 = 3.324e-2
        a12 = 3.624e-3
        a13 = - 1.879e-4
        a21 = - 3.96e-2
        a22 = 1.02e-2
        a23 = - 7.02e-4
        a31 = 1.2378
        a32 = - 1.303e-3
        a33 = 3.06e-6
        a34 = 2.55e-8
        mi_w20_0 = 1.002 # cp
        T_mi = (fprop.T - 273.15)  # K em C
        p_mi = fprop.P*1e-6 # Pa to MPa

        A1 = a11*Csw + a12*Csw**2 + a13*Csw**3
        A2 = a21*Csw + a22*Csw**2 + a23*Csw**3
        S_ai = a31*(20-T_mi)/(96+T_mi) + a32*(20-T_mi)**2/(96+T_mi) + a33*(20-T_mi)**3/(96+T_mi)
        mi_w0 = np.exp(S_ai + np.log(mi_w20_0))
        mi_ast_w = np.exp(A1 + A2*np.log(mi_w0/mi_w20_0) + np.log(mi_w0))

        A0 = 1e-3*(0.8+0.01*(T_mi-90)*np.exp(-0.25*Csw))
        miw = (1+A0*p_mi)*mi_ast_w*1e-3
        mi_phase[:,ctes.n_phases-1,:] = miw

        #import pdb; pdb.set_trace()
        return mi_phase

    def __call__(self, fprop, component_molar_fractions):
        self.component_viscosity(fprop)
        self.phase_viscosity_atm(component_molar_fractions)
        mi_phase = self.phase_viscosity(fprop)
        return mi_phase
