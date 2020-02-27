import numpy as np

## Attention to units: Temperature[K] and Pressure[atm]
# P  = 9.87E-6*P #[Pa] to [atm]
# Still going to add StoneII method so the user can choose
class LorenzBrayClark:

    def __init__(self, n_volumes, fprop, kprop):
        self.n_volumes = n_volumes
        self.P = 9.87E-6*fprop.P
        self.Vc = 1e3*kprop.Vc

    def component_viscosity(self, fprop, kprop):
        mi_components = np.zeros(fprop.Nc)
        Trs = fprop.T / kprop.Tc
        ind_Tr_lower = np.argwhere(Trs <= 1.5)
        ind_Tr_higher = np.argwhere(Trs > 1.5)
        self.pure_components_viscosity_parameters = kprop.Tc ** (1/6) \
                 / (kprop.Mw ** (1/2) * kprop.Pc ** (2/3))

        mi_components[ind_Tr_lower] = 3.4e-4 * Trs[ind_Tr_lower] ** 0.94 / \
                        self.pure_components_viscosity_parameters[ind_Tr_lower]
        mi_components[ind_Tr_higher] = (1.778e-4 * (4.58 * Trs[ind_Tr_higher] -
        1.67) ** (5/8)) / self.pure_components_viscosity_parameters[ind_Tr_higher]

        mi_components = mi_components[:, np.newaxis, np.newaxis]
        self.mi_components = np.ones([fprop.Nc, 1, self.n_volumes]) * mi_components

    def mixture_viscosity(self, fprop, kprop):
        self.component_molar_fractions = np.zeros([fprop.Nc,2,self.n_volumes])
        self.component_molar_fractions[:,0,:] = fprop.x
        self.component_molar_fractions[:,1,:] = fprop.y

        mi_mix = np.sum(self.component_molar_fractions * self.mi_components *
                kprop.Mw[:, np.newaxis, np.newaxis] ** (1/2), axis = 0) \
                /np.sum(self.component_molar_fractions *
                kprop.Mw[:, np.newaxis, np.newaxis] , axis = 0)

        self.mi_mix = mi_mix[np.newaxis,:,:]

    def phase_viscosity(self, fprop, kprop):
        #include in the entry parameters the vc: component critical molar volume
        # mi_phase = np.zeros([1,2,self.n_volumes])
        a = np.array([0.1023, 0.023354, 0.058533, -0.040758, 0.0093324])
        self.phase_mass_densities = np.zeros([1, 2, self.n_volumes])
        self.phase_mass_densities[0,0,:] = fprop.rho_L
        self.phase_mass_densities[0,1,:] = fprop.rho_V

        phase_reduced_rho = self.phase_mass_densities * \
            np.sum(self.component_molar_fractions * kprop.Vc, axis = 0)

        # ind_lower = np.argwhere(phase_reduced_molar_density <= 0.18)
        # ind_higher = np.argwhere(phase_reduced_molar_density > 0.18)

        Xs = a[0]+ a[1] * phase_reduced_rho + a[2] * \
            phase_reduced_rho ** 2 - a[3] * \
            phase_reduced_rho ** 3 + a[4] * \
            phase_reduced_rho ** 4

        neta = (np.sum(self.component_molar_fractions *
                kprop.Tc[:,np.newaxis,np.newaxis] , axis = 0)) \
                / (np.sum(self.component_molar_fractions *
                kprop.Mw[:,np.newaxis,np.newaxis], axis = 0) ** (1/2)
                * np.sum(self.component_molar_fractions *
                kprop.Pc[:,np.newaxis,np.newaxis], axis = 0) ** (2/3))

        # mi_phase[ind_lower] = self.mi_mix[ind_lower] + 2.05e-4 * \
        #                             self.zetas_r[ind_lower] / neta[ind_lower]
        # mi_phase[ind_higher] = self.mi_mix[ind_higher] + (Xs[ind_higher] ** 4 - 1) \
        #                         /(1e4*neta[ind_higher])
        mi_phase = (self.mi_mix + (Xs ** 4 - 1e-4) / neta)*1e-3 #return in Pa.s

        return mi_phase

    def __call__(self, fprop, kprop):
        T = fprop.T
        self.component_viscosity(fprop, kprop)
        self.mixture_viscosity(fprop, kprop)
        mi_phase = self.phase_viscosity(fprop, kprop)
        return mi_phase
