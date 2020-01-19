import numpy as np

## Attention to units: Temperature[K] and Pressure[atm]
# P  = 9.87E-6*P #[Pa] to [atm]
# Still going to add StoneII method so the user can choose
class LorenzBrayClark:
    def __init__(self,n_blocks,fluid_properties):
        self.n_blocks = n_blocks
        self.P = 9.87E-6*fluid_properties.P
        self.Vc = 1e3*fluid_properties.Vc

    def component_viscosity(self,T,fluid_properties):
        mi_components = np.zeros(fluid_properties.Nc)
        Trs = T/Tc
        ind_Tr_lower = np.argwhere(Trs <= 1.5)
        ind_Tr_higher = np.argwhere(Trs > 1.5)
        self.pure_components_viscosity_parameters = fluid_properties.Tc ** (1/6) \
                 / (fluid_properties.Mw ** (1/2) * fluid_properties.Pc ** (2/3))

        mi_components[ind_Tr_lower] = 3.4e-4 * Trs[ind_Tr_lower] ** 0.94 / \
                        self.pure_components_viscosity_parameters[ind_Tr_lower]
        mi_components[ind_Tr_higher] = (1.778e-4 * (4.58 * Trs[ind_Tr_higher] -
        1.67) ** (5/8)) / self.pure_components_viscosity_parameters[ind_Tr_higher]

        mi_components = mi_components[:, np.newaxis, np.newaxis]
        self.mi_components = np.ones([fluid_properties.Nc, 1, self.n_blocks]) * mi_components

    def mixture_viscosity(self,fluid_properties):
        self.component_molar_fractions = np.zeros(fluid_properties.Nc,2,self.n_blocks)
        self.component_molar_fractions[:,0,:] = fluid_properties.x
        self.component_molar_fractions[:,1,:] = fluid_properties.y

        mi_mix = np.sum(component_molar_fractions * self.mi_components *
                fluid_properties.Mw[:, np.newaxis, np.newaxis] ** (1/2), axis = 0) \
                /np.sum(component_molar_fractions *
                fluid_properties.Mw[:, np.newaxis, np.newaxis] , axis = 0)

        self.mi_mix = mi_mix[np.newaxis,:,:]

    def phase_viscosity(self,fluid_properties):
        #include in the entry parameters the vc: component critical molar volume
        # mi_phase = np.zeros([1,2,self.n_blocks])
        a = np.array([0.1023,0.023354,0.058533, -0.040758, 0.0093324])

        phase_reduced_rho = fluid_properties.phase_rho * \
            np.sum(self.component_molar_fractions * fluid_properties.Vc, axis = 0)

        # ind_lower = np.argwhere(phase_reduced_molar_density <= 0.18)
        # ind_higher = np.argwhere(phase_reduced_molar_density > 0.18)

        Xs = a[0]+ a[1] * phase_reduced_rho + a[2] * \
            phase_reduced_rho ** 2 - a[3] * \
            phase_reduced_rho ** 3 + a[4] * \
            phase_reduced_rho ** 4

        neta = (np.sum(self.component_molar_fractions *
                fluid_properties.Tc[:,np.newaxis,np.newaxis] , axis = 0)) \
                / (np.sum(self.component_molar_fractions *
                fluid_properties.Mw[:,np.newaxis,np.newaxis], axis = 0) ** (1/2)
                * np.sum(self.component_molar_fractions *
                fluid_properties.Pc[:,np.newaxis,np.newaxis], axis = 0) ** (2/3))

        # mi_phase[ind_lower] = self.mi_mix[ind_lower] + 2.05e-4 * \
        #                             self.zetas_r[ind_lower] / neta[ind_lower]
        # mi_phase[ind_higher] = self.mi_mix[ind_higher] + (Xs[ind_higher] ** 4 - 1) \
        #                         /(1e4*neta[ind_higher])
        mi_phase = (self.mi_mix + (Xs ** 4 - 1e-4) / neta)*1e-3 #return in Pa.s

        return mi_phase

    def __call__(self,fluid_properties):
        T = fluid_properties.T
        self.component_viscosity(T,fluid_properties)
        self.mixture_viscosity(self,fluid_properties)
        mi_phase = self.phase_viscosity(self,fluid_properties)
        return mi_phase
