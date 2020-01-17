

class LorenzBrayClark:
    def __init__(self,n_blocks):
        self.n_blocks = n_blocks

    def component_viscosity(self,T,fluid_properties):
        mi_components = np.zeros(fluid_properties.Nc)
        Trs = T/Tc
        ind_Tr_lower = np.argwhere(Trs <= 0.15)
        ind_Tr_higher = np.argwhere(Trs > 0,15)
        self.zetas = 5.44 * fluid_properties.Tc ** (1/6) / (fluid_properties.Mw ** (1/2)
                * fluid_properties.Pc ** (2/3))

        mi_components[ind_Tr_lower] = 3.4e-4 * Trs[ind_Tr_lower] ** 0.94 / \
                                        self.zetas[ind_Tr_lower]
        mi_components[ind_Tr_higher] = (1.778e-4 * (4.58 * Trs[ind_Tr_higher] -
                                        1.67) ** (5/8)) / self.zetas[ind_Tr_higher]

        mi_components = mi_components[:, np.newaxis, np.newaxis]
        self.mi_components = np.ones([fluid_properties.Nc, 1, self.n_blocks]) * mi_components

    def mixture_viscosity(self,fluid_properties):
        self.component_molar_fractions = np.zeros(fluid_properties.Nc,2,self.n_blocks)
        self.component_molar_fractions[:,0,:] = fluid_properties.component_molar_fractions[:,0,:]
        self.component_molar_fractions[:,1,:] = fluid_properties.component_molar_fractions[:,1,:]
        mi_mix = np.sum(component_molar_fractions*self.mi_components *
                fluid_properties.Mw[:,np.newaxis,np.newaxis] ** (1/2) , axis = 0) \
                /np.sum(component_molar_fractions *
                fluid_properties.Mw[:,np.newaxis,np.newaxis] , axis = 0)
        self.mi_mix = mi_mix[np.newaxis,:,:]

    def phase_viscosity(self,fluid_properties):
        #include in the entry parameters the vc: component critical molar volume??
        mi_phase = np.zeros([1,2,self.n_blocks])

        phase_reduced_molar_density = fluid_properties.phase_molar_densities \
            *np.sum(self.component_molar_fractions*fluid_properties.vc,axis = 0)

        ind_lower = np.argwhere(phase_reduced_molar_density <= 0.18)
        ind_higher = np.argwhere(phase_reduced_molar_density > 0.18)

        Xs = 1.023 + 0.23364 * phase_reduced_molar_density + 0.58533 * \
            phase_reduced_molar_density ** 2 - 0.40758 * \
            phase_reduced_molar_density ** 3 + 0.093324 * \
            phase_reduced_molar_density ** 4

        neta = 5.44 * (np.sum(self.component_molar_fractions *
                        fluid_properties.Tc[:,np.newaxis,np.newaxis] , axis = 0)) \
                        / (np.sum(self.component_molar_fractions *
                        fluid_properties.Mw[:,np.newaxis,np.newaxis], axis = 0) ** (1/2)
                        * np.sum(self.component_molar_fractions *
                        fluid_properties.Pc[:,np.newaxis,np.newaxis], axis = 0) ** (2/3))

        mi_phase[ind_lower] = self.mi_mix[ind_lower] + 2.05e-4 * \
                                    self.zetas_r[ind_lower] / neta[ind_lower]
        mi_phase[ind_higher] = self.mi_mix[ind_higher] + (Xs[ind_higher] ** 4 - 1) \
                                /(10e4*neta[ind_higher])
        return mi_phase

    def __call__(self,T,fluid_properties):
        self.component_viscosity(T,fluid_properties)
        self.mixture_viscosity(self,fluid_properties)
        mi_phase = self.phase_viscosity(self,fluid_properties)
        return mi_phase
