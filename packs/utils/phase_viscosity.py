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

    ''' Calculo analitico das derivadas da viscosidade em função de P e Nk '''
    def derivative_phase_viscosity_dP_dNk(self, fprop, dx_dP, dy_dP, dCsi_j_dP, dx_dnij, dy_dnij, dCsi_j_dnij, dnildP, dnivdP, dnildNk, dnivdNk):
        self.component_viscosity(fprop)
        self.phase_viscosity_atm(fprop.xkj)
        dxij_dP = np.zeros_like(self.component_molar_fractions)
        dxij_dP[:,0,:] = dx_dP[0:ctes.Nc]
        dxij_dP[:,1,:] = dy_dP[0:ctes.Nc]
        dxij_dnij = np.zeros([ctes.Nc, ctes.Nc, 2, ctes.n_volumes])
        dxij_dnij[:,:,0] = dx_dnij
        dxij_dnij[:,:,1] = dy_dnij
        dnij_dP = np.zeros([ctes.Nc, 2, ctes.n_volumes])
        dnij_dP[:,0] = dnildP
        dnij_dP[:,1] = dnivdP
        dnij_dNk = np.zeros([ctes.Nc, ctes.Nc, 2, ctes.n_volumes])
        dnij_dNk[:,:,0] = dnildNk
        dnij_dNk[:,:,1] = dnivdNk

        d_mi_atm_dP, d_mi_atm_dnij = self.derivative_mi_atm_dP_dnij(dxij_dP, dxij_dnij)
        a = np.array([0.1023, 0.023364, 0.058533, -0.040758, 0.0093324])
        phase_reduced_molar_density = (self.phase_molar_densities[:,0:2,:]) * \
            np.sum(self.component_molar_fractions * self.vc[:,np.newaxis,np.newaxis], axis = 0)
        Xs = a[0]+ a[1] * phase_reduced_molar_density + a[2] * \
            phase_reduced_molar_density ** 2 + a[3] * \
            phase_reduced_molar_density ** 3 + a[4] * \
            phase_reduced_molar_density ** 4
        dXs4_dP, dXs4_dnij = self.dXs4_dP_dnij(dCsi_j_dP, dCsi_j_dnij, dxij_dP, dxij_dnij, a, phase_reduced_molar_density, Xs)
        dneta_dP, dneta_dnij = self.dneta_dP_dnij(dxij_dP, dxij_dnij)
        neta = (np.sum(self.component_molar_fractions *
                ctes.Tc[:,np.newaxis,np.newaxis] , axis = 0)) ** (1/6) \
                / (np.sum(self.component_molar_fractions *
                self.Mw[:,np.newaxis,np.newaxis], axis = 0) ** (1/2)
                * np.sum(self.component_molar_fractions *
                self.Pc[:,np.newaxis,np.newaxis], axis = 0) ** (2/3))

        dmi_phase_dP = (d_mi_atm_dP + ((dXs4_dP*neta - (Xs**4 - 1e-4)*dneta_dP) / \
                         (neta**2))) * 1e-3

        dmi_phase_dnij = (d_mi_atm_dnij + ((dXs4_dnij*neta[np.newaxis] - \
            (Xs**4 - 1e-4)*dneta_dnij) / (neta[np.newaxis]**2))) * 1e-3

        dmi_dP = np.sum(dnij_dP * dmi_phase_dnij, axis = 0) + dmi_phase_dP
        if ctes.load_w: dmi_dP = np.append(dmi_dP, np.zeros([1,1,ctes.n_volumes]), axis=1)

        dmi_dNk_aux = np.sum(dnij_dNk * dmi_phase_dnij[:,np.newaxis], axis = 0)
        dmi_dNk = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        dmi_dNk[0:ctes.Nc, 0:2] = dmi_dNk_aux
        return dmi_dP, dmi_dNk

    def derivative_mi_atm_dP_dnij(self, dxij_dP, dxij_dnij):
        d_mi_atm_dP = (np.sum(dxij_dP * self.mi_components * \
            self.Mw[:, np.newaxis, np.newaxis] ** (1/2), axis = 0)*\
            np.sum(self.component_molar_fractions*self.Mw[:, np.newaxis, np.newaxis] ** (1/2) , axis = 0)\
            - np.sum(self.component_molar_fractions * self.mi_components *\
            self.Mw[:, np.newaxis, np.newaxis] ** (1/2), axis = 0)*\
            np.sum(dxij_dP * self.Mw[:, np.newaxis, np.newaxis] ** (1/2), axis = 0)) / \
            ((np.sum(self.component_molar_fractions * \
            self.Mw[:, np.newaxis, np.newaxis] ** (1/2) , axis = 0))**2)

        d_mi_atm_dnij = (np.sum(dxij_dnij * self.mi_components[:,np.newaxis] * \
               self.Mw[:, np.newaxis, np.newaxis, np.newaxis] ** (1/2), axis = 0)*\
               np.sum(self.component_molar_fractions*self.Mw[:, np.newaxis, np.newaxis] ** \
                      (1/2) , axis = 0)[np.newaxis]\
               - np.sum(self.component_molar_fractions * self.mi_components *\
               self.Mw[:, np.newaxis, np.newaxis] ** (1/2), axis = 0)[np.newaxis]*\
               np.sum(dxij_dnij * self.Mw[:, np.newaxis, np.newaxis, np.newaxis] ** (1/2), axis = 0)) / \
               ((np.sum(self.component_molar_fractions * \
               self.Mw[:, np.newaxis, np.newaxis] ** (1/2) , axis = 0))**2)[np.newaxis]

        return d_mi_atm_dP, d_mi_atm_dnij

    def dXs4_dP_dnij(self, dCsi_j_dP, dCsi_j_dnij, dxij_dP, dxij_dnij, a, phase_reduced_molar_density, Xs):
        dCsi_j_dP = dCsi_j_dP * (10**-6) # mole/m³.Pa to mole/cm³.Pa, fica na mesma unidade do resto do modelo
        dCsi_j_dnij = dCsi_j_dnij * (10**-6) # mole/m³.mol to mole/cm³.mol, fica na mesma unidade do resto do modelo

        dphase_reduced_molar_density_dP = dCsi_j_dP[:,0:2,:]*\
            np.sum(self.component_molar_fractions * self.vc[:,np.newaxis,np.newaxis], axis = 0)\
            + self.phase_molar_densities[:,0:2,:] * np.sum(dxij_dP * \
            self.vc[:,np.newaxis,np.newaxis], axis = 0)

        dphase_reduced_molar_density_dnij = dCsi_j_dnij[:,0:2] * \
            np.sum(self.component_molar_fractions * self.vc[:,np.newaxis,np.newaxis], axis = 0)[np.newaxis]\
            + self.phase_molar_densities[:,0:2,:] * np.sum(dxij_dnij * \
            self.vc[:,np.newaxis,np.newaxis,np.newaxis], axis = 0)

        dXs_dP = a[1]*dphase_reduced_molar_density_dP + \
            a[2]*2*phase_reduced_molar_density*dphase_reduced_molar_density_dP + \
            a[3]*3*(phase_reduced_molar_density**2)*dphase_reduced_molar_density_dP + \
            a[4]*4*(phase_reduced_molar_density**3)*dphase_reduced_molar_density_dP

        dXs_dnij = a[1]*dphase_reduced_molar_density_dnij + \
            a[2]*2*phase_reduced_molar_density*dphase_reduced_molar_density_dnij + \
            a[3]*3*(phase_reduced_molar_density**2)*dphase_reduced_molar_density_dnij + \
            a[4]*4*(phase_reduced_molar_density**3)*dphase_reduced_molar_density_dnij

        dXs4_dP = 4 * (Xs**3) * dXs_dP
        dXs4_dnij = 4 * (Xs**3) * dXs_dnij
        return dXs4_dP, dXs4_dnij

    def dneta_dP_dnij(self, dxij_dP, dxij_dnij):
        aux = (np.sum(self.component_molar_fractions *
            self.Mw[:,np.newaxis,np.newaxis], axis = 0) ** (1/2)
            * np.sum(self.component_molar_fractions *
            self.Pc[:,np.newaxis,np.newaxis], axis = 0) ** (2/3))

        somat_xkj_Tc = np.sum(self.component_molar_fractions *
                ctes.Tc[:,np.newaxis,np.newaxis] , axis = 0)
        somat_xkj_Mw = np.sum(self.component_molar_fractions *
            self.Mw[:,np.newaxis,np.newaxis], axis = 0)
        somat_xkj_Pc = np.sum(self.component_molar_fractions *
            self.Pc[:,np.newaxis,np.newaxis], axis = 0)

        dneta_dP = ((1/6)*((somat_xkj_Tc)**(-5/6))*np.sum(dxij_dP* \
                ctes.Tc[:,np.newaxis,np.newaxis], axis = 0)*aux - (somat_xkj_Tc**(1/6))* \
                ((1/2)*(somat_xkj_Mw**(-1/2))* np.sum(dxij_dP* \
                self.Mw[:,np.newaxis,np.newaxis], axis = 0)*(somat_xkj_Pc**(2/3)) + \
                (somat_xkj_Mw**0.5)*(2/3)*(somat_xkj_Pc**(-1/3))*np.sum(dxij_dP* \
                self.Pc[:,np.newaxis,np.newaxis], axis = 0))) / (aux**2)

        dneta_dnij = ((1/6)*((somat_xkj_Tc)**(-5/6))*np.sum(dxij_dnij* \
                ctes.Tc[:,np.newaxis,np.newaxis,np.newaxis], axis = 0)*aux - (somat_xkj_Tc**(1/6))* \
                ((1/2)*(somat_xkj_Mw**(-1/2))* np.sum(dxij_dnij* \
                self.Mw[:,np.newaxis,np.newaxis,np.newaxis], axis = 0)*(somat_xkj_Pc**(2/3)) + \
                (somat_xkj_Mw**0.5)*(2/3)*(somat_xkj_Pc**(-1/3))*np.sum(dxij_dnij* \
                self.Pc[:,np.newaxis,np.newaxis,np.newaxis], axis = 0))) / (aux**2)

        return dneta_dP, dneta_dnij
