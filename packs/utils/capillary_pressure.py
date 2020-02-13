import numpy as np

class Chang:

    def __init__(self, data_loaded, phase_molar_densities, component_molar_fractions, pretransmissibility):
        self.Swr = data_loaded['compositional_data']['residual_saturations']['Swr']
        self.Sorg = data_loaded['compositional_data']['residual_saturations']['Sorg']
        self.Sorw = data_loaded['compositional_data']['residual_saturations']['Sorw']
        self.Sgr = data_loaded['compositional_data']['residual_saturations']['Sgr']
        self.Sorw = np.array(self.Sorw).astype(float)
        self.Sorg = np.array(self.Sorg).astype(float)
        self.phase_molar_densities = phase_molar_densities
        self.component_molar_fractions = component_molar_fractions
        self.porosity = data_impress['poro']
        self.permeabilities = pretransmissibility #?

    def capillary_pressure(self, sigma, S_term):
        Cpc = 16.125
        Epc = 1.1
        Pcap = - Cpc * sigma * (self.porosity / self.permeabilities) **(1/2) * (S_term) ** Epc
        return Pcap

    def sigma(self, i, j):
        parachor_number = 1
        sigma_ij = (0.001608 * np.sum(parachor_number * (self.phase_molar_densities[0,i,:] * self.component_molar_fractions[:,i,:]
                - self.phase_molar_densities[0,j,:] * self.component_molar_fractions[:,j,:]))) ** (1/0.25)
        return sigma_ij

    def normalized_saturation(self,S, Sr):
        S_norm = (S - Sr) / (1 - self.Swr - self.Sorw - self.Sgr)
        return S_norm

    def saturation_term(self, Sw, Sg, So):
        Sor = self.Sorw * (1 - Sg / (1 - self.Swr - self.Sorg)) + self.Sorg * \
            (Sg / (1 - self.Swr - self.Sorg))
        Sw_norm = self.normalized_saturation(Sw, self.Swr)
        Sg_norm = self.normalized_saturation(Sg, self.Sgr)
        So_norm = self.normalized_saturation(So, Sor)
        S_term_cow = 1 - Sw_norm
        S_term_cog = Sw_norm / (So_norm + Sg_norm)
        return S_term_cog, S_term_cow

    def __call__(self, Sw, Sg, So):
        S_term_cog, S_term_cow = self.saturation_term(Sw, Sg, So)
        oleo = 0; gas = 1; agua = 2
        sigma_wo = self.sigma(agua, oleo)
        sigma_og = self.sigma(oleo, gas)
        Pcow =  self.capillary_pressure(sigma_wo, S_term_cow)
        Pcog = self.capillary_pressure(sigma_og, S_term_cog)
        return Pcow, Pcog
