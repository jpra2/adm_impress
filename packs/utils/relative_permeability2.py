from .. import directories as direc
import numpy as np

class BrooksAndCorey:

    def __init__(self):
        self.Sorw = float(direc.data_loaded['compositional_data']['residual_saturations']['Sorw'])
        self.Sorg = float(direc.data_loaded['compositional_data']['residual_saturations']['Sorg'])
        self.Sgr = float(direc.data_loaded['compositional_data']['residual_saturations']['Sgr'])
        self.Swr = float(direc.data_loaded['compositional_data']['residual_saturations']['Swr'])

        self.n_w = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_w'])
        self.n_o = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_o'])
        self.n_g = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_g'])

        # End-point relative permeability data
        self.krw0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krw0'])
        self.kro0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['kro0'])
        self.krg0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krg0'])

    def relative_permeabilities(self, saturations):
        #saturations = [So,Sg,Sw]
        Sor = self.Sorw * (1 - saturations[1] / (1 - self.Swr - self.Sorg)) + \
                    self.Sorg * saturations[1] / (1 - self.Swr - self.Sorg)

        krw = self.krw0 * ((saturations[2] - self.Swr) / (1 - self.Swr - Sor - self.Sgr)) ** self.n_w
        krw[saturations[2]<self.Swr] = 0
        kro = self.kro0 * ((saturations[0] - Sor) / (1 - self.Swr - Sor - self.Sgr)) ** self.n_o
        #kro[saturations[2]<self.Swr] = self.kro0 * ((saturations[0][saturations[2]<self.Swr] - self.Sor) / (1 - self.Sor - self.Sgr)) ** self.n_o
        kro[saturations[2]<self.Swr] = 1
        krg = self.krg0 * ((saturations[1] - self.Sgr) / (1 - self.Swr - Sor - self.Sgr)) ** self.n_g
        return kro, krg, krw, Sor

    def __call__(self, saturations):
        return self.relative_permeabilities(saturations)

class StoneII:

    def __init__(self):
        self.Sorw = float(direc.data_loaded['compositional_data']['residual_saturations']['Sorw'])
        self.Sorg = float(direc.data_loaded['compositional_data']['residual_saturations']['Sorg'])
        self.Sgr = float(direc.data_loaded['compositional_data']['residual_saturations']['Sgr'])
        self.Swr = float(direc.data_loaded['compositional_data']['residual_saturations']['Swr'])

        self.n_w = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_w'])
        self.n_ow = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_ow'])
        self.n_og = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_og'])
        self.n_g = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_g'])

        # End-point relative permeability data
        self.krw0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krw0'])
        self.krow0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krow0'])
        self.krog0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krog0'])
        self.krg0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krg0'])

    def relative_permeabilities(self, saturations):
        #saturations = [So,Sg,Sw]
        Sor = self.Sorw * (1 - saturations[1] / (1 - self.Swr - self.Sorg)) + \
                    self.Sorg * saturations[1] / (1 - self.Swr - self.Sorg)

        krw = self.krw0 * ((saturations[2] - self.Swr) / (1 - self.Swr - self.Sorw)) ** self.n_w
        krw[saturations[2]<self.Swr] = 0
        krg = self.krg0 * ((saturations[1] - self.Sgr) / (1 - self.Swr - self.Sorg - self.Sgr)) ** self.n_g

        krow = self.krow0 * ((1 - saturations[2] - self.Sorw) / (1 - self.Swr - self.Sorw)) ** self.n_ow
        krog = self.krog0 * ((1 - saturations[1] - self.Sorg - self.Swr) / (1 - self.Swr - self.Sgr - self.Sorg))

        kro = self.krow0 *((krow/self.krow0 + krw) * (krog/self.krow0 + krg) - (krw + krg))
        #kro[saturations[2]<self.Swr] = self.kro0 * ((saturations[0][saturations[2]<self.Swr] - self.Sor) / (1 - self.Sor - self.Sgr)) ** self.n_o
        kro[saturations[2]<self.Swr] = 1

        return kro, krg, krw, Sor

    def __call__(self, saturations):
        return self.relative_permeabilities(saturations)
