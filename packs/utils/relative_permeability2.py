from .. import directories as direc
import numpy as np

class BrooksAndCorey:

    def __init__(self):
        self.Sorw = float(direc.data_loaded['compositional_data']['Sorw'])
        self.Sorg = float(direc.data_loaded['compositional_data']['Sorg'])
        self.Sgr = float(direc.data_loaded['compositional_data']['Sgr'])
        self.Swr = float(direc.data_loaded['compositional_data']['Swr'])

        self.n_w = float(direc.data_loaded['compositional_data']['n_w'])
        self.n_o = float(direc.data_loaded['compositional_data']['n_o'])
        self.n_g = float(direc.data_loaded['compositional_data']['n_g'])

        self.krw0 = float(direc.data_loaded['compositional_data']['krw0'])
        self.kro0 = float(direc.data_loaded['compositional_data']['kro0'])
        self.krg0 = float(direc.data_loaded['compositional_data']['krg0'])

    def relative_permeabilities(self, Sw, So, Sg):
        self.Sor = self.Sorw * (1 - Sg / (1 - self.Swr - self.Sorg)) + \
                    self.Sorg * Sg / (1 - self.Swr - self.Sorg)
        kr = np.zeros(saturations.shape)
        krw = self.krw0 * ((Sw - self.Swr) / (1 - self.Swr - self.Sorw - self.Sgr)) ** self.n_w
        kro = self.kro0 * ((So - self.Sor) / (1 - self.Swr - self.Sorw - self.Sgr)) ** self.n_o
        krg = self.krg0 * ((Sg - self.Sgr) / (1 - self.Swr - self.Sorw - self.Sgr)) ** self.n_g
        return krw, kro, krg

    def __call__(self, So, Sg, Sw):
        return self.relative_permeabilities(Sw, So, Sg)
