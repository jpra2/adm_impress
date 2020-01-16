from .. import directories as direc
import numpy as np

class BrooksAndCorey:

    def __init__(self):
        self.Sorw = float(direc.data_loaded['compositional_data']['Sorw'])
        self.Sorg = float(direc.data_loaded['compositional_data']['Sorg'])
        self.Sgr = float(direc.data_loaded['compositional_data']['Sgr'])
        self.Swr = float(direc.data_loaded['compositional_data']['Swr'])

        self.e_w = float(direc.data_loaded['compositional_data']['e_w'])
        self.e_o = float(direc.data_loaded['compositional_data']['e_o'])
        self.e_g = float(direc.data_loaded['compositional_data']['e_g'])

        self.krw0 = float(direc.data_loaded['compositional_data']['krw0'])
        self.kro0 = float(direc.data_loaded['compositional_data']['kro0'])
        self.krg0 = float(direc.data_loaded['compositional_data']['krg0'])

    def relative_permeabilities(self, Sw, So, Sg):
        self.Sor = self.Sorw * (1 - Sg / (1 - self.Swr - self.Sorg)) + \
                    self.Sorg * Sg / (1 - self.Swr - self.Sorg)
        kr = np.zeros(saturations.shape)
        kr[] = self.krw0 * ((Sw - self.Swr) / (1 - self.Swr - self.Sorw - self.Sgr)) ** self.e_w
        kro = self.kro0 * ((saturations[0,:] - self.Sor) / (1 - self.Swr - self.Sorw - self.Sgr)) ** self.e_o
        krg = self.krg0 * ((Sg - self.Sgr) / (1 - self.Swr - self.Sorw - self.Sgr)) ** self.e_g
        return krw, kro, krg

    def __call__(self, saturations):
        return self.relative_permeabilities(Sw, So, Sg)
    # def _stemp(self, S):
    #     return (S - self.Swc) / (1 - self.Swc - self.Sor)
    #
    # def _krw(self, S_temp):
    #     return self.krw0*(np.power(S_temp, self.n_w))
    #
    # def _kro(self, S_temp):
    #     return self.kro0*(np.power(1 - S_temp, self.n_o))
    #
    # def calculate(self, saturations):
    #
    #     n = len(saturations)
    #     ids = np.arange(n)
    #     ids_fora = ids[(saturations < 0) | (saturations > 1)]
    #     if len(ids_fora) > 0:
    #         raise ValueError('saturacao errada')
    #
    #     krw = np.zeros(n)
    #     kro = krw.copy()
    #     ids1 = ids[saturations > (1 - self.Sor)]
    #     ids2 = ids[saturations < self.Swc]
    #     ids_r = np.setdiff1d(ids, np.union1d(ids1, ids2))
    #
    #     krw[ids1] = np.repeat(self.krw0, len(ids1))
    #     kro[ids2] = np.repeat(self.kro0, len(ids2))
    #     Stemp = self._stemp(saturations[ids_r])
    #     krw[ids_r] = self._krw(Stemp)
    #     kro[ids_r] = self._kro(Stemp)
    #
    #     return krw, kro
