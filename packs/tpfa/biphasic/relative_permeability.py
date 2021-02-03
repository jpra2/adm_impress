import numpy as np
import pdb

class BrooksAndCorey:

    def __init__(self, data_loaded):
        self.Sor = float(data_loaded['Sor'])
        self.Swc = float(data_loaded['Swc'])
        self.n_w = float(data_loaded['n_w'])
        self.n_o = float(data_loaded['n_o'])
        self.krw0 = float(data_loaded['krw0'])
        self.kro0 = float(data_loaded['kro0'])
        self.ds = float(data_loaded['ds'])

    def _stemp(self, S):

        # inds = np.arange(len(S))
        #
        # ite = inds[S > 1 - self.Sor + self.ds]
        # ite2 = inds[S < self.Swc - self.ds]
        #
        # if len(ite) > 0 or len(ite2) > 0:
        #     raise ValueError('Saturation Not supported')
        #
        # S[ite] = 1 - self.Sor
        # S[ite2] = self.Swc
        # return (S - self.Swc) / (1 - self.Swc - self.Sor)

        return S

    def _krw(self, S_temp):
        # return self.krw0*(np.power(S_temp, self.n_w))
        return np.power(S_temp, self.n_w)

    def _kro(self, S_temp):
        # return self.kro0*(np.power(1 - S_temp, self.n_o))
        return np.power(1 - S_temp, self.n_o)

    def calculate(self, saturations):
        n = len(saturations)
        ids = np.arange(n)
        ids_fora = ids[(saturations < 0) | (saturations > 1)]
        if len(ids_fora) > 0:
            raise ValueError('saturacao errada')

        krw = np.zeros(n)
        kro = krw.copy()
        ids1 = ids[saturations > (1 - self.Sor)]
        ids2 = ids[saturations < self.Swc]
        ids_r = np.setdiff1d(ids, np.union1d(ids1, ids2))

        krw[ids1] = np.repeat(self.krw0, len(ids1))
        kro[ids2] = np.repeat(self.kro0, len(ids2))
        Stemp = self._stemp(saturations[ids_r])
        krw[ids_r] = self._krw(Stemp)
        kro[ids_r] = self._kro(Stemp)

        return krw, kro

    def __call__(self, saturations):
        return self.calculate(saturations)
