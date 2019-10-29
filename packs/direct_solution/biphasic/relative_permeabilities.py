import pdb
from ... import directories as direc



class UpdatingDatas:

    def __init__(self, M):
        self.Sor = direc.data_loaded['biphasic_data']['Sor']
        self.Swc = direc.data_loaded['biphasic_data']['Swc']
        self.n_w = direc.data_loaded['biphasic_data']['n_w']
        self.n_o = direc.data_loaded['biphasic_data']['n_o']
        self.gama_w = direc.data_loaded['biphasic_data']['gama_w']
        self.gama_o = direc.data_loaded['biphasic_data']['gama_o']
        self.mesh = M
        M.relative_permeability = self

    def brooks_and_corey(self, S):

        if S > (1 - self.Sor) and S <= 1:
            krw = 1.0
            kro = 0.0
        elif S < self.Swc and S >= 0:
            krw = 0.0
            kro = 1.0
        else:
            S_temp = (S - self.Swc) / (1 - self.Swc - self.Sor)
            if S_temp < 0 or S_temp > 1:
                print('erro S_temp')
                pdb.set_trace()
            krw = (S_temp) ** (n_w)
            kro = (1 - S_temp) ** (n_o)

        return krw, kro