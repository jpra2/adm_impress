from ..monophasic.monophasic1 import Monophasic
from . import directories_biph
from ...preprocess.preprocess1 import set_saturation_regions
from ... import directories as direc
from . import relative_permeability
import numpy as np

class Biphasic(Monophasic):

    def __init__(self, M, load=False):
        super().__init__(M)
        name_relative_permeability = direc.data_loaded['biphasic_data']['relative_permeability']
        self.relative_permeability = getattr(relative_permeability, name_relative_permeability)
        self.name_datas = directories_biph.name_datas
        self.mi_w = direc.data_loaded['biphasic_data']['mi_w']
        self.mi_o = direc.data_loaded['biphasic_data']['mi_o']

        if not load:
            set_saturation_regions(M)
            self.set_gama()
            self.update_relative_permeability()
            self.update_mobilities()
            M.data.update_variables_to_mesh()
            M.data.export_variables_to_npz()
        else:
            self.load_gama()

    def set_gama(self):
        M = self.mesh

        gama_w = direc.data_loaded['biphasic_data']['gama_w']
        gama_o = direc.data_loaded['biphasic_data']['gama_o']

        gama_w = np.repeat(gama_w, self.n_volumes)
        gama_o = np.repeat(gama_o, self.n_volumes)

        saturations = M.data.get_variable('saturation')

        self.gama = gama_w*saturations + gama_o*(1-saturations)
        M.data.set_variable('gama', self.gama)

    def load_gama(self):
        self.gama = self.mesh.data.get_variable('gama')

    def load_infos_biphasic(self):

        pass

    def update_relative_permeability(self):
        M = self.mesh
        krw, kro = self.relative_permeability(M.data.get_variable('saturation'))
        M.data.set_variable('krw', krw)
        M.data.set_variable('kro', kro)

    def update_mobilities(self):
        M = self.mesh
        n = self.n_volumes
        lambda_w = M.data.get_variable('krw')/self.mi_w
        lambda_o = M.data.get_variable('kro')/self.mi_o
        lambda_t = lambda_w + lambda_o
        fw_vol = lambda_w/lambda_t

        M.data.set_variable('lambda_w', lambda_w)
        M.data.set_variable('lambda_o', lambda_o)
        M.data.set_variable('lambda_t', lambda_t)
        M.data.set_variable('fw_vol', fw_vol)
