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
            self.update_transmissibility_ini()
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

        saturations = M.data['saturation']

        self.gama = gama_w*saturations + gama_o*(1-saturations)
        M.data['gama'] = self.gama

    def load_gama(self):
        self.gama = self.mesh.data['gama']

    def load_infos_biphasic(self):

        pass

    def update_relative_permeability(self):
        M = self.mesh
        krw, kro = self.relative_permeability(M.data['saturation'])
        M.data['krw'] = krw
        M.data['kro'] = kro

    def update_mobilities(self):
        M = self.mesh
        n = self.n_volumes
        lambda_w = M.data['krw']/self.mi_w
        lambda_o = M.data['kro']/self.mi_o
        lambda_t = lambda_w + lambda_o
        fw_vol = lambda_w/lambda_t

        M.data['lambda_w'] = lambda_w
        M.data['lambda_o'] = lambda_o
        M.data['lambda_t'] = lambda_t
        M.data['fw_vol'] = fw_vol

    def update_transmissibility_ini(self):
        M = self.mesh

        pretransmissibility = M.data['pretransmissibility'].copy()
        internal_faces = M.data.elements_lv0['internal_faces']
        b_faces = M.data.elements_lv0['boundary_faces']
        t_internal = pretransmissibility[internal_faces]
        vols_viz_internal_faces = M.data.elements_lv0['vols_viz_internal_faces']
        vols_viz_boundary_faces = M.data.elements_lv0['vols_viz_boundary_faces'].flatten()
        fw_vol = M.data['fw_vol']
        lamda_t = M.data['lambda_t']

        ws_inj = M.contours.datas['ws_inj']

        v0 = vols_viz_internal_faces[:, 0]
        v1 = vols_viz_internal_faces[:, 1]
        ids = np.arange(len(v0))

        idv0 = np.array([], dtype=np.int64)
        idv1 = idv0.copy()
        for w in ws_inj:
            vv0 = ids[v0 == w]
            vv1 = ids[v1 == w]
            idv0 = np.append(idv0, vv0)
            idv1 = np.append(idv1, vv1)

        idv0 = np.array(idv0).flatten()
        idv1 = np.array(idv1).flatten()
        idsv = np.union1d(idv0, idv1)
        ids_fora = np.setdiff1d(ids, idsv)

        total_mobility_internal_faces = np.zeros(len(internal_faces))
        fw_internal_faces = total_mobility_internal_faces.copy()

        total_mobility_internal_faces[idv0] = lambda_t[v0[idv0]]
        total_mobility_internal_faces[idv1] = lambda_t[v1[idv1]]
        total_mobility_internal_faces[ids_fora] = (lambda_t[v0[ids_fora]] + lambda_t[v1[ids_fora]])/2

        fw_internal_faces[idv0] = fw_vol[v0[idv0]]
        fw_internal_faces[idv1] = fw_vol[v1[idv1]]
        fw_internal_faces[ids_fora] = (fw_vol[v0[ids_fora]] + fw_vol[v1[ids_fora]])/2

        total_mobility_faces = transmissibility.copy()
        fw_faces = transmissibility.copy()

        total_mobility_faces[internal_faces] = total_mobility_internal_faces
        total_mobility_faces[b_faces] = lamda_t[vols_viz_boundary_faces]

        fw_faces[internal_faces] = fw_internal_faces
        fw_faces[b_faces] = fw_vol[vols_viz_boundary_faces]

        transmissibility = pretransmissibility*total_mobility_faces

        M.data['transmissibility'] = transmissibility
        M.data['lambda_t_faces'] = total_mobility_faces
        M.data['fw_faces'] = fw_faces

















        import pdb; pdb.set_trace()
