import yaml
import numpy as np
import pdb
import scipy.sparse as sp

class PhisicalProperties:

    def __init__(self, input_file=''):
        if input_file == '':
            input_file = 'input_cards/physic/physic_card.yml'

        with open(input_file, 'r') as f:
            data_loaded = yaml.safe_load(f)

        # self._gravity = data_loaded['gravity']
        self._gravity_vector = np.array(data_loaded['gravity_vector'], dtype=float)
        PhisicalProperties._gravity_vector = np.array(data_loaded['gravity_vector'], dtype=float)
        soma = np.sum(np.absolute(self._gravity_vector))
        if soma > 0:
            self._gravity = True
        else:
            self._gravity = False

    @property
    def gravity_vector(self):
        return self._gravity_vector

    def get_nkga(self, keq_faces, u_normal_faces, areas):

        assert u_normal_faces.shape[1] == 3
        assert len(keq_faces) == len(u_normal_faces) == len(areas)
        ni = len(keq_faces)

        nkg = u_normal_faces*self.gravity_vector
        nkg = nkg.sum(axis=1)
        nkg = nkg*keq_faces
        nkg = nkg*areas

        return nkg

    def get_g_source_w_o_internal_faces(self, nkga_internal_faces, mob_w_internal_faces, mob_o_internal_faces, rho_w, rho_o):

        assert len(nkga_internal_faces) ==  len(mob_w_internal_faces) == len(mob_o_internal_faces)

        g_source_w_internal_faces = nkga_internal_faces*mob_w_internal_faces*rho_w
        g_source_o_internal_faces = nkga_internal_faces*mob_o_internal_faces*rho_o

        return g_source_w_internal_faces, g_source_o_internal_faces

    @staticmethod
    def get_total_g_source_volumes(volumes, volumes_adj_internal_faces, g_source_total_internal_faces):

        v0 = volumes_adj_internal_faces
        n_volumes = len(volumes)
        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([g_source_total_internal_faces, -g_source_total_internal_faces]).flatten()
        flux_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
        return flux_volumes
