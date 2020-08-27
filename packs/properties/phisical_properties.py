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

        self._gravity = data_loaded['gravity']
        self._gravity_vector = np.array(data_loaded['gravity_vector'], dtype=float)

    @property
    def gravity_vector(self):
        return self._gravity_vector

    def get_nkga(self, keq_faces, u_normal_faces, areas):

        assert u_normal_faces.shape[1] == 3
        assert len(keq_faces) == len(u_normal_faces) == len(areas)

        nkg = np.empty(len(keq_faces))

        # u2 = u_normal_faces.copy()
        #
        # pdb.set_trace()
        # yy = np.tile(self.gravity_vector,(u_normal_faces.shape[0], 1))
        # # yy = yy.reshape((1, u_normal_faces.shape[0], 3))
        #
        # vv = np.tensordot(u2, yy, axes=(0,1))
        #
        # bb = np.tensordot(u_normal_faces, yy, axes=1)


        for i in range(len(keq_faces)):

            nkg[i] = np.dot(u_normal_faces[i], self.gravity_vector)*keq_faces[i]*areas[i]

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
