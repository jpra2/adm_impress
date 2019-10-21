import directories as direc
import scipy.sparse as sp
import numpy as np

class Monophasic:

    def __init__(self, M):
        self.mesh = M
        self.n_nodes = len(M.data.elements_lv0[direc.entities_lv0[0]])
        self.n_faces = len(M.data.elements_lv0[direc.entities_lv0[1]])
        self.n_edges = len(M.data.elements_lv0[direc.entities_lv0[2]])
        self.n_volumes = len(M.data.elements_lv0[direc.entities_lv0[3]])
        self.Tini = None

    def get_transmissibility_matrix_without_contours(self):
        M = self.mesh
        vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[2]]
        v0 = vols_viz_internal_faces
        internal_faces = M.data.elements_lv0[direc.entities_lv0_0[0]]
        transmissibility_faces = M.data.variables[M.data.variables_impress['transmissibility']]
        transmissibility_internal_faces = transmissibility_faces[internal_faces]
        t0 = transmissibility_internal_faces

        lines = np.array([v0[:,0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:,1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0[:], t0[:], -t0[:], -t0[:]]).flatten()

        T = sp.csc_matrix((data,(lines,cols)),shape=(self.n_volumes,self.n_volumes))

        self.Tini = T

    def get_RHS_term(self):
        M = self.mesh
        contours = M.contours.datas

        b = np.zeros(self.n_volumes)

        ws_p = contours['ws_p']
        values_p = contours['values_p']
        ws_q = contours['ws_q']
        values_q = contours['values_q']

        b[ws_p] = values_p
        b[ws_q] = values_q

        self.b = b

    def get_transmissibility_matrix(self):

        M = self.mesh
        contours = M.contours.datas

        if self.Tini:
            pass
        else:
            self.get_transmissibility_matrix_without_contours()

        ws_p = contours['ws_p'].flatten()

        T = self.Tini.tolil(copy=True)
        T[ws_p] = sp.lil_matrix((len(ws_p), T.shape[0]))
        T[ws_p, ws_p] = np.ones(len(ws_p))

        self.T = T

    def get_solution(self, x):
        self.x = x
        M.data.variables[M.data.variables_impress['pressure']] = x
        M.data.update_variables_to_mesh(['pressure'])
