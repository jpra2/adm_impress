import directories as direc
import scipy.sparse as sp

class Monophasic:

    def __init__(self, M):
        self.mesh = M
        self.n_nodes = len(M.data.elements_lv0[drec.elements_lv0[0]])
        self.n_faces = len(M.data.elements_lv0[drec.elements_lv0[1]])
        self.n_edges = len(M.data.elements_lv0[drec.elements_lv0[2]])
        self.n_volumes = len(M.data.elements_lv0[drec.elements_lv0[3]])

    def get_transmissibility_matrix(self):
        M = self.mesh
        vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[1]]
        v0 = vols_viz_internal_faces
        internal_faces = M.data.elements_lv0[direc.entities_lv0_0[0]]
        transmissibility_faces = M.data.variables[M.data.variables_impress['transmissibility']]
        transmissibility_internal_faces = transmissibility_faces[internal_faces]
        t0 = transmissibility_internal_faces

        # csc_matrix((data,(lines,cols)),shape=(len(M1.all_volumes),len(M1.all_volumes)))

        lines = np.array([v0[:,0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:,1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0[:], t0[:], -t0[:], -t0[:]]).flatten()

        T = sp.csc_matrix((data,(lines,cols)),shape=(self.n_volumes,self.n_volumes))

        self.T = T
