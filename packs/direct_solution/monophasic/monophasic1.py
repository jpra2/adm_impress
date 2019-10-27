from packs import directories as direc
import scipy.sparse as sp
import numpy as np

class Monophasic:

    def __init__(self, M):
        print('\n Monophasic object \n')
        self.mesh = M
        self.gravity = direc.data_loaded['gravity']
        self.n_nodes = len(M.data.elements_lv0[direc.entities_lv0[0]])
        self.n_faces = len(M.data.elements_lv0[direc.entities_lv0[1]])
        self.n_edges = len(M.data.elements_lv0[direc.entities_lv0[2]])
        self.n_volumes = len(M.data.elements_lv0[direc.entities_lv0[3]])
        self.datas = dict()

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

        self.datas['Tini'] = T

    def get_RHS_term(self):
        M = self.mesh
        contours = M.contours.datas

        b = np.zeros(self.n_volumes)

        ws_p = contours['ws_p'] # pocos de pressao prescrita
        values_p = contours['values_p'] # valores de pressao prescrita
        ws_q = contours['ws_q'] # pocos de vazao prescrita
        values_q = contours['values_q'] # valor da vazao prescrita

        b[ws_p] = values_p
        b[ws_q] = values_q

        self.datas['b'] = b

    def get_transmissibility_matrix(self):

        M = self.mesh
        contours = M.contours.datas

        ws_p = contours['ws_p'].flatten()

        T = self.datas['Tini'].tolil().copy()
        T[ws_p] = sp.lil_matrix((len(ws_p), T.shape[0]))
        T[ws_p, ws_p] = np.ones(len(ws_p))

        self.datas['T'] = T

    def get_solution(self, x):
        M = self.mesh
        self.datas['x'] = x
        M.data.variables[M.data.variables_impress['pressure']] = x
        # M.data.update_variables_to_mesh(['pressure'])

    def get_flux_faces_and_volumes(self):

        M = self.mesh

        vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[2]]
        v0 = vols_viz_internal_faces
        internal_faces = M.data.elements_lv0[direc.entities_lv0_0[0]]
        transmissibility_faces = M.data.variables[M.data.variables_impress['transmissibility']]
        transmissibility_internal_faces = transmissibility_faces[internal_faces]
        t0 = transmissibility_internal_faces
        area_faces = M.data.variables[M.data.variables_impress['area']]
        area_internal_faces = area_faces[internal_faces]
        a0 = area_internal_faces
        velocicty_faces = M.data.variables[M.data.variables_impress['velocity_faces']]

        x = self.datas['x']

        ps0 = x[v0[:,0]]
        ps1 = x[v0[:,1]]

        flux_internal_faces = (ps1-ps0)*t0 + self.get_gravity_source_term()
        velocity = flux_internal_faces/a0
        velocicty_faces[internal_faces] = velocity
        flux_volumes = M.data.variables[M.data.variables_impress['flux_volumes']]
        flux_faces = M.data.variables[M.data.variables_impress['flux_faces']]

        flux_faces[internal_faces] = flux_internal_faces

        lines = np.array([v0[:,0], v0[:,1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_internal_faces, -flux_internal_faces]).flatten()
        flux_volumes = sp.csc_matrix((data,(lines,cols)),shape=(self.n_volumes,1)).toarray()

        M.data.variables[M.data.variables_impress['flux_volumes']] = flux_volumes
        M.data.variables[M.data.variables_impress['flux_faces']] = flux_faces
        M.data.variables[M.data.variables_impress['velocicty_faces']] = velocity_faces
        # M.data.update_variables_to_mesh(['flux_volumes', 'flux_faces', 'velocicty_faces'])

    def get_gravity_source_term(self):
        M = self.mesh
        centroids = M.data.centroids[direc.entities_lv0[3]]
        source_term_volumes = np.zeros(len(centroids))
        transmissibility_faces = M.data.variables[M.data.variables_impress['transmissibility']]
        internal_faces = M.data.elements_lv0[direc.entities_lv0_0[0]]

        if self.gravity:

            gamma = direc.data_loaded['monophasic_data']['gama']

            vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[2]]
            v0 = vols_viz_internal_faces
            transmissibility_faces = M.data.variables[M.data.variables_impress['transmissibility']]
            transmissibility_internal_faces = transmissibility_faces[internal_faces]
            t0 = transmissibility_internal_faces
            zs = centroids[:,2]
            delta_zs = zs[v0[:,1]] - zs[v0[:,0]] # tamanho de internal faces, delta_z para cada face
            source_term_internal_faces = delta_zs*gamma

            source_term_faces = np.zeros(len(transmissibility_faces))
            source_term_faces[internal_faces] = source_term_internal_faces


            source_term_volumes[v0[:,0]] += source_term_internal_faces
            source_term_volumes[v0[:,1]] -= source_term_internal_faces

        else:
            source_term_volumes = np.zeros(len(centroids))
            source_term_faces = np.zeros(len(transmissibility_faces))

        M.data.variables[M.data.variables_impress['flux_grav_volumes']] = source_term_volumes
        M.data.variables[M.data.variables_impress['flux_grav_faces']] = source_term_faces

        return source_term_volumes
