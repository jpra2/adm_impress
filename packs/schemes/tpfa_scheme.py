
from ..data_class.data_manage import dataManager
from .. import directories as direc
import numpy as np
import scipy.sparse as sp

class tpfaScheme:

    def __init__(self, M: 'mesh object', data_name: 'nome do arquivo para ser gravado') -> None:

        self.mesh = M
        self.gravity = direc.data_loaded['gravity']
        self.n_nodes = len(M.data.elements_lv0[direc.entities_lv0[0]])
        self.n_faces = len(M.data.elements_lv0[direc.entities_lv0[1]])
        self.n_edges = len(M.data.elements_lv0[direc.entities_lv0[2]])
        self.n_volumes = len(M.data.elements_lv0[direc.entities_lv0[3]])
        self.data = dataManager(data_name)
        M.scheme = self

    def export_data(self):
        self.data.export_to_npz()

    def get_flux_faces_and_volumes(self) -> None:

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
        velocity_faces = M.data.variables[M.data.variables_impress['velocity_faces']]
        u_normal = M.data[M.data.variables_impress['u_normal']]

        x = self.data['p']

        ps0 = x[v0[:, 0]]
        ps1 = x[v0[:, 1]]

        flux_internal_faces = -((ps1 - ps0) * t0 + M.data.variables[M.data.variables_impress['flux_grav_faces']][internal_faces])
        velocity = (flux_internal_faces / a0).reshape([len(internal_faces), 1])
        velocity = velocity * u_normal[internal_faces]
        velocity_faces[internal_faces] = velocity
        flux_faces = M.data.variables[M.data.variables_impress['flux_faces']]

        flux_faces[internal_faces] = flux_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_internal_faces, -flux_internal_faces]).flatten()
        flux_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        M.data.variables[M.data.variables_impress['flux_volumes']] = flux_volumes
        M.data.variables[M.data.variables_impress['flux_faces']] = flux_faces
        M.data.variables[M.data.variables_impress['velocity_faces']] = velocity_faces
        # M.data.update_variables_to_mesh(['flux_volumes', 'flux_faces', 'velocicty_faces'])

    def get_gravity_source_term(self):
        M = self.mesh
        centroids = M.data['centroid_volumes']
        source_term_volumes = np.zeros(len(centroids))
        transmissibility_faces = M.data.variables[M.data.variables_impress['transmissibility']]
        source_term_faces = np.zeros(len(transmissibility_faces))
        internal_faces = M.data.elements_lv0[direc.entities_lv0_0[0]]

        if self.gravity:

            gamma = M.data['gama']

            vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[2]]
            v0 = vols_viz_internal_faces
            transmissibility_faces = M.data.variables[M.data.variables_impress['transmissibility']]
            transmissibility_internal_faces = transmissibility_faces[internal_faces]
            t0 = transmissibility_internal_faces
            zs = centroids[:, 2]

            # for i, f in enumerate(internal_faces):
            #     q_g_f = -1*(zs[v0[i][1]]*gamma[v0[i][1]] - zs[v0[i][0]]*gamma[v0[i][0]])*t0[i]
            #     source_term_faces[f] = q_g_f
            #     source_term_volumes[v0[i][0]] += q_g_f
            #     source_term_volumes[v0[i][1]] -= q_g_f
            #
            # import pdb; pdb.set_trace()

            source_term_internal_faces = -1*(zs[v0[:, 1]]*gamma[v0[:, 1]] - zs[v0[:, 0]]*gamma[v0[:, 0]])*t0
            source_term_faces[internal_faces] = source_term_internal_faces

            lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
            cols = np.zeros(len(lines), dtype=np.int32)
            data = np.array([source_term_internal_faces, -source_term_internal_faces]).flatten()
            source_term_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        M.data.variables[M.data.variables_impress['flux_grav_volumes']] = source_term_volumes
        M.data.variables[M.data.variables_impress['flux_grav_faces']] = source_term_faces

    def get_RHS_term(self):
        M = self.mesh
        wells = M.contours.datas
        b = np.zeros(self.n_volumes)
        self.get_gravity_source_term()
        b[:] = M.data.variables[M.data.variables_impress['flux_grav_volumes']]

        if self.gravity:

            M.contours.add_gravity(M, M.data.variables[M.data.variables_impress['gama']])

        ws_p = wells['ws_p']  # pocos de pressao prescrita
        values_p = wells['values_p']  # valores de pressao prescrita
        ws_q = wells['ws_q']  # pocos de vazao prescrita
        values_q = wells['values_q']  # valor da vazao prescrita

        b[ws_p] = values_p
        b[ws_q] += values_q

        self.data['b'] = b

    def get_transmissibility_matrix(self):
        '''
        transmissibility matrix with boundary conditions
        '''

        M = self.mesh
        contours = M.contours.datas

        ws_p = contours['ws_p']

        T = self.data['Tini'].tolil().copy()
        T[ws_p] = np.zeros((len(ws_p), T.shape[0]))
        T[ws_p, ws_p] = np.ones(len(ws_p))

        self.data['T'] = T.tocsc()

    def get_transmissibility_matrix_without_boundary_conditions(self) -> None:
        '''
        transmissibility matrix without boundary conditions
        '''

        M = self.mesh
        vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[2]]
        v0 = vols_viz_internal_faces
        internal_faces = M.data.elements_lv0[direc.entities_lv0_0[0]]
        transmissibility_faces = M.data.variables[M.data.variables_impress['transmissibility']]
        transmissibility_internal_faces = transmissibility_faces[internal_faces]
        t0 = transmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, self.n_volumes))

        self.data['Tini'] = T

    def load_data(self):
        self.data.load_from_npz()

    def update_initial_pressure_guess(self, p) -> None:
        '''
        update initial pressure guess
        '''
        M = self.mesh
        self.data['p'] = p
        M.data.variables[M.data.variables_impress['pressure']] = p
