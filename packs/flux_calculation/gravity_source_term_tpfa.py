import numpy as np
import scipy.sparse as sp

class GravitySourceTermTpfa:

    def get_gravity_source_term(self):

        centroids = self.data_impress['centroid_volumes']
        source_term_volumes = np.zeros(len(centroids))
        transmissibility_faces = self.data_impress[self.data_impress.variables_impress['transmissibility']]
        source_term_faces = np.zeros(len(transmissibility_faces))
        internal_faces = self.elements_lv0['internal_faces']

        if self.gravity:

            gamma = self.data_impress['gama']

            vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
            v0 = vols_viz_internal_faces
            transmissibility_internal_faces = transmissibility_faces[internal_faces]
            t0 = transmissibility_internal_faces
            zs = centroids[:, 2]

            source_term_internal_faces = -1*(zs[v0[:, 1]]*gamma[v0[:, 1]] - zs[v0[:, 0]]*gamma[v0[:, 0]])*t0
            source_term_faces[internal_faces] = source_term_internal_faces

            lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
            cols = np.zeros(len(lines), dtype=np.int32)
            data = np.array([source_term_internal_faces, -source_term_internal_faces]).flatten()
            source_term_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        self.data_impress[self.data_impress.variables_impress['flux_grav_volumes']] = source_term_volumes.copy()
        self.data_impress[self.data_impress.variables_impress['flux_grav_faces']] = source_term_faces.copy()
