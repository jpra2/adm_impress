import numpy as np
import scipy.sparse as sp
from ..directories import data_loaded
import pdb

class TpfaFlux:

    def get_gravity_source_term(self):

        biphasic = data_loaded['biphasic']
        monophasic = data_loaded['monophasic']

        if biphasic == monophasic:
            raise ValueError('biphasic == monophasic')

        if monophasic:
            self.get_gravity_source_term_mono()
            return 0

        centroids = self.data_impress['centroid_volumes']
        source_term_volumes = np.zeros(len(centroids))
        transmissibility_faces = self.data_impress[self.data_impress.variables_impress['transmissibility']]
        source_term_faces = np.zeros(len(transmissibility_faces))
        internal_faces = self.elements_lv0['internal_faces']
        areas_internal_faces = self.data_impress['area'][internal_faces]
        k_harm_internal_faces = self.data_impress['k_harm'][internal_faces]
        dh_internal_faces = self.data_impress['dist_cent'][internal_faces]
        u_normal_internal_faces = self.data_impress['u_normal'][internal_faces]
        self._data['grav_source_term_water_volumes'] = source_term_volumes.copy()
        self._data['grav_source_term_water_faces'] = source_term_faces.copy()

        # upwind_grav = np.full((len(internal_faces), 2), False, dtype=bool)

        # import pdb; pdb.set_trace()

        if self.gravity:

            vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
            v0 = vols_viz_internal_faces
            source_term_internal_faces = self.flux_grav_total_faces[internal_faces]
            # source_term_internal_faces = np.rint(source_term_internal_faces)
            source_term_faces[internal_faces] = source_term_internal_faces

            # grav_source_term_water_internal_faces_2 = -1*(zs[v0[:, 1]] - zs[v0[:, 0]])*(self.lambda_w_internal_faces*gama_w)*areas_internal_faces*k_harm_internal_faces/dh_internal_faces
            # grav_source_term_water_internal_faces = -1*(self.grad_z_internal_faces)*(lambda_w_internal_faces*gama_w)*areas_internal_faces*k_harm_internal_faces
            grav_source_term_water_internal_faces = self.flux_grav_w_faces[internal_faces]
            # grav_source_term_water_internal_faces = 1*(zs[v0[:, 1]] - zs[v0[:, 0]])*(lambda_w_internal_faces*gama_w)*areas_internal_faces*k_harm_internal_faces/dh_internal_faces



            lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
            cols = np.zeros(len(lines), dtype=np.int32)
            data = np.array([source_term_internal_faces, -source_term_internal_faces]).flatten()
            source_term_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

            # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
            # cols = np.zeros(len(lines), dtype=np.int32)
            data = np.array([grav_source_term_water_internal_faces, -grav_source_term_water_internal_faces]).flatten()
            self._data['grav_source_term_water_volumes'] = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
            self._data['grav_source_term_water_faces'][internal_faces] = grav_source_term_water_internal_faces

        self.data_impress[self.data_impress.variables_impress['flux_grav_volumes']] = source_term_volumes.copy()
        self.data_impress[self.data_impress.variables_impress['flux_grav_faces']] = source_term_faces.copy()

        self.data_impress['flux_grav_w_faces_vec'][internal_faces] = self._data['grav_source_term_water_faces'][internal_faces].reshape(len(internal_faces), 1)*u_normal_internal_faces
        self.data_impress['flux_grav_o_faces_vec'][internal_faces] = (self.data_impress['flux_grav_faces'][internal_faces] - self._data['grav_source_term_water_faces'][internal_faces]).reshape(len(internal_faces), 1)*u_normal_internal_faces

        return 0

    def get_gravity_source_term_mono(self):

        centroids = self.data_impress['centroid_volumes']
        source_term_volumes = np.zeros(len(centroids))
        transmissibility_faces = self.data_impress[self.data_impress.variables_impress['transmissibility']]
        source_term_faces = np.zeros(len(transmissibility_faces))
        internal_faces = self.elements_lv0['internal_faces']
        areas_internal_faces = self.data_impress['area'][internal_faces]
        k_harm_internal_faces = self.data_impress['k_harm'][internal_faces]
        dh_internal_faces = self.data_impress['dist_cent'][internal_faces]
        u_normal_internal_faces = self.data_impress['u_normal'][internal_faces]
        self._data['grav_source_term_water_volumes'] = source_term_volumes.copy()
        self._data['grav_source_term_water_faces'] = source_term_faces.copy()

        # upwind_grav = np.full((len(internal_faces), 2), False, dtype=bool)

        # import pdb; pdb.set_trace()

        if self.gravity:

            gamma = self.data_impress['gama']
            gama_faces = self.data_impress['gama_faces']

            vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
            v0 = vols_viz_internal_faces
            transmissibility_internal_faces = transmissibility_faces[internal_faces]
            t0 = transmissibility_internal_faces
            zs = centroids[:, 2]

            source_term_internal_faces = -1*(zs[v0[:, 1]] - zs[v0[:, 0]])*t0*gama_faces[internal_faces]

            lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
            cols = np.zeros(len(lines), dtype=np.int32)
            data = np.array([source_term_internal_faces, -source_term_internal_faces]).flatten()
            source_term_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        self.data_impress[self.data_impress.variables_impress['flux_grav_volumes']] = source_term_volumes.copy()
        self.data_impress[self.data_impress.variables_impress['flux_grav_faces']] = source_term_faces.copy()

    def get_flux_faces_and_volumes(self) -> None:

        self.data_impress['flux_volumes'] = self.flux_volumes
        self.data_impress['flux_faces'] = self.flux_faces
        self.data_impress['velocity_faces'] = self.velocity_faces_vec

    def get_velocity_faces(self):

        internal_faces = self.elements_lv0['internal_faces']
        flux_internal_faces = self.data_impress['flux_faces'][internal_faces]
        areas_internal_faces = self.data_impress['area'][internal_faces]
        u_normal_internal_faces = self.data_impress['u_normal'][internal_faces]
        area_faces = self.data_impress['area']
        area_internal_faces = area_faces[internal_faces]
        a0 = area_internal_faces
        velocity_faces = np.zeros(self.data_impress['velocity_faces'].shape)

        velocity = (flux_internal_faces / a0).reshape([len(internal_faces), 1])
        velocity = velocity * u_normal_internal_faces
        velocity_faces[internal_faces] = velocity
        self.data_impress['velocity_faces'] = velocity_faces

    def get_flux_volumes(self):

        vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        v0 = vols_viz_internal_faces
        internal_faces = self.elements_lv0['internal_faces']
        flux_internal_faces = self.data_impress['flux_faces'][internal_faces]


        # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        lines=np.concatenate(v0.T)
        # cols = np.repeat(0, len(lines))
        # data = np.array([flux_internal_faces, -flux_internal_faces]).flatten()
        data=np.concatenate([flux_internal_faces, -flux_internal_faces])
        # flux_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        flux_volumes=np.bincount(lines,weights=data)

        self.data_impress['flux_volumes'] = flux_volumes

class TpfaFlux2:

    def get_flux_volumes_and_velocity(self) -> None:

        vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        v0 = vols_viz_internal_faces
        n_volumes = len(self.elements_lv0['volumes'])
        internal_faces = self.elements_lv0['internal_faces']
        area_faces = self.data_impress['area']
        area_internal_faces = area_faces[internal_faces]
        a0 = area_internal_faces
        velocity_faces = np.zeros(self.data_impress['velocity_faces'].shape)
        u_normal = self.data_impress['u_normal']
        flux_faces = self.data_impress['flux_faces']

        flux_internal_faces = flux_faces[internal_faces]
        velocity = (flux_internal_faces / a0).reshape([len(internal_faces), 1])
        velocity = velocity * u_normal[internal_faces]
        velocity_faces[internal_faces] = velocity

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_internal_faces, -flux_internal_faces]).flatten()
        flux_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()

        self.data_impress['flux_volumes'] = flux_volumes
        self.data_impress['velocity_faces'] = velocity_faces

def get_flux_faces(p0: 'left side pressure',
    p1: 'right side pressure',
    t0: 'transmissibility faces',
    flux_grav_faces: 'gravity flux of faces'
):

    return -((p1 - p0) * t0 - flux_grav_faces)
