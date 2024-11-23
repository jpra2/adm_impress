import numpy as np
import scipy.sparse as sp
from ..data_class.structured_mesh_properties import StructuredMeshProperties
from ..directories import data_loaded
import pdb


class biphasicProperties(StructuredMeshProperties):

    def gama_volumes_average():
        doc = "The gama_volumes_average property."
        def fget(self):
            gama_w = self.biphasic_data['gama_w']
            gama_o = self.biphasic_data['gama_o']
            lambda_w = self.lambda_w_volumes
            lambda_o = self.lambda_o_volumes
            return np.repeat(gama_w, len(lambda_w))*lambda_w + np.repeat(gama_o, len(lambda_o))*lambda_o
        return locals()
    gama_volumes_average = property(**gama_volumes_average())

    def gama_w():
        doc = "The gama_w property."
        def fget(self):
            gama_w = self.biphasic_data['gama_w']
            return gama_w
        return locals()
    gama_w = property(**gama_w())

    def gama_o():
        doc = "The gama_o property."
        def fget(self):
            gama_o = self.biphasic_data['gama_o']
        return locals()
    gama_o = property(**gama_o())

    def lambda_w_volumes():
        doc = "The lambda_w_volumes property."
        def fget(self):
            return self.data_impress['krw']/self.biphasic_data['mi_w']
        return locals()
    lambda_w_volumes = property(**lambda_w_volumes())

    def lambda_o_volumes():
        doc = "The lambda_o_volumes property."
        def fget(self):
            return self.data_impress['kro']/self.biphasic_data['mi_o']
        return locals()
    lambda_o_volumes = property(**lambda_o_volumes())

    def lambda_t_volumes():
        doc = "The lambda_t_volumes property."
        def fget(self):
            return self.lambda_w_volumes + self.lambda_o_volumes
        return locals()
    lambda_t_volumes = property(**lambda_t_volumes())

    def fw_volumes():
        doc = "The fw_volumes property."
        def fget(self):
            return self.lambda_w_volumes/self.lambda_t_volumes
        return locals()
    fw_volumes = property(**fw_volumes())

    def lambda_w_internal_faces():
        doc = "The lambda_w_internal_faces property."
        def fget(self):
            # return self.data_impress['lambda_w'][self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            lambda_w_internal_faces = self.lambda_w_volumes[self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            faces_normal_with_z = self.faces_normal_with_z
            if len(faces_normal_with_z) > 0:
                internal_faces = self.elements_lv0['internal_faces']
                faces_with_lambda = np.intersect1d(faces_normal_with_z, internal_faces)
                remaped_faces_with_lambda = self.rmap_internal_faces(faces_with_lambda)
                lambda_w_internal_faces2 = np.zeros(len(internal_faces))
                lambda_w_internal_faces2[remaped_faces_with_lambda] = lambda_w_internal_faces[remaped_faces_with_lambda]
                lambda_w_internal_faces = lambda_w_internal_faces2

            return lambda_w_internal_faces
        return locals()
    lambda_w_internal_faces = property(**lambda_w_internal_faces())

    def lambda_o_internal_faces():
        doc = "The lambda_o_internal_faces property."
        def fget(self):
            # return self.data_impress['lambda_o'][self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            lambda_o_internal_faces = self.lambda_o_volumes[self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate_o']]]
            faces_normal_with_z = self.faces_normal_with_z
            if len(faces_normal_with_z) > 0:
                internal_faces = self.elements_lv0['internal_faces']
                faces_with_lambda = np.intersect1d(faces_normal_with_z, internal_faces)
                remaped_faces_with_lambda = self.rmap_internal_faces(faces_with_lambda)
                lambda_o_internal_faces2 = np.zeros(len(internal_faces))
                lambda_o_internal_faces2[remaped_faces_with_lambda] = lambda_o_internal_faces[remaped_faces_with_lambda]
                lambda_o_internal_faces = lambda_o_internal_faces2

            return lambda_o_internal_faces
        return locals()
    lambda_o_internal_faces = property(**lambda_o_internal_faces())

    def lambda_t_internal_faces():
        doc = "The lambda_t_internal_faces property."
        def fget(self):
            # return self.data_impress['lambda_t'][self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            # return self.lambda_t_volumes[self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            return self.lambda_w_internal_faces + self.lambda_o_internal_faces
        return locals()
    lambda_t_internal_faces = property(**lambda_t_internal_faces())

    def fw_internal_faces():
        doc = "The fw_internal_faces property."
        def fget(self):
            # return self.data_impress['fw_vol'][self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            # return self.fw_volumes[self.elements_lv0['neig_internal_faces'][self._data['upwind_identificate']]]
            return self.lambda_w_internal_faces/self.lambda_t_internal_faces
        return locals()
    fw_internal_faces = property(**fw_internal_faces())

    def lambda_w_faces():
        doc = "The lambda_w_faces property."
        def fget(self):
            b_faces = self.elements_lv0['boundary_faces']
            neig_b_faces = self.elements_lv0['neig_boundary_faces']
            internal_faces = self.elements_lv0['internal_faces']
            n_faces = len(self.elements_lv0['faces'])
            lambda_w_faces = np.zeros(n_faces)
            lambda_w_faces[internal_faces] = self.lambda_w_internal_faces
            lambda_w_faces[b_faces] = self.lambda_w_volumes[neig_b_faces]
            return lambda_w_faces
        return locals()
    lambda_w_faces = property(**lambda_w_faces())

    def lambda_o_faces():
        doc = "The lambda_o_faces property."
        def fget(self):
            b_faces = self.elements_lv0['boundary_faces']
            neig_b_faces = self.elements_lv0['neig_boundary_faces']
            internal_faces = self.elements_lv0['internal_faces']
            n_faces = len(self.elements_lv0['faces'])
            lambda_o_faces = np.zeros(n_faces)
            lambda_o_faces[internal_faces] = self.lambda_o_internal_faces
            lambda_o_faces[b_faces] = self.lambda_o_volumes[neig_b_faces]
            return lambda_o_faces
        return locals()
    lambda_o_faces = property(**lambda_o_faces())

    def lambda_t_faces():
        doc = "The lambda_t_faces property."
        def fget(self):
            return self.lambda_w_faces + self.lambda_o_faces
        return locals()
    lambda_t_faces = property(**lambda_t_faces())

    def fw_faces():
        doc = "The fw_faces property."
        def fget(self):
            return self.lambda_w_faces/self.lambda_t_faces
        return locals()
    fw_faces = property(**fw_faces())

    def flux_press_w_internal_faces():
        doc = "The flux_press_w_internal_faces property."
        def fget(self):
            internal_faces = self.elements_lv0['internal_faces']
            areas_internal_faces = self.data_impress['area'][internal_faces]
            k_harm_internal_faces = self.data_impress['k_harm'][internal_faces]
            flux_press_w_internal_faces = -self.grad_p_internal_faces*areas_internal_faces*k_harm_internal_faces*self.lambda_w_internal_faces

            return flux_press_w_internal_faces
        return locals()
    flux_press_w_internal_faces = property(**flux_press_w_internal_faces())

    def flux_press_o_internal_faces():
        doc = "The flux_press_o_internal_faces property."
        def fget(self):
            internal_faces = self.elements_lv0['internal_faces']
            areas_internal_faces = self.data_impress['area'][internal_faces]
            k_harm_internal_faces = self.data_impress['k_harm'][internal_faces]

            flux_press_o_internal_faces = -self.grad_p_internal_faces*areas_internal_faces*k_harm_internal_faces*self.lambda_o_internal_faces
            return flux_press_o_internal_faces
        return locals()
    flux_press_o_internal_faces = property(**flux_press_o_internal_faces())

    def g_w_internal_faces():
        doc = "The g_w_internal_faces property."
        def fget(self):
            internal_faces = self.elements_lv0['internal_faces']
            if not self.gravity:
                return np.zeros(len(internal_faces))
            gama_w = self.biphasic_data['gama_w']
            # return -(self.grad_z_internal_faces)*gama_w
            return -(self.delta_z_internal_faces)*gama_w
        return locals()
    g_w_internal_faces = property(**g_w_internal_faces())

    def g_o_internal_faces():
        doc = "The g_o_internal_faces property."
        def fget(self):
            internal_faces = self.elements_lv0['internal_faces']
            if not self.gravity:
                return np.zeros(len(internal_faces))
            gama_o = self.biphasic_data['gama_o']
            return -(self.delta_z_internal_faces)*gama_o
        return locals()
    g_o_internal_faces = property(**g_o_internal_faces())

    def flux_grav_w_internal_faces():
        doc = "The flux_grav_w_internal_faces property."
        def fget(self):
            internal_faces = self.elements_lv0['internal_faces']
            if not self.gravity:
                return np.zeros(len(internal_faces))
            areas_internal_faces = self.data_impress['area'][internal_faces]
            k_harm_internal_faces = self.data_impress['k_harm'][internal_faces]
            # up_g = self.up_g
            # lambda_w_internal_faces = self.lambda_w_volumes[up_g]
            lambda_w_internal_faces = self.lambda_w_internal_faces
            # gama_w = self.biphasic_data['gama_w']

            # flux_grav_w_internal_faces = -1*(self.grad_z_internal_faces)*(lambda_w_internal_faces*gama_w)*areas_internal_faces*k_harm_internal_faces
            flux_grav_w_internal_faces = self.g_w_internal_faces*lambda_w_internal_faces*areas_internal_faces*k_harm_internal_faces
            return flux_grav_w_internal_faces
        return locals()
    flux_grav_w_internal_faces = property(**flux_grav_w_internal_faces())

    def flux_grav_o_internal_faces():
        doc = "The flux_grav_o_internal_faces property."
        def fget(self):
            internal_faces = self.elements_lv0['internal_faces']
            if not self.gravity:
                return np.zeros(len(internal_faces))
            areas_internal_faces = self.data_impress['area'][internal_faces]
            k_harm_internal_faces = self.data_impress['k_harm'][internal_faces]
            # up_g = self.up_g
            # lambda_o_internal_faces = self.lambda_o_volumes[up_g]
            lambda_o_internal_faces = self.lambda_o_internal_faces
            # gama_o = self.biphasic_data['gama_o']

            # flux_grav_o_internal_faces = -1*(self.grad_z_internal_faces)*(lambda_o_internal_faces*gama_o)*areas_internal_faces*k_harm_internal_faces
            flux_grav_o_internal_faces = self.g_o_internal_faces*lambda_o_internal_faces*areas_internal_faces*k_harm_internal_faces
            return flux_grav_o_internal_faces
        return locals()
    flux_grav_o_internal_faces = property(**flux_grav_o_internal_faces())

    def flux_grav_total_internal_faces():
        doc = "The flux_grav_total_internal_faces property."
        def fget(self):
            return self.flux_grav_o_internal_faces + self.flux_grav_w_internal_faces
        return locals()
    flux_grav_total_internal_faces = property(**flux_grav_total_internal_faces())

    def flux_grav_w_faces():
        doc = "The flux_grav_w_faces property."
        def fget(self):
            flux_grav_w_faces = np.zeros(len(self.elements_lv0['faces']))
            if self.gravity:
                pass
            else:
                return flux_grav_w_faces
            internal_faces = self.elements_lv0['internal_faces']
            flux_grav_w_faces[internal_faces] = self.flux_grav_w_internal_faces
            return flux_grav_w_faces
        return locals()
    flux_grav_w_faces = property(**flux_grav_w_faces())

    def flux_grav_o_faces():
        doc = "The flux_grav_o_faces property."
        def fget(self):
            flux_grav_o_faces = np.zeros(len(self.elements_lv0['faces']))
            if self.gravity:
                pass
            else:
                return flux_grav_o_faces
            internal_faces = self.elements_lv0['internal_faces']
            flux_grav_o_internal_faces = self.flux_grav_o_internal_faces
            flux_grav_o_faces[internal_faces] = flux_grav_o_internal_faces
            return flux_grav_o_faces
        return locals()
    flux_grav_o_faces = property(**flux_grav_o_faces())

    def flux_grav_total_faces():
        doc = "The flux_grav_total_faces property."
        def fget(self):
            return self.flux_grav_w_faces + self.flux_grav_o_faces
        return locals()
    flux_grav_total_faces = property(**flux_grav_total_faces())

    def flux_w_internal_faces():
        doc = "The flux_w_internal_faces property."
        def fget(self):
            return self.flux_press_w_internal_faces + self.flux_grav_w_internal_faces
        return locals()
    flux_w_internal_faces = property(**flux_w_internal_faces())

    def flux_o_internal_faces():
        doc = "The flux_o_internal_faces property."
        def fget(self):
            return self.flux_press_o_internal_faces + self.flux_grav_o_internal_faces
        return locals()
    flux_o_internal_faces = property(**flux_o_internal_faces())

    @property
    def flux_w_internal_faces2(self):
        internal_faces = self.elements_lv0['internal_faces']
        ak = self.data_impress['k_harm'][internal_faces]*self.data_impress['area'][internal_faces]
        fw_internal_faces = self.fw_internal_faces
        lambda_o_internal_faces = self.lambda_o_internal_faces
        flux_internal_faces = self.flux_internal_faces
        g_o_internal_faces = self.g_o_internal_faces
        g_w_internal_faces = self.g_w_internal_faces

        flux_w_internal_faces = -ak*fw_internal_faces*lambda_o_internal_faces*(g_o_internal_faces - g_w_internal_faces) + fw_internal_faces*flux_internal_faces
        return flux_w_internal_faces

    @property
    def flux_o_internal_faces2(self):
        internal_faces = self.elements_lv0['internal_faces']
        ak = self.data_impress['k_harm'][internal_faces]*self.data_impress['area'][internal_faces]
        fo_internal_faces = (1 - self.fw_internal_faces)
        lambda_w_internal_faces = self.lambda_w_internal_faces
        flux_internal_faces = self.flux_internal_faces
        g_o_internal_faces = self.g_o_internal_faces
        g_w_internal_faces = self.g_w_internal_faces

        flux_o_internal_faces = ak*fo_internal_faces*lambda_w_internal_faces*(g_o_internal_faces - g_w_internal_faces) + fo_internal_faces*flux_internal_faces
        return flux_o_internal_faces

    def flux_internal_faces():
        doc = "The flux_internal_faces property."
        def fget(self):
            return self.flux_w_internal_faces + self.flux_o_internal_faces
        return locals()
    flux_internal_faces = property(**flux_internal_faces())

    def flux_o_faces():
        doc = "The flux_o_faces property."
        def fget(self):
            internal_faces = self.elements_lv0['internal_faces']
            flux_o_faces = np.zeros(len(self.elements_lv0['faces']))
            flux_o_faces[internal_faces] = self.flux_o_internal_faces
            return flux_o_faces
        return locals()
    flux_o_faces = property(**flux_o_faces())

    def flux_w_faces():
        doc = "The flux_w_faces property."
        def fget(self):
            internal_faces = self.elements_lv0['internal_faces']
            flux_w_faces = np.zeros(len(self.elements_lv0['faces']))
            flux_w_faces[internal_faces] = self.flux_w_internal_faces
            return flux_w_faces
        return locals()
    flux_w_faces = property(**flux_w_faces())

    def flux_faces():
        doc = "The flux_faces property."
        def fget(self):
            return self.flux_w_faces + self.flux_o_faces
        return locals()
    flux_faces = property(**flux_faces())

    def flux_volumes():
        doc = "The flux_volumes property."
        def fget(self):
            flux_internal_faces = self.flux_internal_faces
            v0 = self.elements_lv0['neig_internal_faces']
            n_volumes = len(self.data_impress['GID_0'])
            lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
            cols = np.repeat(0, len(lines))
            data = np.array([flux_internal_faces, -flux_internal_faces]).flatten()
            flux_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
            return flux_volumes
        return locals()
    flux_volumes = property(**flux_volumes())

    def flux_w_volumes():
        doc = "The flux_w_volumes property."
        def fget(self):
            v0 = self.elements_lv0['neig_internal_faces']
            n_volumes = len(self.data_impress['GID_0'])
            flux_w_internal_faces = self.flux_w_internal_faces
            lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
            cols = np.repeat(0, len(lines))
            data = np.array([flux_w_internal_faces, -flux_w_internal_faces]).flatten()
            flux_w_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
            ws_prod = self.wells['ws_prod']
            ws_inj = self.wells['ws_inj']
            flux_volumes = self.flux_volumes

            if len(ws_prod) > 0:
                flux_w_volumes[ws_prod] -= flux_volumes[ws_prod]*self.fw_volumes[ws_prod]

            if len(ws_inj) > 0:
                flux_w_volumes[ws_inj] -= flux_volumes[ws_inj]*self.fw_volumes[ws_inj]

            return flux_w_volumes
        return locals()
    flux_w_volumes = property(**flux_w_volumes())

    def flux_o_volumes():
        doc = "The flux_o_volumes property."
        def fget(self):
            v0 = self.elements_lv0['neig_internal_faces']
            n_volumes = len(self.data_impress['GID_0'])
            flux_volumes = self.flux_volumes
            flux_o_internal_faces = self.flux_o_internal_faces
            lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
            cols = np.repeat(0, len(lines))
            data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
            flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
            ws_prod = self.wells['ws_prod']
            # flux_o_volumes[ws_prod] -= flux_volumes[ws_prod]*(1 - self.fw_volumes[ws_prod])
            return flux_o_volumes
        return locals()
    flux_o_volumes = property(**flux_o_volumes())

    @property
    def update_flux_w_volumes_for_wells(self):
        internal_faces = self.elements_lv0['internal_faces']
        flux_volumes = self.data_impress['flux_volumes'].copy()
        flux_w_internal_faces = self.data_impress['flux_w_faces'][internal_faces]
        v0 = self.elements_lv0['neig_internal_faces']
        n_volumes = len(self.data_impress['GID_0'])
        ws_prod = self.wells['ws_prod']
        ws_inj = self.wells['ws_inj']

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_w_internal_faces, -flux_w_internal_faces]).flatten()
        flux_w_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        if self.pare1 == True:
            pdb.set_trace()
            self.pare1 = False

        flux_w_volumes[ws_prod] -= flux_volumes[ws_prod]*self.fw_volumes[ws_prod]
        flux_w_volumes[ws_inj] -= flux_volumes[ws_inj]*self.fw_volumes[ws_inj]
        return flux_w_volumes

    @property
    def update_flux_o_volumes_for_wells(self):
        internal_faces = self.elements_lv0['internal_faces']
        flux_volumes = self.data_impress['flux_volumes'].copy()
        flux_o_internal_faces = self.data_impress['flux_o_faces'][internal_faces]
        v0 = self.elements_lv0['neig_internal_faces']
        n_volumes = len(self.data_impress['GID_0'])
        ws_prod = self.wells['ws_prod']
        ws_inj = self.wells['ws_inj']

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
        flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        flux_o_volumes[ws_prod] -= flux_volumes[ws_prod]*(1 - self.fw_volumes[ws_prod])
        flux_o_volumes[ws_inj] -= flux_volumes[ws_inj]*(1 - self.fw_volumes[ws_inj])
        return flux_o_volumes

    @property
    def update_saturation_wells(self):
        wells = self.wells['all_wells']
        flux_w_volumes_new_wells = self.update_flux_w_volumes_for_wells[wells]

        sats = self.data_impress['saturation_last'][wells]

        fw_volumes = -flux_w_volumes_new_wells
        volumes = self.data_impress['volume'][wells]
        phis = self.data_impress['poro'][wells]

        sats += (fw_volumes*self.delta_t)/(phis*volumes)

        return sats

    def velocity_w_faces():
        doc = "The velocity_w_faces property."
        def fget(self):
            areas = self.data_impress['area']
            return self.flux_w_faces/areas
        return locals()
    velocity_w_faces = property(**velocity_w_faces())

    def velocity_o_faces():
        doc = "The velocity_o_faces property."
        def fget(self):
            areas = self.data_impress['area']
            return self.flux_o_faces/areas
        return locals()
    velocity_o_faces = property(**velocity_o_faces())

    def velocity_faces():
        doc = "The velocity_faces property."
        def fget(self):
            return self.velocity_w_faces + self.velocity_o_faces
        return locals()
    velocity_faces = property(**velocity_faces())

    def velocity_w_faces_vec():
        doc = "The velocity_w_faces_vec property."
        def fget(self):
            u_normal = self.data_impress['u_normal']
            return self.velocity_w_faces.reshape([len(self.elements_lv0['faces']), 1]) * u_normal
        return locals()
    velocity_w_faces_vec = property(**velocity_w_faces_vec())

    def velocity_o_faces_vec():
        doc = "The velocity_o_faces_vec property."
        def fget(self):
            u_normal = self.data_impress['u_normal']
            return self.velocity_o_faces.reshape([len(self.elements_lv0['faces']), 1]) * u_normal
        return locals()
    velocity_o_faces_vec = property(**velocity_o_faces_vec())

    def velocity_faces_vec():
        doc = "The velocity_faces_vec property."
        def fget(self):
            return self.velocity_w_faces_vec + self.velocity_o_faces_vec
        return locals()
    velocity_faces_vec = property(**velocity_faces_vec())

    def flux_w_faces_vec():
        doc = "The flux_w_faces_vec property."
        def fget(self):
            areas = self.data_impress['area'].reshape(len(self.data_impress['area']), 1)
            return self.velocity_w_faces_vec*areas
        return locals()
    flux_w_faces_vec = property(**flux_w_faces_vec())

    def flux_o_faces_vec():
        doc = "The flux_o_faces_vec property."
        def fget(self):
            areas = self.data_impress['area'].reshape(len(self.data_impress['area']), 1)
            return self.velocity_o_faces_vec*areas
        return locals()
    flux_o_faces_vec = property(**flux_o_faces_vec())

    def flux_faces_vec():
        doc = "The flux_faces_vec property."
        def fget(self):
            return self.flux_w_faces_vec + self.flux_o_faces_vec
        return locals()
    flux_faces_vec = property(**flux_faces_vec())

    def change_upwind_o(self):
        self._data['upwind_identificate_o'] = ~self._data['upwind_identificate_o']

    def change_upwind_w(self):
        self._data['upwind_identificate'] = ~self._data['upwind_identificate']

    def faces_normal_with_z():
        doc = "The faces_normal_with_z property."
        def fget(self):
            # import pdb; pdb.set_trace()
            test = data_loaded['test_segregation']
            if test == True:
                try:
                    return self._faces_normal_with_z
                except AttributeError:
                    faces = self.elements_lv0['faces']
                    u_normal_faces = self.data_impress['u_normal']
                    tt = np.absolute(u_normal_faces[:,2])==1
                    faces_normal_with_z = faces[tt]
                    self._faces_normal_with_z = faces_normal_with_z
                    return self._faces_normal_with_z
            else:
                return np.array([], dtype=int)

        return locals()
    faces_normal_with_z = property(**faces_normal_with_z())

    def flux_sigma_internal_faces():
        doc = "The flux_sigma_internal_faces property."
        def fget(self):
            internal_faces = self.elements_lv0['internal_faces']
            area_int = self.data_impress['area'][internal_faces]
            k_harm_internal_faces = self.data_impress['k_harm'][internal_faces]
            f1 = -(self.grad_p_internal_faces*area_int*k_harm_internal_faces*self.lambda_t_internal_faces - self.flux_grav_total_internal_faces)
            return f1
        return locals()
    flux_sigma_internal_faces = property(**flux_sigma_internal_faces())

    def upwind_w():
        doc = "The upwind_w property."
        def fget(self):
            v0 = self.elements_lv0['neig_internal_faces']
            upwind_w = self._data['upwind_identificate']
            upwind_w = v0[upwind_w]
            return upwind_w
        return locals()
    upwind_w = property(**upwind_w())

    def upwind_o():
        doc = "The upwind_o property."
        def fget(self):
            v0 = self.elements_lv0['neig_internal_faces']
            upwind_o = self._data['upwind_identificate_o']
            upwind_o = v0[upwind_w]
            return upwind_o
        return locals()
    upwind_o = property(**upwind_o())
