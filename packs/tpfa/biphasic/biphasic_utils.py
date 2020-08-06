
class BiphasicUtils:

    def mobility_w_volumes(self, krw, mi_w):
        return krw/mi_w

    def mobility_o_volumes(self, kro, mi_o):
        return kro/mi_o

    def total_mobility_volumes(self, kro, krw, mi_w, mi_o):
        return self.mobility_w(krw, mi_w) + self.mobility_o(kro, mi_o)

    def fw_volumes(self, krw, kro, mi_w, mi_o):
        return self.mobility_w_volumes(krw, mi_w)/(self.total_mobility_volumes(kro, krw, mi_w, mi_o))

    def mobility_rho_w(self, krw, mi_w, rho_w):
        return self.mobility_w_volumes(krw, mi_w)*rho_w

    def mobility_rho_o(self, kro, mi_o, rho_o):
        return self.mobility_o_volumes(kro, mi_o)*rho_o

    def mobility_rho_sum(self, krw, kro, mi_w, mi_o, rho_w, rho_o):
        return self.mobility_rho_w(krw, mi_w, rho_w) + self.mobility_rho_o(kro, mi_o, rho_o)

    def g_w(self, krw, mi_w, rho_w, gravity_vector):
        return self.mobility_rho_w(krw, mi_w, rho_w)*gravity_vector

    def g_o(self, kro, mi_o, rho_o, gravity_vector):
        return self.mobility_rho_o(kro, mi_o, rho_o)*gravity_vector

    def g_sum(self, krw, kro, mi_w, mi_o, rho_w, rho_o, gravity_vector):
        return self.g_w(krw, mi_w, rho_w, gravity_vector) + self.g_o(kro, mi_o, rho_o, gravity_vector)

    def T_volumes(self, permeability, hi_volume):
        return permeability/hi_volume

    def M_volumes_w(self, permeability, hi_volume, krw, mi_w):
        return self.T_volumes(permeability, hi_volume)*self.mobility_w_volumes(krw, mi_w)

    def M_volumes_o(self, permeability, hi_volume, kro, mi_o):
        return self.T_volumes(permeability, hi_volume)*self.mobility_o_volumes(kro, mi_o)

    def M_volumes_sum(self, permeability, hi_volume, krw, kro, mi_w, mi_o):
        return self.M_volumes_w(permeability, hi_volume, krw, mi_w) + self.M_volumes_o(permeability, hi_volume, kro, mi_o)

    def gravity_source_internal_faces_w_cons(self, krw, volumes, permeability, u_normal_internal_faces, gravity_vector, areas_internal_faces, volumes_adj_internal_faces, hi):

        abs_u_normal_faces = np.absolute(u_normal_internal_faces)
        abs_u_normal_faces = abs_u_normal_faces.reshape(len(abs_u_normal_faces), 3, 1)

        pdb.set_trace()

        perm_adj_volumes = permeability[volumes_adj_internal_faces]

        # k0 = np.tensordot(perm_adj_volumes,



        # krw_adj_volumes = krw[volumes_adj_internal_faces]
        # permeability_adj_volumes = permeability[volumes_adj_internal_faces]
        #
        # for i in range(len(hi)):
        #     perms_adjs = permeability_adj_volumes[i]
        #     u_normal = u_normal_internal_faces[i]
        #     area = areas_internal_faces[i]
        #     Mr0 = self.

        pass
