from packs.data_class import DataManager
from packs.utils.utils_old import get_box
import numpy as np


class WellsData(DataManager):

    def add_gravity_2(self, volumes, gravity_vector, centroid_volumes, volumes_adj_volumes_by_faces, centroid_nodes, saturation, rho_w, rho_o):

        # M = self.mesh
        # gama = M.data['gama']
        ## gama_vector = (rho_w + rho)
        cent_nodes = centroid_nodes
        # self.Lz = cent_nodes.max(axis=0)[2]
        # Lsmax = cent_nodes.max(axis=0)
        Lsmax = centroid_volumes.max(axis=0)
        rho_volumes = rho_w*saturation + rho_o*(1-saturation)
        gama = -gravity_vector[2]*rho_volumes

        ws_p = self._data['ws_p']
        centroids = centroid_volumes

        if len(ws_p) < 1:
            return 0
        values_p_ini = self._data['values_p_ini']

        lim = 1e-10
        delta = centroids.min()/8
        ws_p_sep = self['ws_p_sep']
        values_p_ini_sep = self['values_p_ini_sep']


        for ps, values in zip(ws_p_sep, values_p_ini_sep):
            self['presc_pressure'][ps] = values
            centroids_ps = centroids[ps]
            # x_max_ws_p = centroids_ps[:,0].max()
            # y_max_ws_p = centroids_ps[:,1].max()
            # z_max_ws_p = centroids_ps[:,2].max()
            xmax, ymax, zmax = centroids_ps.max(axis=0)
            xmin, ymin, zmin = centroids_ps.min(axis=0)
            box = np.array([np.array([xmax-delta, ymin-delta, zmin-delta]), np.array([xmax+delta, ymax+delta, zmax+delta])])
            ps_xmax = get_box(centroids_ps, box)
            box = np.array([np.array([xmin-delta, ymax-delta, zmin-delta]), np.array([xmax+delta, ymax+delta, zmax+delta])])
            ps_ymax = get_box(centroids_ps, box)
            box = np.array([np.array([xmin-delta, ymin-delta, zmax-delta]), np.array([xmax+delta, ymax+delta, zmax+delta])])
            ps_zmax = get_box(centroids_ps, box)
            ps_zmax = ps[ps_zmax]
            ps_ymax = ps[ps_ymax]
            ps_xmax = ps[ps_xmax]
            value = values[0]
            vols_2_all = []
            for ppp in ps_zmax:
                # viz_zmax = volumes[self.elements_lv0['adj_matrix_volumes_volumes'][ppp]]
                viz_zmax = volumes_adj_volumes_by_faces[ppp]
                centroids_viz_zmax = centroids[viz_zmax]
                delta_z = centroids[ppp][2] - centroids_viz_zmax[:,2]
                ident_lim = delta_z > 0
                vols_2 = viz_zmax[ident_lim]
                delta_z = delta_z[ident_lim]

                for i, deltz in enumerate(delta_z):
                    if deltz > 0:
                        if vols_2[i] in ps:
                            vols_2_all.append(vols_2[i])
                            self['presc_pressure'][vols_2[i]] = self['presc_pressure'][ppp] + gama[vols_2[i]]*deltz

            while len(vols_2_all) > 0:
                vols_2_all = self.calc_presc_pressure(vols_2_all, ps, volumes, centroids, gama, volumes_adj_volumes_by_faces)

        # zs_ws_p = M.data['centroid_volumes'][ws_p][:,2]
        # gama_ws_p = gama[ws_p]
        #
        # dz = gama_ws_p*(-zs_ws_p + self.Lz)
        # values_p = values_p_ini + dz
        # self._data['values_p'] = values_p

        self['values_p'] = self['presc_pressure'][self['ws_p']]

    def calc_presc_pressure(self, vols_2_all, ps, volumes, centroids, gama, volumes_adj_volumes_by_faces):
        vols_2_all_2 = []
        for ppp in vols_2_all:
            viz_zmax = volumes_adj_volumes_by_faces[ppp]
            centroids_viz_zmax = centroids[viz_zmax]
            delta_z = centroids[ppp][2] - centroids_viz_zmax[:,2]
            ident_lim = delta_z > 0
            vols_2 = viz_zmax[ident_lim]
            delta_z = delta_z[ident_lim]

            for i, deltz in enumerate(delta_z):
                if deltz > 0:
                    if vols_2[i] in ps:
                        vols_2_all_2.append(vols_2[i])
                        self['presc_pressure'][vols_2[i]] = self['presc_pressure'][ppp] + gama[vols_2[i]]*deltz

        return vols_2_all_2
