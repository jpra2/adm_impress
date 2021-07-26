from .. import directories as direc
from ..utils.utils_old import get_box, getting_tag
from pymoab import types
import numpy as np
from ..data_class.data_manager import DataManager
import collections
from ..utils.utils_old import get_box
import pdb

class Wells(DataManager):

    def __init__(self, M, elements_lv0, load: bool=False, data_name: str='wells.npz'):
        super().__init__(data_name, load=load)
        self._gravity = direc.data_loaded['gravity']
        self.tags = dict()
        self.tags_to_infos = dict()
        self.names = ['ws_p', 'ws_q', 'ws_inj', 'ws_prod', 'values_p', 'values_q', 'all_wells']
        self.mesh = M
        self.elements_lv0 = elements_lv0
        if not load:
            self.run()
        else:
            # self.load_tags()
            self.create_tags()
            self.set_infos()
            self._loaded = True

        self['presc_pressure'] = np.zeros(len(self.elements_lv0['volumes']))

    def add_gravity(self):
        assert direc.data_loaded['gravity'] == True
        volumes = self.elements_lv0['volumes']

        M = self.mesh
        gama = M.data['gama']
        cent_nodes = M.data['centroid_nodes']
        self.Lz = cent_nodes.max(axis=0)[2]

        ws_p = self._data['ws_p']
        centroids = M.data['centroid_volumes']

        if len(ws_p) < 1:
            return 0
        values_p_ini = self._data['values_p_ini']

        lim = 1e-10
        delta = centroids.min()/4
        ws_p_sep = self['ws_p_sep']
        values_p_ini_sep = self['values_p_ini_sep']

        for ps, values in zip(ws_p_sep, values_p_ini_sep):
            self['presc_pressure'][ps] = values
            centroids_ps = centroids[ps]
            z_max_ws_p = centroids_ps[:,2].max()
            xmax, ymax, zmax = centroids_ps.max(axis=0)
            xmin, ymin, zmin = centroids_ps.min(axis=0)
            box = np.array([np.array([xmin-delta, ymin-delta, z_max_ws_p-delta]), np.array([xmax+delta, ymax+delta, z_max_ws_p+delta])])
            ps_zmax = get_box(centroids_ps, box)
            ps_zmax = ps[ps_zmax]
            value = values[0]
            vols_2_all = []
            for ppp in ps_zmax:
                viz_zmax = volumes[self.elements_lv0['adj_matrix_volumes_volumes'][ppp]]
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
                vols_2_all = self.calc_presc_pressure(vols_2_all, ps, volumes, centroids, gama)

        # zs_ws_p = M.data['centroid_volumes'][ws_p][:,2]
        # gama_ws_p = gama[ws_p]
        #
        # dz = gama_ws_p*(-zs_ws_p + self.Lz)
        # values_p = values_p_ini + dz
        # self._data['values_p'] = values_p

        self['values_p'] = self['presc_pressure'][self['ws_p']]

    def add_gravity_2(self, volumes, gravity_vector, centroid_volumes, volumes_adj_volumes_by_faces, centroid_nodes, saturation, rho_w, rho_o):

        # M = self.mesh
        # gama = M.data['gama']
        ## gama_vector = (rho_w + rho)
        cent_nodes = centroid_nodes
        # self.Lz = cent_nodes.max(axis=0)[2]
        Lsmax = cent_nodes.max(axis=0)
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
                vols_2_all = self.calc_presc_pressure(vols_2_all, ps, volumes, centroids, gama)

        # zs_ws_p = M.data['centroid_volumes'][ws_p][:,2]
        # gama_ws_p = gama[ws_p]
        #
        # dz = gama_ws_p*(-zs_ws_p + self.Lz)
        # values_p = values_p_ini + dz
        # self._data['values_p'] = values_p

        self['values_p'] = self['presc_pressure'][self['ws_p']]

    def calc_presc_pressure(self, vols_2_all, ps, volumes, centroids, gama):
        vols_2_all_2 = []
        for ppp in vols_2_all:
            viz_zmax = volumes[self.elements_lv0['adj_matrix_volumes_volumes'][ppp]]
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

    def create_tags(self):
        assert not self._loaded
        M = self.mesh
        mb = M.core.mb

        l = ['P', 'Q']
        for name in l:
            n = 1
            tipo = 'double'
            entitie = 'volumes'
            t1 = types.MB_TYPE_DOUBLE
            t2 = types.MB_TAG_SPARSE
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        l = ['INJ', 'PROD']
        for name in l:
            n = 1
            tipo = 'integer'
            entitie = 'volumes'
            t1 = types.MB_TYPE_INTEGER
            t2 = types.MB_TAG_SPARSE
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

    def get_wells(self):
        assert not self._loaded
        M = self.mesh

        data_wells = direc.data_loaded['Wells']
        centroids = M.data['centroid_volumes']
        gravity = direc.data_loaded['gravity']

        ws_p = [] ## pocos com pressao prescrita
        ws_q = [] ## pocos com vazao prescrita
        ws_inj = [] ## pocos injetores
        ws_prod = [] ## pocos produtores
        values_p = [] ## valor da pressao prescrita
        values_q = [] ## valor da vazao prescrita

        for p in data_wells:

            well = data_wells[p]
            type_region = well['type_region']
            tipo = well['type']
            prescription = well['prescription']
            value = well['value']

            if type_region == direc.types_region_data_loaded[1]: #box

                p0 = well['p0']
                p1 = well['p1']
                limites = np.array([p0, p1])
                vols = get_box(centroids, limites)
                nv = len(vols)

                if prescription == 'Q':
                    try:
                        val = value/nv
                    except:
                        print("Nenhum volumes corresponde ao poço")
                    
                    if tipo == 'Injector':
                        val *= -1

                    ws_q.append(vols)
                    values_q.append(np.repeat(val, nv))

                elif prescription == 'P':
                    val = value
                    ws_p.append(vols.astype(np.int32))
                    values_p.append(np.repeat(val, nv))

                if tipo == 'Injector':
                    ws_inj.append(vols)
                elif tipo == 'Producer':
                    ws_prod.append(vols)

        self['ws_q_sep'] = ws_q
        self['ws_p_sep'] = ws_p
        self['values_p_sep'] = np.array(values_p)
        self['values_p_ini_sep'] = self['values_p_sep'].copy()
        self['values_q_sep'] = np.array(values_q)

        # import pdb; pdb.set_trace()

        if len(ws_q) == 0:
            ws_q = np.array(ws_q, dtype=int)
        else:
            ws_q = np.concatenate(ws_q)

        if len(ws_p) == 0:
            ws_p = np.array(ws_p, dtype=int)
        else:
            ws_p = np.concatenate(ws_p)


        # ws_q = np.concatenate(ws_q)
        # ws_p = np.concatenate(ws_p)
        # values_p = np.concatenate(values_p)

        if len(values_p) == 0:
            values_p = np.array(values_p)
        else:
            values_p = np.concatenate(values_p)

        if len(values_q) == 0:
            values_q = np.array(values_q)
        else:
            values_q = np.concatenate(values_q)

        # values_q = np.concatenate(values_q)

        if len(ws_inj) == 0:
            ws_inj = np.array(ws_inj, dtype=int)
        else:
            ws_inj = np.concatenate(ws_inj)

        if len(ws_prod) == 0:
            ws_prod = np.array(ws_prod, dtype=int)
        else:
            ws_prod = np.concatenate(ws_prod)

        # ws_inj = np.concatenate(ws_inj)
        # ws_prod = np.concatenate(ws_prod)

        self['ws_p'] = ws_p.astype(np.int64)
        self['ws_q'] = ws_q.astype(np.int64)
        self['ws_inj'] = ws_inj.astype(np.int64)
        self['ws_prod'] = ws_prod.astype(np.int64)
        self['values_p'] = values_p
        self['values_q'] = values_q
        self['all_wells'] = np.union1d(ws_inj, ws_prod).astype(int)
        self['values_p_ini'] = values_p.copy()

    def set_infos(self):
        assert not self._loaded
        M = self.mesh

        mb = M.core.mb

        all_volumes = np.array(M.core.all_volumes)
        ws_p = all_volumes[self['ws_p']]
        ws_q = all_volumes[self['ws_q']]
        ws_prod = all_volumes[self['ws_prod']]
        ws_inj = all_volumes[self['ws_inj']]
        values_p = self['values_p']
        values_q = self['values_q']

        mb.tag_set_data(self.tags['INJ'], ws_inj, np.repeat(1, len(ws_inj)))
        mb.tag_set_data(self.tags['PROD'], ws_prod, np.repeat(1, len(ws_prod)))
        mb.tag_set_data(self.tags['P'], ws_p, values_p)

        if len(ws_q) > 0:
            mb.tag_set_data(self.tags['Q'], ws_q, values_q.sum(axis=0))

    def load_tags(self):
        assert not self._loaded
        M = self.mesh

        mb = M.core.mb

        names = ['P', 'Q', 'INJ', 'PROD']
        for name in names:
            self.tags[name] = mb.tag_get_handle(name)

    def save_mesh(self):
        M = self.mesh
        M.state = 4
        np.save(direc.state_path, np.array([M.state]))
        np.save(direc.path_local_last_file_name, np.array([direc.names_outfiles_steps[4]]))
        M.core.print(file=direc.output_file+str(M.state))

    def update_values_to_mesh(self):
        M = self.mesh

        self.mesh.core.mb.tag_set_data(self.tags['P'], M.core.all_volumes[self['ws_p']], self['values_p'])
        self.mesh.core.mb.tag_set_data(self.tags['Q'], M.core.all_volumes[self['ws_q']], self['values_q'])

    def correct_wells(self):
        if len(self['ws_q']) == 0:
            return 0

        M = self.mesh
        wells_q = self['ws_q']

        facs_nn = self['facs_nn']
        k_harm_faces = M.data['k_harm'].copy()
        k_max = k_harm_faces.max()

        k_harm_faces[facs_nn] = np.repeat(k_max, len(facs_nn))

        areas = M.data['area']
        dist_cent = M.data['dist_cent']
        pretransmissibility_faces = (areas*k_harm_faces)/dist_cent

        M.data['k_harm'] = k_harm_faces
        M.data[M.data.variables_impress['pretransmissibility']] = pretransmissibility_faces

    def get_facs_nn(self):

        assert not self._loaded
        M = self.mesh
        wells_q = self['ws_q']

        if len(wells_q) > 0:
            fc_n = M.volumes.bridge_adjacencies(wells_q, 3, 2).flatten()
            contador = collections.Counter(fc_n)
            facs_nn = np.array([k for k, v in contador.items() if v > 1], dtype=np.int64)
            self['facs_nn'] = facs_nn
        else:
            self['facs_nn'] = np.array([], dtype=np.int64)

    def loaded(self):
        assert not self._loaded
        self._loaded = True

    def run(self):
        self.create_tags()
        self.get_wells()
        self.set_infos()
        # self.get_facs_nn()
        # self.correct_wells()
        self.loaded()
        pass
