from .. import directories as direc
from ..utils.utils_old import get_box, getting_tag
from pymoab import types
import numpy as np
from ..data_class.data_manage import dataManager
import collections

class Contours:

    def __init__(self, M):
        self._gravity = direc.data_loaded['gravity']
        self._loaded = False
        # self.ws_p = [] # pocos de pressao prescrita
        # self.ws_q = []  # pocos de vazao prescrita
        # self.ws_inj = [] # pocos injetores
        # self.ws_pro = [] # pocos produtores
        # self.values_p = [] # valores de pressao prescrita
        # self.values_q = [] # valores de vazao prescrita
        self.name_datas = direc.names_datas_contour
        self.name_file = direc.names_outfiles_steps[4]
        # self.datas = dict()
        self.datas = dataManager('contours.npz')
        self.tags = dict()
        self.tags_to_infos = dict()
        self.names = ['ws_p', 'ws_q', 'ws_inj', 'ws_prod', 'values_p', 'values_q', 'all_wells']
        M.contours = self
        self.mesh = M

    def add_gravity(self, M, gama):
        assert direc.data_loaded['gravity'] == True

        cent_nodes = M.data.variables[M.data.variables_impress['NODES']]
        self.Lz = cent_nodes.max(axis=0)[2]

        ws_p = self.datas['ws_p']
        values_p_ini = self.datas['values_p_ini']

        zs_ws_p = M.data['centroid_volumes'][ws_p][:,2]
        gama_ws_p = gama[ws_p]

        dz = gama_ws_p*(-zs_ws_p + self.Lz)
        values_p = values_p_ini + dz
        self.datas['values_p'] = values_p

    def create_tags(self, M):
        assert not self._loaded

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

    def get_wells(self, M):
        assert not self._loaded

        data_wells = direc.data_loaded['Wells']
        centroids = M.data['centroid_volumes']
        gravity = direc.data_loaded['gravity']

        ws_p = []
        ws_q = []
        ws_inj = []
        ws_prod = []
        values_p = []
        values_q = []

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

                    val = value/nv
                    if tipo == 'Injector':
                        val *= -1

                    ws_q.append(vols)
                    values_q.append(np.repeat(val, nv))

                elif prescription == 'P':
                    val = value
                    ws_p.append(vols)
                    values_p.append(np.repeat(val, nv))

                if tipo == 'Injector':
                    ws_inj.append(vols)
                elif tipo == 'Producer':
                    ws_prod.append(vols)

        # self.ws_q = np.array(self.ws_q).flatten()
        # self.ws_p = np.array(self.ws_p).flatten()
        # self.values_p = np.array(self.values_p).flatten()
        # self.values_q = np.array(self.values_q).flatten()
        # self.ws_inj = np.array(self.ws_inj).flatten()
        # self.ws_pro = np.array(self.ws_pro).flatten()

        ws_q = np.array(ws_q).flatten()
        ws_p = np.array(ws_p).flatten()
        values_p = np.array(values_p).flatten()
        values_q = np.array(values_q).flatten()
        ws_inj = np.array(ws_inj).flatten()
        ws_prod = np.array(ws_prod).flatten()

        self.datas['ws_p'] = ws_p
        self.datas['ws_q'] = ws_q
        self.datas['ws_inj'] = ws_inj
        self.datas['ws_prod'] = ws_prod
        self.datas['values_p'] = values_p
        self.datas['values_q'] = values_q
        self.datas['all_wells'] = np.union1d(ws_inj, ws_prod)
        self.datas['values_p_ini'] = values_p.copy()

    def set_infos(self, M):
        assert not self._loaded

        mb = M.core.mb

        all_volumes = np.array(M.core.all_volumes)
        ws_p = all_volumes[self.datas['ws_p'].flatten()]
        ws_q = all_volumes[self.datas['ws_q'].flatten()]
        ws_prod = all_volumes[self.datas['ws_prod'].flatten()]
        ws_inj = all_volumes[self.datas['ws_inj'].flatten()]
        values_p = self.datas['values_p'].flatten()
        values_q = self.datas['values_q'].flatten()

        mb.tag_set_data(self.tags['INJ'], ws_inj, np.repeat(1, len(ws_inj)))
        mb.tag_set_data(self.tags['PROD'], ws_prod, np.repeat(1, len(ws_prod)))
        mb.tag_set_data(self.tags['P'], ws_p, values_p)
        mb.tag_set_data(self.tags['Q'], ws_q, values_q)

    def load_tags(self, M):
        assert not self._loaded

        mb = M.core.mb

        names = ['P', 'Q', 'INJ', 'PROD']
        for name in names:
            self.tags[name] = mb.tag_get_handle(name)

    def export_to_npz_dep0(self):
        assert not self._loaded

        file_name = self.name_datas

        np.savez(file_name, **self.datas)

    def export_to_npz(self):

        self.datas.export_to_npz()

    def load_from_npz_dep0(self):
        assert not self._loaded
        file_name = self.name_datas

        arq = np.load(file_name)

        for name, values in arq.items():
            self.datas[name] = values

    def load_from_npz(self):

        self.datas.load_from_npz()

    def save_mesh(self, M):
        M.state = 4
        np.save(direc.state_path, np.array([M.state]))
        np.save(direc.path_local_last_file_name, np.array([direc.names_outfiles_steps[4]]))
        M.core.print(file=direc.output_file+str(M.state))

    def update_values(self):
        M = self.mesh

        self.mesh.core.mb.tag_set_data(self.tags['P'], M.core.all_volumes[self.datas['ws_p']], self.datas['values_p'])
        self.mesh.core.mb.tag_set_data(self.tags['Q'], M.core.all_volumes[self.datas['ws_q']], self.datas['values_q'])

    def correct_wells(self):
        M = self.mesh
        wells_q = self.datas['ws_q']

        fc_n = M.volumes.bridge_adjacencies(wells_q, 3, 2).flatten()
        contador = collections.Counter(fc_n)
        facs_nn = np.array([k for k, v in contador.items() if v > 1])
        k_harm_faces = M.data['k_harm'].copy()
        k_max = k_harm_faces.max()
        k_harm_faces[facs_nn] = np.repeat(k_max, len(facs_nn))

        areas = M.data['area']
        dist_cent = M.data['dist_cent']
        pretransmissibility_faces = (areas*k_harm_faces)/dist_cent

        M.data['k_harm'] = k_harm_faces
        M.data[M.data.variables_impress['pretransmissibility']] = pretransmissibility_faces
        M.data.update_variables_to_mesh(['k_harm', 'pretransmissibility'])

    def loaded(self):
        assert not self._loaded
        self._loaded = True
