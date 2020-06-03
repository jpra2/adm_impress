from .. import directories as direc
from ..utils.utils_old import get_box, getting_tag
from pymoab import types
import numpy as np
from ..data_class.data_manager import DataManager
import collections

class Wells(DataManager):

    def __init__(self, M, load: bool=False, data_name: str='wells.npz'):
        super().__init__(data_name, load=load)
        self._gravity = direc.data_loaded['gravity']
        self.tags = dict()
        self.tags_to_infos = dict()
        self.names = ['ws_p', 'ws_q', 'ws_inj', 'ws_prod', 'values_p', 'values_q', 'all_wells']
        self.mesh = M
        if not load:
            self.run()
        else:
            # self.load_tags()
            self.create_tags()
            self.set_infos()
            self._loaded = True

    def add_gravity(self):
        assert direc.data_loaded['gravity'] == True
        M = self.mesh
        gama = M.data['gama']
        cent_nodes = M.data['centroid_nodes']
        self.Lz = cent_nodes.max(axis=0)[2]

        ws_p = self._data['ws_p']
        if len(ws_p) < 1:
            return 0
        values_p_ini = self._data['values_p_ini']

        zs_ws_p = M.data['centroid_volumes'][ws_p][:,2]
        gama_ws_p = gama[ws_p]

        dz = gama_ws_p*(-zs_ws_p + self.Lz)
        values_p = values_p_ini + dz
        self._data['values_p'] = values_p

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
        values_q_type = []

        for p in data_wells:

            well = data_wells[p]
            type_region = well['type_region']
            tipo = well['type']
            prescription = well['prescription']
            value = np.array(well['value']).astype(float)

            if type_region == direc.types_region_data_loaded[1]: #box

                p0 = well['p0']
                p1 = well['p1']
                limites = np.array([p0, p1])
                vols = get_box(centroids, limites)
                if len(values_q)==0 and prescription=='Q':
                    values_q = np.zeros([len(vols), len(value)])
                nv = len(vols)
                i = 0
                if prescription == 'Q':
                    val = value/nv
                    if tipo == 'Injector':
                        val *= -1
                    ws_q.append(vols)
                    values_type = np.repeat(well['value_type'], nv)
                    values_q[i:nv,:] = val
                    values_q_type.append(values_type)
                    i = nv

                elif prescription == 'P':
                    val = value
                    ws_p.append(vols)
                    values_p.append(np.repeat(val, nv))

                if tipo == 'Injector':
                    ws_inj.append(vols)
                elif tipo == 'Producer':
                    ws_prod.append(vols)

        ws_q = np.array(ws_q).flatten()
        ws_p = np.array(ws_p).flatten()
        values_p = np.array(values_p).flatten()
        values_q = np.array(values_q)#.flatten()
        ws_inj = np.array(ws_inj).flatten()
        ws_prod = np.array(ws_prod).flatten()

        self['ws_p'] = ws_p.astype(int)
        self['ws_q'] = ws_q.astype(int)
        self['ws_inj'] = ws_inj.astype(int)
        self['ws_prod'] = ws_prod.astype(int)
        self['values_p'] = values_p
        self['values_q'] = values_q
        self['all_wells'] = np.union1d(ws_inj, ws_prod)
        self['values_p_ini'] = values_p.copy()
        self['value_type'] = values_q_type

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
        if len(values_q>0):
            for i in range(len(values_q[0,:])):
                mb.tag_set_data(self.tags['Q'], ws_q, values_q[:,i])


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

        # fc_n = M.volumes.bridge_adjacencies(wells_q, 3, 2).flatten()
        # contador = collections.Counter(fc_n)
        # facs_nn = np.array([k for k, v in contador.items() if v > 1])
        # self['facs_nn'] = facs_nn

        self['facs_nn'] = []

    def loaded(self):
        assert not self._loaded
        self._loaded = True

    def run(self):
        self.create_tags()
        self.get_wells()
        self.set_infos()
        self.get_facs_nn()
        self.correct_wells()
        self.loaded()
        pass
