import directories as direc
from utils.utils_old import get_box, getting_tag
from pymoab import types

class Contours:

    def __init__(self, M):
        self.ws_p = []
        self.ws_q = []
        self.ws_inj = []
        self.ws_pro = []
        self.values_p = []
        self.values_q = []
        self.tags = dict()
        self.tags_to_infos = dict()
        M.contours = self

    def create_tags(self, M):

        mb = M.core.mb

        l = ['P', 'Q']
        for name in l:
            n = 1
            tipo = 'double'
            entitie = 'volumes'
            t1 = types.MB_TYPE_DOUBLE
            t2 = types.MB_TAG_SPARSE
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        l = ['P', 'Q']
        for name in l:
            n = 1
            tipo = 'integer'
            entitie = 'volumes'
            t1 = types.MB_TYPE_INTEGER
            t2 = types.MB_TAG_SPARSE
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        pass

    def load_wells(self, M):

        data_wells = direc.data_loaded['Wells']
        centroids = M.data.centroids[direc.entities_lv0[3]]
        volumes = M.data.elements_lv0[direc.entities_lv0[3]]

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
                vols = get_box(cenroids, limites)
                nv = len(vols)
                if prescription == 'Q':

                    val = value/n
                    if tipo == 'Producer':
                        val *= -/1

                    self.ws_q.append(vols)
                    self.values_q.append(np.repeat(val, nv))

                if prescription == 'P':
                    val = value
                    self.ws_p.append(vols)
                    self.values_p.append(np.repeat(val, nv))

                if tipo == 'Injector':
                    self.ws_inj.append(vols)
                elif tipo == 'Producer':
                    self.ws_pro.append(vols)
