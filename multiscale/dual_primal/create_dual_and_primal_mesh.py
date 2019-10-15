import directories as direc
from utils.utils_old import getting_tag
# from pymoab import core, types, rng, topo_util
from pymoab import types

class DualPrimalMesh1:

    def __init__(self, M):
        self._loaded = False
        self.tags = dict()
        self._internals = dict()
        self._faces = dict()
        self._edges = dict()
        self._vertex = dict()

    def create_tags(self, M):
        if self._loaded:
            return 1

        mb = M.core.mb

        l = ['D1', 'D2', 'FINE_TO_PRIMAL_CLASSIC_1', 'FINE_TO_PRIMAL_CLASSIC_2']
        for name in l:
            n = 1
            tipo = 'integer'
            entitie = 'volumes'
            t1 = types.MB_TYPE_INTEGER
            t2 = types.MB_TAG_SPARSE
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        l = ['PRIMAL_ID_1']
        for name in l:
            n = 1
            tipo = 'integer'
            entitie = 'coarse_volumes_lv1'
            t1 = types.MB_TYPE_INTEGER
            t2 = types.MB_TAG_SPARSE
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        l = ['PRIMAL_ID_2']
        for name in l:
            n = 1
            tipo = 'integer'
            entitie = 'coarse_volumes_lv2'
            t1 = types.MB_TYPE_INTEGER
            t2 = types.MB_TAG_SPARSE
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        l = ['L2_MESHSET']
        for name in l:
            n = 1
            tipo = 'meshset'
            entitie = 'root_set'
            t1 = types.MB_TYPE_HANDLE
            t2 = types.MB_TAG_MESH
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        return 0

    def loaded(self):
        self._loaded = True

    def run(self, M):

        if self._loaded:
            return 1

        M.dualprimal = self
        self.create_tags(M)
        self.set_primal_level_l_meshsets(M)
        self.loaded()

    def set_primal_level_l_meshsets(self, M):
        if self._loaded:
            return 1

        mb = M.core.mb
        l2_meshset = mb.create_meshset()
        mb.tag_set_data(self.tags['L2_MESHSET'], 0, l2_meshset)

        coarse_volumes = M.coarse.elements

        import pdb; pdb.set_trace()


        return 0


    def internals():
        doc = "The internals property."
        def fget(self, n):
            return self._internals[n]
        def fset(self, n, value):
            assert  not self._loaded, 'nao pode alterar internals'
            self._internals[n] = value
        def fdel(self):
            return 0
            # del self._internals
        return locals()
    internals = property(**internals())

    def faces():
        doc = "The internals property."
        def fget(self, n):
            return self._faces[n]
        def fset(self, n, value):
            assert  not self._loaded, 'nao pode alterar internals'
            self._faces[n] = value
        def fdel(self):
            return 0
            # del self._internals
        return locals()
    faces = property(**internals())
