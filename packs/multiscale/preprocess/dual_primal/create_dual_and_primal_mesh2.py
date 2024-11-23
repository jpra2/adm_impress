from .... import directories as direc
from ....utils.utils_old import getting_tag, get_box, Min_Max, add_topology
# from pymoab import core, types, rng, topo_util
from pymoab import types, rng
import numpy as np
import pdb
from ....data_class.data_manager import DataManager
import scipy.sparse as sp
from ....directories import data_loaded
from ....errors.err import DualStructureError
from ....adm.adm_method import get_levelantids_levelids

from packs.multiscale.preprocess.dual_primal.paralell import paralell_dual_and_primal

import time


class MultilevelData(DataManager):

    def __init__(self, data_impress, M, data_name: str='MultilevelData.npz', load=False):
        carregar = load
        super().__init__(data_name=data_name, load=load)
        self.tags = dict()
        self.tags_to_infos = dict()
        self._carregar = carregar
        M.multilevel_data = self
        self.mesh = M

        self.levels = data_loaded['n_levels']
        self.l1 = 1
        self.data_impress = data_impress

        self.fine_primal_id = 'fine_primal_id_level_'
        self.coarse_volumes = 'coarse_volumes_level_'
        self.coarse_primal_id = 'coarse_primal_id_level_'
        self.fine_dual_id = 'fine_dual_id_level_'
        self.interns = 'interns_level_'
        self.faces = 'faces_level_'
        self.edges = 'edges_level_'
        self.vertex = 'vertex_level_'
        self.meshset_vertices = 'meshset_vertices_level_'
        self.reordered_id = 'reordered_id_'
        self.faces_boundary_meshset_level = 'FACES_BOUNDARY_MESHSETS_LEVEL_'
        self.name_mesh = 'flying/multilevel_data'
        self.coarse_neig_face = 'coarse_neig_face_level_'
        self.coarse_id_neig_face = 'coarse_id_neig_face_level_'
        self.restriction = 'restriction_level_'
        self.coarse_faces = 'coarse_faces_level_'
        self.coarse_internal_faces = 'coarse_internal_faces_level_'
        self.coarse_intersect_faces = 'coarse_intersect_faces_level_'
        self.fine_vertex_coarse_volumes = 'fine_vertex_coarse_volumes_level_'
        self.neig_intersect_faces = 'neig_intersect_faces_level_'
        self.internal_boundary_fine_volumes = 'internal_boundary_fine_volumes_level_'
        self.dual_structure = 'dual_structure_level_'
        self.centroids_name = 'centroids_level_'
        self.volumes_without_grav = 'volumes_without_grav_level_'

    def run(self):
        M = self.mesh

        assert not self._loaded
        import time

        t0=time.time()
        self.create_tags()
        print("Time to create tags: {} seconds".format(time.time()-t0))
        t0=time.time()
        self.generate_dual_and_primal_any_D(M)
        print("Time to create dual: {} seconds".format(time.time()-t0))
        t0=time.time()
        self.get_elements(M)
        print("Time to get elements: {} seconds".format(time.time()-t0))
        t0=time.time()
        self.get_boundary_coarse_faces(M)
        print("Time to get boundary: {} seconds".format(time.time()-t0))
        t0=time.time()
        # self.get_dual_structure()
        # print("Time to dual structure: {} seconds".format(time.time()-t0))
        # t0=time.time()
        self.get_dual_structure_with_graph()
        print("Time to dual structure with graph: {} seconds".format(time.time()-t0))
        t0=time.time()
        self.set_volumes_without_gravity_source_term()
        print("Time to set volumes gravity: {} seconds".format(time.time()-t0))
        t0=time.time()
        self.export_to_npz()
        print("Time to export npz: {} seconds".format(time.time()-t0))
        t0=time.time()
        self.loaded()
        print("Time to loaded: {} seconds".format(time.time()-t0))

    def create_tags(self):
        assert not self._loaded
        M = self.mesh

        mb = M.core.mb

        l = ['D1', 'D2', 'FINE_TO_PRIMAL_CLASSIC_1', 'FINE_TO_PRIMAL_CLASSIC_2', 'reordered_id_1', 'reordered_id_2',
        'local_id_internos', 'local_fac_internos']
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
            entitie = 'coarse_volumes'
            t1 = types.MB_TYPE_INTEGER
            t2 = types.MB_TAG_SPARSE
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        l = ['NEIG_FACE']
        for name in l:
            n = 1
            tipo = 'handle'
            entitie = 'coarse_volumes_lv1'
            t1 = types.MB_TYPE_HANDLE
            t2 = types.MB_TAG_SPARSE
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        l = ['L2_MESHSET', 'MV_1', 'MV_2']
        for name in l:
            n = 1
            tipo = 'handle'
            entitie = 'root_set'
            t1 = types.MB_TYPE_HANDLE
            t2 = types.MB_TAG_MESH
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        return 0

    def load_tags(self):
        assert not self._loaded
        M = self.mesh

        tags0 = ['D', 'FINE_TO_PRIMAL_CLASSIC_', 'PRIMAL_ID_', 'MV_', 'reordered_id_']
        tags1 = ['L2_MESHSET', 'local_id_internos', 'local_fac_internos', 'NEIG_FACE']
        name_tag_faces_boundary_meshsets = self.faces_boundary_meshset_level
        n_levels = 2

        mb = M.core.mb

        for name in tags0:
            for i in range(2):
                j = i + 1
                name2 = name + str(j)
                tag = mb.tag_get_handle(name2)
                self.tags[name2] = tag

        for name in tags1:
            self.tags[name] = mb.tag_get_handle(name)

        for i in range(n_levels):
            name_tag = name_tag_faces_boundary_meshsets + str(i+1)
            tag_boundary = mb.tag_get_handle(name_tag)
            self.tags[name_tag] = tag_boundary

    def loaded(self):
        assert not self._loaded
        self._loaded = True

    def generate_dual_and_primal_any_D(self, M):
        M1=M.core
        M1.all_centroids=M.data["centroid_volumes"]
        M1.primal_id_tag1=self.tags["PRIMAL_ID_1"]
        M1.primal_id_tag2=self.tags["PRIMAL_ID_2"]
        M1.fine_to_primal1_classic_tag=self.tags["FINE_TO_PRIMAL_CLASSIC_1"]
        M1.fine_to_primal2_classic_tag=self.tags["FINE_TO_PRIMAL_CLASSIC_2"]
        M1.D1_tag=self.tags["D1"]
        M1.D2_tag=self.tags["D2"]

        coord_nodes = M.data['centroid_nodes']
        cent_volumes = M.data['centroid_volumes']

        t0=time.time()
        print("creating dual mesh")
        paralell_dual_and_primal.DualPrimal(M1, coord_nodes, cent_volumes, external_vertex_on_boundary=True)
        print(time.time()-t0,"tempo para criar a dual")

    def get_elements(self, M):
        assert not self._loaded

        mb = M.core.mb
        mtu = M.core.mtu
        tags_fine = ['D', 'FINE_TO_PRIMAL_CLASSIC_']
        tags_coarse = ['PRIMAL_ID_']
        coarse_id_impress = 'GID_'
        tag_mv = ['MV_']
        tag_reordered_id = ['reordered_id_']
        all_volumes = M.core.all_volumes
        dict_volumes = dict(zip(all_volumes, M.volumes.all))
        mvs = [M.core.root_set]
        fine_centroids = self.data_impress['centroid_volumes']

        self._data[self.centroids_name + str(0)] = fine_centroids.copy()

        for i in range(2):
            n = i + 1
            level = n
            name_tag_c = tags_coarse[0] + str(n)
            dual_fine_name = tags_fine[0] + str(n)
            primal_fine_name = tags_fine[1] + str(n)
            tag_reord_id = tag_reordered_id[0] + str(n)
            mv = mvs[i]
            n_reord = 0

            interns = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]),
                                                      np.array([0]))
            if n not in [1]:
                mb.tag_set_data(self.tags[tag_reord_id], interns, np.arange(n_reord, len(interns)))
            n_reord += len(interns)
            if n == 1:
                interns = np.array([dict_volumes[k] for k in interns])
            else:
                interns = mb.tag_get_data(self.tags[tags_fine[1] + str(n-1)], interns, flat=True)

            faces = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]),
                                                    np.array([1]))
            if n not in [1]:
                mb.tag_set_data(self.tags[tag_reord_id], faces, np.arange(n_reord, n_reord + len(faces)))
            n_reord += len(faces)
            if n == 1:
                faces = np.array([dict_volumes[k] for k in faces])
            else:
                faces = mb.tag_get_data(self.tags[tags_fine[1] + str(n-1)], faces, flat=True)

            edges = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]),
                                                    np.array([2]))
            if n not in [1]:
                mb.tag_set_data(self.tags[tag_reord_id], edges, np.arange(n_reord, n_reord + len(edges)))
            n_reord += len(edges)
            if n == 1:
                edges = np.array([dict_volumes[k] for k in edges])
            else:
                edges = mb.tag_get_data(self.tags[tags_fine[1] + str(n-1)], edges, flat=True)

            vertex = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]),
                                                     np.array([3]))
            if n not in [1]:
                mb.tag_set_data(self.tags[tag_reord_id], vertex, np.arange(n_reord, n_reord + len(vertex)))
            n_reord += len(vertex)
            ids_fine_vertexes = np.array([dict_volumes[k] for k in vertex])
            if n == 1:
                # vertexes = np.array([dict_volumes[k] for k in vertex])
                vertexes = ids_fine_vertexes
            else:
                vertexes = mb.tag_get_data(self.tags[tags_fine[1] + str(n-1)], vertex, flat=True)

            self._data[self.fine_vertex_coarse_volumes+str(level)] = ids_fine_vertexes
            self._data[self.interns + str(level)] = interns
            self._data[self.faces + str(level)] = faces
            self._data[self.edges + str(level)] = edges
            self._data[self.vertex + str(level)] = vertexes
            wire_num = np.array([len(interns), len(faces), len(edges), len(vertexes)])
            nv = wire_num[3]

            coarse_volumes = []
            coarse_primal_ids = []
            coarse_neig_face = []
            coarse_id_neig_face = []
            lines_r = []
            cols_r = []

            contador_coarse_gids = 0
            centroids_coarse = np.zeros((len(vertex), 3), dtype=float)

            for i, vert in enumerate(vertex):
                neigs = []
                neigs_ids = []
                primal_id = mb.tag_get_data(self.tags[primal_fine_name], vert, flat=True)[0]
                centroid_vert = fine_centroids[dict_volumes[vert]]
                centroids_coarse[primal_id] = centroid_vert

                coarse_volume = \
                mb.get_entities_by_type_and_tag(0, types.MBENTITYSET, np.array([self.tags[name_tag_c]]),
                                                np.array([primal_id]))[0]
                coarse_volumes.append(coarse_volume)
                coarse_primal_ids.append(primal_id)
                elems_in_meshset = mb.get_entities_by_handle(coarse_volume)
                n_elems = len(elems_in_meshset)
                gggids = np.array([dict_volumes[k] for k in elems_in_meshset])
                local_id = np.arange(n_elems)
                coarse_global_id = np.arange(contador_coarse_gids, contador_coarse_gids + n_elems)
                self.data_impress['COARSE_GID_'+str(level)][gggids] = coarse_global_id
                self.data_impress['COARSE_LOCAL_ID_'+str(level)][gggids] = local_id
                contador_coarse_gids += n_elems
                if n == 1:
                    gids = gggids
                    # gids = mb.tag_get_data(self.tags[tag_reord_id], elems_in_meshset, flat=True)
                else:
                    gids = np.unique(mb.tag_get_data(self.tags[tags_fine[1] + str(n-1)], elems_in_meshset, flat=True))
                elems_fora = mtu.get_bridge_adjacencies(elems_in_meshset, 2, 3)
                elems_fora = rng.subtract(elems_fora, elems_in_meshset)
                ids_meshsets_vizinhos = np.unique(mb.tag_get_data(self.tags[primal_fine_name], elems_fora, flat=True))                
                for j in ids_meshsets_vizinhos:
                    m2 = mb.get_entities_by_type_and_tag(M.core.root_set, types.MBENTITYSET, np.array([self.tags[name_tag_c]]), np.array([j]))[0]
                    neigs.append(m2)
                    neigs_ids.append(j)

                neigs = np.array(neigs)
                neigs_ids = np.array(neigs_ids)
                coarse_neig_face.append(neigs)
                coarse_id_neig_face.append(neigs_ids)
                if level == 1:
                    d2 = mb.tag_get_data(self.tags[tags_fine[0] + str(level+1)], vert, flat=True)[0]
                    mb.tag_set_data(self.tags[tags_fine[0] + str(level+1)], elems_in_meshset, np.repeat(d2, len(elems_in_meshset)))

            coarse_neig_face = np.array(coarse_neig_face)
            coarse_id_neig_face = np.array(coarse_id_neig_face)
            # centroids_coarse = np.array(centroids_coarse)

            self._data[self.coarse_neig_face + str(level)] = coarse_neig_face
            self._data[self.coarse_id_neig_face + str(level)] = coarse_id_neig_face
            self._data[self.coarse_volumes + str(level)] = np.array(coarse_volumes)
            self._data[self.coarse_primal_id + str(level)] = np.array(coarse_primal_ids)
            self._data[self.centroids_name + str(level)] = centroids_coarse
            # dtype = [('elements', np.uint64), ('id', np.uint64)]
            # structured_array = np.zeros(len(coarse_volumes), dtype=dtype)
            # structured_array['elements'] = np.array(coarse_volumes)
            # structured_array['id'] = np.array(coarse_primal_ids)

            for volume, vizinhos in zip(coarse_volumes, coarse_neig_face):
                m = mb.create_meshset()
                mb.add_entities(m, vizinhos)
                mb.tag_set_data(self.tags['NEIG_FACE'], volume, m)

            nnn = tag_mv[0] + str(level)
            # if not self._carregar:
            mv1 = mb.create_meshset()
            mb.add_entities(mv1, vertex)
            mb.tag_set_data(self.tags[nnn], M.core.root_set, mv1)
            # else:
            # mv1 = mb.tag_get_data(self.tags[nnn], M.core.root_set, flat=True)[0]

            # self.mvs[level] = mv1
            mvs.append(mv1)
            self._data[self.meshset_vertices + str(level)] = mv1

            fine_primal_id = mb.tag_get_data(self.tags[primal_fine_name], all_volumes, flat=True)
            self._data[self.fine_primal_id + str(level)] = fine_primal_id

            fine_dual_id = mb.tag_get_data(self.tags[dual_fine_name], all_volumes, flat=True)
            self._data[self.fine_dual_id + str(level)] = fine_dual_id

        for m in self._data[self.coarse_volumes + str(1)]:
            elements = mb.get_entities_by_handle(m)
            ne = len(elements)
            ids = np.arange(ne)
            dual_info = mb.tag_get_data(self.tags[tags_fine[0] + str(1)], elements, flat=True)
            id_vert = ids[dual_info == 3]
            vertex = elements[id_vert]
            reord_id_2 = mb.tag_get_data(self.tags[tag_reordered_id[0] + str(2)], vertex, flat=True)[0]
            mb.tag_set_data(self.tags[tag_reordered_id[0] + str(2)], elements, np.repeat(reord_id_2, ne))

        self.data_impress['DUAL_1'] = mb.tag_get_data(self.tags['D1'], all_volumes, flat=True)
        self.data_impress['DUAL_2'] = mb.tag_get_data(self.tags['D2'], all_volumes, flat=True)
        self.data_impress[coarse_id_impress + str(2)] = mb.tag_get_data(self.tags['FINE_TO_PRIMAL_CLASSIC_2'], all_volumes, flat=True)
        self.data_impress[coarse_id_impress + str(1)] = mb.tag_get_data(self.tags['FINE_TO_PRIMAL_CLASSIC_1'], all_volumes, flat=True)
        self.data_impress[coarse_id_impress + str(0)] = M.volumes.all

    def get_boundary_coarse_faces(self, M):
        assert not self._loaded
        # meshsets_nv1 = self._coarse_volumes[1]
        # meshsets_nv2 = self._coarse_volumes[2]
        meshsets_nv1 = self[self.coarse_volumes + str(1)]
        meshsets_nv2 = self[self.coarse_volumes + str(2)]

        mb = M.core.mb
        mtu = M.core.mtu
        n_levels = 2
        all_volumes = M.core.all_volumes
        dict_volumes = dict(zip(all_volumes, M.volumes.all))

        name_tag_faces_boundary_meshsets = self.faces_boundary_meshset_level
        all_meshsets = [meshsets_nv1, meshsets_nv2]
        d_faces = dict(zip(M.core.all_faces, M.faces.all))
        b_faces_all = M.faces.boundary

        from ....utils import pymoab_utils as utpy

        for i in range(n_levels):
            level = i+1
            name = name_tag_faces_boundary_meshsets + str(i + 1)
            meshsets = all_meshsets[i]
            n = 1
            tipo = 'handle'
            entitie = 'root_set'
            t1 = types.MB_TYPE_HANDLE
            t2 = types.MB_TAG_MESH
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)
            tag_boundary = self.tags[name]
            utpy.set_faces_in_boundary_by_meshsets(mb, mtu, meshsets, tag_boundary, M)
            faces_boundary = mb.tag_get_data(tag_boundary, M.core.root_set, flat=True)[0]
            faces_boundary = mb.get_entities_by_handle(faces_boundary)
            faces_boundary = np.array([d_faces[k] for k in faces_boundary])
            self._data[self.faces_boundary_meshset_level + str(i+1)] = faces_boundary
            if len(faces_boundary)>0:
                self._data[self.neig_intersect_faces+str(level)] = M.faces.bridge_adjacencies(faces_boundary, 2, 3)


            coarse_faces = []
            coarse_internal_faces = []
            coarse_intersect_faces = []
            coarse_internal_boundary_volumes = []

            for m in meshsets:
                # primal_id = mb.tag_get_data(tag_coarse_id, m, flat=True)[0]
                # assert primal_id == cids[cont]
                elements = mb.get_entities_by_handle(m)
                volumes = np.array([dict_volumes[k] for k in elements])
                faces = mtu.get_bridge_adjacencies(elements, 3, 2)
                faces = np.array([d_faces[k] for k in faces])
                coarse_faces.append(faces)
                internal_faces = np.setdiff1d(faces, faces_boundary)
                internal_faces = np.setdiff1d(internal_faces, b_faces_all)
                coarse_internal_faces.append(internal_faces)
                intersect_faces = np.intersect1d(faces, faces_boundary)
                coarse_intersect_faces.append(intersect_faces)
                boundary_faces = np.setdiff1d(faces, internal_faces)
                # internal_boundary_volumes = np.concatenate(M.faces.bridge_adjacencies(boundary_faces, 2, 3))
                internal_boundary_volumes = np.concatenate(M.faces.bridge_adjacencies(intersect_faces, 2, 3))
                internal_boundary_volumes = np.intersect1d(internal_boundary_volumes, volumes)
                coarse_internal_boundary_volumes.append(internal_boundary_volumes)

                # cont += 1

            coarse_faces = np.array(coarse_faces)
            coarse_internal_faces = np.array(coarse_internal_faces)
            coarse_intersect_faces = np.array(coarse_intersect_faces)
            coarse_internal_boundary_volumes = np.array(coarse_internal_boundary_volumes)
            self._data[self.internal_boundary_fine_volumes+str(level)] = coarse_internal_boundary_volumes
            self._data[self.coarse_faces+str(level)] = coarse_faces
            self._data[self.coarse_intersect_faces+str(level)] = coarse_intersect_faces
            self._data[self.coarse_internal_faces+str(level)] = coarse_internal_faces

    def get_elements_2(self, M):
        assert not self._loaded
        meshsets_nv1 = self._coarse_volumes[1]
        meshsets_nv2 = self._coarse_volumes[2]
        all_meshsets = [meshsets_nv1, meshsets_nv2]

        dict_all_faces = M.data.dict_elements[direc.entities_lv0[2]]
        dict_all_volumes = M.data.dict_elements[direc.entities_lv0[3]]
        dict_all_edges = M.data.dict_elements[direc.entities_lv0[1]]
        dict_all_nodes = M.data.dict_elements[direc.entities_lv0[0]]

        mb = M.core.mb
        mtu = M.core.mtu
        n_levels = 2

        name_tag_faces_boundary_meshsets = self.faces_boundary_meshset_level

        for i in range(n_levels):
            coarse_volumes_property = []
            level = i+1
            name = name_tag_faces_boundary_meshsets + str(level)
            meshsets = all_meshsets[i]
            tag_boundary = self.tags[name]
            boundary_faces_elements = mb.tag_get_data(tag_boundary, M.core.root_set, flat=True)[0]
            boundary_faces_elements = mb.get_entities_by_handle(boundary_faces_elements)
            boundary_faces_all = np.array([dict_all_faces[f] for f in boundary_faces_elements])

            for m in meshsets:
                loc = self._coarse_volumes[i+1] == m
                primal_id = self._coarse_primal_id[level]['id'][loc][0]

                volumes_element = mb.get_entities_by_handle(m)
                volumes = np.array([dict_all_volumes[k] for k in volumes_element])
                faces_element = mtu.get_bridge_adjacencies(volumes_element, 3, 2)
                faces = np.array([dict_all_faces[f] for f in faces_element])
                edges_element = mtu.get_bridge_adjacencies(volumes_element, 3, 1)
                edges = np.array([dict_all_edges[f] for f in edges_element])
                nodes_element = mtu.get_bridge_adjacencies(volumes_element, 3, 0)
                nodes = np.array([dict_all_nodes[f] for f in nodes_element])
                boundary_faces = np.intersect1d(boundary_faces_all, faces)

                dados = {
                    'volumes': volumes, 'volumes_element': volumes_element,
                    'faces': faces, 'faces_element': faces_element,
                    'edges': edges, 'edges_element': edges_element,
                    'nodes': nodes, 'nodes_element': nodes_element,
                    'primal_id': primal_id, 'boundary_faces': boundary_faces
                }

                coarse_volume = CoarseVolume(dados)
                coarse_volumes_property.append(coarse_volume)

            self._coarse_volumes_property[level] = coarse_volumes_property

    def get_dual_structure(self):
        M = self.mesh
        mb = M.core.mb
        gids = self.data_impress['GID_0']
        dt = [('volumes', np.dtype(int)), ('dual_id', np.dtype(int)), ('primal_id', np.dtype(int))]
        all_volumes = M.core.all_volumes

        for level in range(1, self.levels):
            structure = []
            gid_level = self.data_impress['GID_'+str(level-1)]
            coarse_id_level = self.data_impress['GID_'+str(level)]
            dual_ids = self.data_impress['DUAL_'+str(level)]
            set_interns = set(gids[dual_ids==dual_ids.min()])

            while set_interns:
                intern0 = [set_interns.pop()]
                inter = M.volumes.bridge_adjacencies(intern0, 0, 3)
                dif = set(inter) - set(intern0)

                while dif & set_interns:
                    intern0 = np.setdiff1d(inter, gids[dual_ids!=dual_ids.min()])
                    try:
                        inter = np.unique(np.concatenate(M.volumes.bridge_adjacencies(intern0, 0, 3)))
                    except:
                        inter = np.unique(M.volumes.bridge_adjacencies(intern0, 0, 3))
                    dif = set(inter) - set(intern0)

                if level == 1:
                    primais1 = coarse_id_level[inter]
                    all_primal_ids = np.unique(primais1)
                    all_vertex = []
                    for gidc in all_primal_ids:
                        vertex_all = gid_level[dual_ids==3]
                        vols_in_coarse_id = gid_level[coarse_id_level==gidc]
                        vertex = np.intersect1d(vertex_all, vols_in_coarse_id)
                        all_vertex.append(vertex)

                    all_vertex = np.concatenate(all_vertex)
                    vertex_in_inter = np.intersect1d(inter, all_vertex)
                    all_vertex = np.setdiff1d(all_vertex, vertex_in_inter)
                    if len(all_vertex) > 0:
                        inter = np.concatenate([inter, all_vertex])

                _gids1 = gid_level[inter]
                _primais = coarse_id_level[inter]
                _duais = dual_ids[inter]

                if level > 1:
                    yy1 = dict(zip(_gids1, _primais))
                    yy2 = dict(zip(_gids1, _duais))
                    test1 = np.array(list(yy1.keys()))
                    test2 = np.array(list(yy2.keys()))
                    if not np.allclose(test1, test2):
                        print('erro')
                        import pdb; pdb.set_trace()
                    gids2 = test1
                    primais = np.array(list(yy1.values()))
                    duais = np.array(list(yy2.values()))
                    # gids2, duais = get_levelantids_levelids(_gids1, _duais)
                    # gids2, primais = get_levelantids_levelids(_gids1, _primais)
                else:
                    gids2 = _gids1
                    duais = _duais
                    primais = _primais

                # ####################################
                # ## teste
                # if level == 1:
                #     vertices = gids2[duais==3]
                #     if len(vertices) < 8:
                #         av = mb.create_meshset()
                #         mb.add_entities(av, all_volumes[gids2])
                #         mb.write_file('teste.vtk', [av])
                #         import pdb; pdb.set_trace()
                # else:
                #     vertices = gids2[duais==3]
                #     av = mb.create_meshset()
                #     for vv in gids2:
                #         mb.add_entities(av, all_volumes[gids[gid_level==vv]])
                #     mb.write_file('teste.vtk', [av])
                #     import pdb; pdb.set_trace()
                # ####################################

                sarray = np.zeros(len(gids2), dtype=dt)
                sarray['volumes'] = gids2
                sarray['dual_id'] = duais
                sarray['primal_id'] = primais
                structure.append(sarray)
                set_interns = set_interns - set(inter)


            self._data[self.dual_structure+str(level)] = np.array(structure)

    def get_dual_structure_with_graph(self):
        M = self.mesh
        mb = M.core.mb

        gids = self.data_impress['GID_0']
        dt = [('volumes', np.dtype(int))]
        all_volumes = M.core.all_volumes

        for level in range(1, self.levels):
            dual_flags = self.data_impress['DUAL_'+str(level)]
            primal_ids = self.data_impress['GID_'+str(level)]
            interns = gids[dual_flags==dual_flags.min()]
            intern_definitor = -np.ones(len(M.volumes.all),dtype=int)
            intern_definitor[interns]=interns

            adjs=M.faces.bridge_adjacencies(M.faces.internal,2,3)
            adjs=intern_definitor[adjs]
            adjs=adjs[(adjs>-1).sum(axis=1)==2]

            adjs0=adjs[:,0]
            adjs1=adjs[:,1]

            mapd=np.arange(len(M.volumes.all))
            mapd[interns]=np.arange(len(interns))
            lines=np.concatenate([mapd[adjs0],mapd[adjs1]])
            cols=np.concatenate([mapd[adjs1],mapd[adjs0]])
            data=np.ones(len(lines))

            from scipy.sparse import csc_matrix, csgraph
            graph=csc_matrix((data,(lines,cols)),shape=(len(interns),len(interns)))
            n_l,labels=csgraph.connected_components(graph,connection='strong')
            conjs_interns=[interns[labels==l] for l in range(n_l)]

            structure = [np.unique(np.concatenate(M.volumes.bridge_adjacencies(intern0, 0, 3))) for intern0 in conjs_interns]
            self._data[self.dual_structure+str(level)]=structure

    def set_volumes_without_gravity_source_term(self):

        lim = self.data_impress['hs'].min()*(0.2)

        # all_centroids = self.data_impress['centroid_volumes']

        # for level in range(1, self.levels):
        # fazendo apenas para o nivel 1
        for level in range(1, 2):
            structures = self._data[self.dual_structure+str(level)]
            all_vols_without_grav = []
            all_centroids = self._data[self.centroids_name + str(level-1)]
            dual_flags = self.data_impress["DUAL_"+str(level)]
            for structure in structures:
                # volumes = structure['volumes']
                # dual_id = structure['dual_id']
                volumes = structure
                dual_id = dual_flags[volumes]

                local_centroids = all_centroids[volumes]
                xmin, ymin, zmin = local_centroids.min(axis=0)
                xmax, ymax, zmax = local_centroids.max(axis=0)
                b1 = np.array([np.array([xmin-lim, ymin-lim, zmax-lim]), np.array([xmax+lim, ymax+lim, zmax+lim])])
                b2 = np.array([np.array([xmin-lim, ymin-lim, zmin-lim]), np.array([xmax+lim, ymax+lim, zmin+lim])])

                vols_without_grav = np.concatenate([volumes[get_box(local_centroids, b1)], volumes[get_box(local_centroids, b2)]])
                vols_without_grav = np.setdiff1d(vols_without_grav, volumes[dual_id==3])
                all_vols_without_grav.append(vols_without_grav)

            try:
                all_vols_without_grav = np.unique(np.concatenate(all_vols_without_grav))
            except Exception as e:
                all_vols_without_grav = np.unique(all_vols_without_grav)


            self._data[self.volumes_without_grav + str(level-1)] = all_vols_without_grav

    def save_mesh(self):
        M = self.mesh

        self.data_impress.update_variables_to_mesh()

        # M.core.print(file=self.name_mesh, config_input='input_cards/print_settings.yml')
        M.save_variables('multiscale_data')
