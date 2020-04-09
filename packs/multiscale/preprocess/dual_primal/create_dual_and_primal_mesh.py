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
from .paralell import create_dual_and_primal
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
        self.levels = 3
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

        self.create_tags()
        # self.set_primal_level_l_meshsets(M)
        self.generate_dual_and_primal_any_D(M)
        self.get_elements(M)
        self.get_boundary_coarse_faces(M)
        self.get_dual_structure()
        self.set_volumes_without_gravity_source_term()
        # self.get_elements_2(M)
        self.export_to_npz()
        self.loaded()

        # if not self._carregar:
        #     self.save_mesh(M)
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

    def generate_dual_and_primal_3D(self, M):
        assert not self._loaded

        def get_hs(M, coord_nodes):

            unis = np.array([np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])])
            nos0 = M.volumes.bridge_adjacencies(0, 2, 0)[0]
            n0 = coord_nodes[nos0[0]]
            hs = np.zeros(3)

            for i in range(1, len(nos0)):
                n1 = coord_nodes[nos0[i]]
                v = n0 - n1
                norma = np.linalg.norm(v)
                uni = np.absolute(v / norma)
                if np.allclose(uni, unis[0]):
                    hs[0] = norma
                elif np.allclose(uni, unis[1]):
                    hs[1] = norma
                elif np.allclose(uni, unis[2]):
                    hs[2] = norma

            return hs

        cr1 = direc.data_loaded['Crs']['Cr1']
        cr2 = direc.data_loaded['Crs']['Cr2']

        coord_nodes = M.data['centroid_nodes']

        mb = M.core.mb
        mtu = M.core.mtu

        Lx, Ly, Lz = coord_nodes.max(axis=0)
        xmin, ymin, zmin = coord_nodes.min(axis=0)
        xmax, ymax, zmax = Lx, Ly, Lz

        # lx, ly, lz = get_hs(M, coord_nodes)
        lx, ly, lz = M.data['hs'][0]

        dx0 = lx
        dy0 = ly
        dz0 = lz

        nx = int(round(Lx / lx))
        ny = int(round(Ly / ly))
        nz = int(round(Lz / lz))

        l1 = [cr1[0] * lx, cr1[1] * ly, cr1[2] * lz]
        l2 = [cr2[0] * lx, cr2[1] * ly, cr2[2] * lz]

        x1 = nx * lx
        y1 = ny * ly
        z1 = nz * lz

        L2_meshset = mb.create_meshset()
        mb.tag_set_data(self.tags['L2_MESHSET'], M.core.root_set, L2_meshset)

        lx2, ly2, lz2 = [], [], []
        # O valor 0.01 é adicionado para corrigir erros de ponto flutuante
        for i in range(int(round(Lx / l2[0]))):    lx2.append(xmin + i * l2[0])
        for i in range(int(round(Ly / l2[1]))):    ly2.append(ymin + i * l2[1])
        for i in range(int(round(Lz / l2[2]))):    lz2.append(zmin + i * l2[2])
        lx2.append(Lx)
        ly2.append(Ly)
        lz2.append(Lz)

        lx1, ly1, lz1 = [], [], []
        for i in range(int(round(l2[0] / l1[0]))):   lx1.append(i * l1[0])
        for i in range(int(round(l2[1] / l1[1]))):   ly1.append(i * l1[1])
        for i in range(int(round(l2[2] / l1[2]))):   lz1.append(i * l1[2])

        D_x = max(Lx - int(round(Lx / l1[0])) * l1[0], Lx - int(round(Lx / l2[0])) * l2[0])
        D_y = max(Ly - int(round(Ly / l1[1])) * l1[1], Ly - int(round(Ly / l2[1])) * l2[1])
        D_z = max(Lz - int(round(Lz / l1[2])) * l1[2], Lz - int(round(Lz / l2[2])) * l2[2])
        nD_x = int((D_x + 0.001) / l1[0])
        nD_y = int((D_y + 0.001) / l1[1])
        nD_z = int((D_z + 0.001) / l1[2])

        lxd1 = [xmin + dx0 / 100]
        for i in range(int(round(Lx / l1[0])) - 2 - nD_x):
            lxd1.append(l1[0] / 2 + (i + 1) * l1[0])
        lxd1.append(xmin + Lx - dx0 / 100)

        lyd1 = [ymin + dy0 / 100]
        for i in range(int(round(Ly / l1[1])) - 2 - nD_y):
            lyd1.append(l1[1] / 2 + (i + 1) * l1[1])
        lyd1.append(ymin + Ly - dy0 / 100)

        lzd1 = [zmin + dz0 / 100]

        for i in range(int(round(Lz / l1[2])) - 2 - nD_z):
            lzd1.append(l1[2] / 2 + (i + 1) * l1[2])
        lzd1.append(xmin + Lz - dz0 / 100)

        # print("definiu planos do nível 1")
        lxd2 = [lxd1[0]]
        for i in range(1, int(len(lxd1) * l1[0] / l2[0]) - 1):
            lxd2.append(lxd1[int(i * l2[0] / l1[0] + 0.0001) + 1])
        lxd2.append(lxd1[-1])

        lyd2 = [lyd1[0]]
        for i in range(1, int(len(lyd1) * l1[1] / l2[1]) - 1):
            lyd2.append(lyd1[int(i * l2[1] / l1[1] + 0.00001) + 1])
        lyd2.append(lyd1[-1])

        lzd2 = [lzd1[0]]
        for i in range(1, int(len(lzd1) * l1[2] / l2[2]) - 1):
            lzd2.append(lzd1[int(i * l2[2] / l1[2] + 0.00001) + 1])
        lzd2.append(lzd1[-1])

        # print("definiu planos do nível 2")

        centroids = M.data['centroid_volumes']
        all_volumes = np.array(M.core.all_volumes)
        dict_volumes = dict(zip(all_volumes, M.volumes.all))

        D1_tag = self.tags['D1']
        D2_tag = self.tags['D2']
        primal_id_tag1 = self.tags['PRIMAL_ID_1']
        primal_id_tag2 = self.tags['PRIMAL_ID_2']
        fine_to_primal1_classic_tag = self.tags['FINE_TO_PRIMAL_CLASSIC_1']
        fine_to_primal2_classic_tag = self.tags['FINE_TO_PRIMAL_CLASSIC_2']

        nc1 = 0
        nc2 = 0

        for i in range(len(lx2) - 1):
            # t1=time.time()
            if i == len(lx2) - 2:
                sx = D_x
            sy = 0

            #################################################
            x0 = lx2[i]
            x1 = lx2[i + 1]
            box_x = np.array([[x0 - 0.01, ymin, zmin], [x1 + 0.01, ymax, zmax]])
            inds_vols_x = get_box(centroids, box_x)
            vols_x = all_volumes[inds_vols_x]
            x_centroids = centroids[inds_vols_x]
            map_x_centroids = dict(zip(range(len(x_centroids)), x_centroids))
            ######################################

            for j in range(len(ly2) - 1):
                if j == len(ly2) - 2:
                    sy = D_y
                sz = 0
                #########################
                y0 = ly2[j]
                y1 = ly2[j + 1]
                box_y = np.array([[x0 - 0.01, y0 - 0.01, zmin], [x1 + 0.01, y1 + 0.01, zmax]])
                inds_vols_y = get_box(x_centroids, box_y)
                vols_y = vols_x[inds_vols_y]
                y_centroids = np.array([map_x_centroids[k] for k in inds_vols_y])
                map_y_centroids = dict(zip(range(len(y_centroids)), y_centroids))
                ###############
                for k in range(len(lz2) - 1):
                    if k == len(lz2) - 2:
                        sz = D_z
                    ########################################
                    z0 = lz2[k]
                    z1 = lz2[k + 1]
                    # tb=time.time()
                    box_dual_1 = np.array([[x0 - 0.01, y0 - 0.01, z0 - 0.01], [x1 + 0.01, y1 + 0.01, z1 + 0.01]])
                    inds_vols = get_box(y_centroids, box_dual_1)
                    vols = vols_y[inds_vols]
                    ####################
                    l2_meshset = mb.create_meshset()
                    cont = 0
                    elem_por_L2 = vols
                    mb.add_entities(l2_meshset, elem_por_L2)
                    ## alterar
                    # centroid_p2=np.array([self.mesh.mtu.get_average_position([np.uint64(v)]) for v in elem_por_L2])
                    inds_glob_p2 = np.array([dict_volumes[k] for k in elem_por_L2])
                    centroid_p2 = centroids[inds_glob_p2]
                    cx, cy, cz = centroid_p2[:, 0], centroid_p2[:, 1], centroid_p2[:, 2]
                    try:
                        posx = np.where(abs(cx - lxd2[i]) <= l1[0] / 1.9)[0]
                        posy = np.where(abs(cy - lyd2[j]) <= l1[1] / 1.9)[0]
                        posz = np.where(abs(cz - lzd2[k]) <= l1[2] / 1.9)[0]
                    except:
                        import pdb; pdb.set_trace()
                    f1a2v3 = np.zeros(len(elem_por_L2), dtype=int)
                    f1a2v3[posx] += 1
                    f1a2v3[posy] += 1
                    f1a2v3[posz] += 1
                    mb.tag_set_data(D2_tag, elem_por_L2, f1a2v3)
                    mb.tag_set_data(fine_to_primal2_classic_tag, elem_por_L2, np.repeat(nc2, len(elem_por_L2)))
                    # mb.add_parent_child(L2_meshset, l2_meshset)
                    # sg = mb.get_entities_by_handle(l2_meshset)
                    # print(k, len(sg), time.time()-t1)
                    # t1=time.time()
                    mb.tag_set_data(primal_id_tag2, l2_meshset, nc2)
                    centroids_primal2 = centroid_p2
                    nc2 += 1
                    s1x = 0
                    for m in range(len(lx1)):
                        a = int(l2[0] / l1[0]) * i + m
                        if Lx - D_x == lx2[i] + lx1[m] + l1[0]:  # and D_x==Lx-int(Lx/l1[0])*l1[0]:
                            s1x = D_x
                        s1y = 0
                        for n in range(len(ly1)):
                            b = int(l2[1] / l1[1]) * j + n
                            if Ly - D_y == ly2[j] + ly1[n] + l1[1]:  # and D_y==Ly-int(Ly/l1[1])*l1[1]:
                                s1y = D_y
                            s1z = 0

                            for o in range(len(lz1)):
                                c = int(l2[2] / l1[2]) * k + o
                                if Lz - D_z == lz2[k] + lz1[o] + l1[2]:
                                    s1z = D_z
                                l1_meshset = mb.create_meshset()
                                box_primal1 = np.array([np.array([lx2[i] + lx1[m], ly2[j] + ly1[n], lz2[k] + lz1[o]]),
                                                        np.array([lx2[i] + lx1[m] + l1[0] + s1x,
                                                                  ly2[j] + ly1[n] + l1[1] + s1y,
                                                                  lz2[k] + lz1[o] + l1[2] + s1z])])
                                # elem_por_L1 = get_box(elem_por_L2, centroids_primal2, box_primal1, False)
                                inds_elem_por_l1 = get_box(centroids_primal2, box_primal1)
                                elem_por_L1 = elem_por_L2[inds_elem_por_l1]
                                mb.add_entities(l1_meshset, elem_por_L1)
                                cont1 = 0
                                values_1 = []
                                for e in elem_por_L1:
                                    cont1 += 1
                                    f1a2v3 = 0
                                    M_M = Min_Max(e, M.core)
                                    if (M_M[0] < lxd1[a] and M_M[1] >= lxd1[a]):
                                        f1a2v3 += 1
                                    if (M_M[2] < lyd1[b] and M_M[3] >= lyd1[b]):
                                        f1a2v3 += 1
                                    if (M_M[4] < lzd1[c] and M_M[5] >= lzd1[c]):
                                        f1a2v3 += 1
                                    values_1.append(f1a2v3)
                                mb.tag_set_data(D1_tag, elem_por_L1, values_1)
                                mb.tag_set_data(fine_to_primal1_classic_tag, elem_por_L1,
                                                np.repeat(nc1, len(elem_por_L1)))
                                mb.tag_set_data(primal_id_tag1, l1_meshset, nc1)
                                nc1 += 1
                                # mb.add_parent_child(l2_meshset, l1_meshset)
    def generate_dual_and_primal_any_D(self, M):
        M1=M.core
        M1.all_centroids=M.data["centroid_volumes"]
        M1.primal_id_tag1=self.tags["PRIMAL_ID_1"]
        M1.primal_id_tag2=self.tags["PRIMAL_ID_2"]
        M1.fine_to_primal1_classic_tag=self.tags["FINE_TO_PRIMAL_CLASSIC_1"]
        M1.fine_to_primal2_classic_tag=self.tags["FINE_TO_PRIMAL_CLASSIC_2"]
        M1.D1_tag=self.tags["D1"]
        M1.D2_tag=self.tags["D2"]
        def get_hs(M, coord_nodes):

            unis = np.array([np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])])
            nos0 = M.volumes.bridge_adjacencies(0, 2, 0)[0]
            n0 = coord_nodes[nos0[0]]
            hs = np.zeros(3)

            for i in range(1, len(nos0)):
                n1 = coord_nodes[nos0[i]]
                v = n0 - n1
                norma = np.linalg.norm(v)
                uni = np.absolute(v / norma)
                if np.allclose(uni, unis[0]):
                    hs[0] = norma
                elif np.allclose(uni, unis[1]):
                    hs[1] = norma
                elif np.allclose(uni, unis[2]):
                    hs[2] = norma

            return hs

        cr1 = direc.data_loaded['Crs']['Cr1']
        cr2 = direc.data_loaded['Crs']['Cr2']

        coord_nodes = M.data['centroid_nodes']

        mb = M.core.mb
        mtu = M.core.mtu

        Lx, Ly, Lz = coord_nodes.max(axis=0)
        xmin, ymin, zmin = coord_nodes.min(axis=0)
        xmax, ymax, zmax = Lx, Ly, Lz

        # lx, ly, lz = get_hs(M, coord_nodes)
        lx, ly, lz = M.data['hs'][0]

        dx0 = lx
        dy0 = ly
        dz0 = lz

        nx = int(round(Lx / lx))
        ny = int(round(Ly / ly))
        nz = int(round(Lz / lz))

        l1 = [cr1[0] * lx, cr1[1] * ly, cr1[2] * lz]
        l2 = [cr2[0] * lx, cr2[1] * ly, cr2[2] * lz]

        x1 = nx * lx
        y1 = ny * ly
        z1 = nz * lz

        L2_meshset = mb.create_meshset()
        mb.tag_set_data(self.tags['L2_MESHSET'], M.core.root_set, L2_meshset)

        lx2, ly2, lz2 = [], [], []
        # O valor 0.01 é adicionado para corrigir erros de ponto flutuante
        for i in range(int(round(Lx / l2[0]))):    lx2.append(xmin + i * l2[0])
        for i in range(int(round(Ly / l2[1]))):    ly2.append(ymin + i * l2[1])
        for i in range(int(round(Lz / l2[2]))):    lz2.append(zmin + i * l2[2])
        lx2.append(Lx)
        ly2.append(Ly)
        lz2.append(Lz)

        lx1, ly1, lz1 = [], [], []
        for i in range(int(round(l2[0] / l1[0]))):   lx1.append(i * l1[0])
        for i in range(int(round(l2[1] / l1[1]))):   ly1.append(i * l1[1])
        for i in range(int(round(l2[2] / l1[2]))):   lz1.append(i * l1[2])

        D_x = max(Lx - int(round(Lx / l1[0])) * l1[0], Lx - int(round(Lx / l2[0])) * l2[0])
        D_y = max(Ly - int(round(Ly / l1[1])) * l1[1], Ly - int(round(Ly / l2[1])) * l2[1])
        D_z = max(Lz - int(round(Lz / l1[2])) * l1[2], Lz - int(round(Lz / l2[2])) * l2[2])
        nD_x = int((D_x + 0.001) / l1[0])
        nD_y = int((D_y + 0.001) / l1[1])
        nD_z = int((D_z + 0.001) / l1[2])

        lxd1 = [xmin + dx0 / 100]
        for i in range(int(round(Lx / l1[0])) - 2 - nD_x):
            lxd1.append(l1[0] / 2 + (i + 1) * l1[0])
        lxd1.append(xmin + Lx - dx0 / 100)

        lyd1 = [ymin + dy0 / 100]
        for i in range(int(round(Ly / l1[1])) - 2 - nD_y):
            lyd1.append(l1[1] / 2 + (i + 1) * l1[1])
        lyd1.append(ymin + Ly - dy0 / 100)

        lzd1 = [zmin + dz0 / 100]

        for i in range(int(round(Lz / l1[2])) - 2 - nD_z):
            lzd1.append(l1[2] / 2 + (i + 1) * l1[2])
        lzd1.append(xmin + Lz - dz0 / 100)

        # print("definiu planos do nível 1")
        lxd2 = [lxd1[0]]
        for i in range(1, int(len(lxd1) * l1[0] / l2[0]) - 1):
            lxd2.append(lxd1[int(i * l2[0] / l1[0] + 0.0001) + 1])
        lxd2.append(lxd1[-1])

        lyd2 = [lyd1[0]]
        for i in range(1, int(len(lyd1) * l1[1] / l2[1]) - 1):
            lyd2.append(lyd1[int(i * l2[1] / l1[1] + 0.00001) + 1])
        lyd2.append(lyd1[-1])

        lzd2 = [lzd1[0]]
        for i in range(1, int(len(lzd1) * l1[2] / l2[2]) - 1):
            lzd2.append(lzd1[int(i * l2[2] / l1[2] + 0.00001) + 1])
        lzd2.append(lzd1[-1])

        D_x=max(Lx-int(Lx/l1[0])*l1[0],Lx-int(Lx/l2[0])*l2[0])
        D_y=max(Ly-int(Ly/l1[1])*l1[1],Ly-int(Ly/l2[1])*l2[1])
        D_z=max(Lz-int(Lz/l1[2])*l1[2],Lz-int(Lz/l2[2])*l2[2])
        class input_dual():
            def __init__(self):
                self.l = l1, l2
                self.L = Lx, Ly, Lz
                self.d0 = dx0, dy0, dz0
                self.ncs = 0, 0
                self.lims = xmin, xmax, ymin, ymax, zmin, zmax
                self.meshsets = L2_meshset
                self.D = D_x, D_y, D_z

        class partition():
            def __init__(self):
                self.l1 = lx1, ly1, lz1
                self.ld1 = lxd1, lyd1, lzd1
                self.l2 = lx2, ly2, lz2
                self.ld2 = lxd2, lyd2, lzd2

        class SubDomain(input_dual, partition):
            def __init__(self, input_dual, partition):
                self.input_dual_and_primal = input_dual
                self.partition = partition
                self.M1 = M1
        dual_primal = create_dual_and_primal.paralell_dual_and_primal(SubDomain(input_dual(),partition()),nworker=10,first=True)
        primal_1=dual_primal.v1
        primal_2=dual_primal.v2
        prim1=[]
        for p2 in primal_1:
            for p1 in p2:
                for p0 in p1:
                    id1=M1.mb.tag_get_data(M1.fine_to_primal1_classic_tag,p0[0],flat=True)
                    ms=M1.mb.create_meshset()
                    M1.mb.add_entities(ms,p0)
                    M1.mb.tag_set_data(M1.primal_id_tag1,ms,id1)

                v2=np.concatenate(p1)
                id2=M1.mb.tag_get_data(M1.fine_to_primal2_classic_tag,v2[0],flat=True)
                ms=M1.mb.create_meshset()
                M1.mb.add_entities(ms,v2)
                M1.mb.tag_set_data(M1.primal_id_tag2,ms,id2)


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
            self._data[self.neig_intersect_faces+str(level)] = M.faces.bridge_adjacencies(faces_boundary, 2, 3)


            coarse_faces = []
            coarse_internal_faces = []
            coarse_intersect_faces = []
            coarse_internal_boundary_volumes = []

            # tag_coarse_id = self.tags['PRIMAL_ID_' + str(level)]
            # cids = self._data[self.coarse_primal_id + str(level)]
            # cont = 0

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
                    inter = np.unique(np.concatenate(M.volumes.bridge_adjacencies(intern0, 0, 3)))
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

    def set_volumes_without_gravity_source_term(self):

        lim = self.data_impress['hs'].min()*(0.2)

        # all_centroids = self.data_impress['centroid_volumes']

        # for level in range(1, self.levels):
        # fazendo apenas para o nivel 1
        for level in range(1, 2):
            structures = self._data[self.dual_structure+str(level)]
            all_vols_without_grav = []
            all_centroids = self._data[self.centroids_name + str(level-1)]

            for structure in structures:
                volumes = structure['volumes']
                dual_id = structure['dual_id']
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
