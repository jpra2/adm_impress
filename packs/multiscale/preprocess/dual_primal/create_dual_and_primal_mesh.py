from .... import directories as direc
from ....utils.utils_old import getting_tag, get_box, Min_Max, add_topology
# from pymoab import core, types, rng, topo_util
from pymoab import types, rng
import numpy as np
import pdb
from ....data_class.data_manager import DataManager


class MultilevelData(DataManager):

    def __init__(self, M, data_name: str='MultilevelData.npz', load=False):
        carregar = load

        super().__init__(data_name=data_name, load=load)
        self.tags = dict()
        self.tags_to_infos = dict()
        self._carregar = carregar
        M.multilevel_data = self
        self.mesh = M
        self.levels = 3
        self.l1 = 1

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
            entitie = 'coarse_volumes_lv2'
            t1 = types.MB_TYPE_INTEGER
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
        tags1 = ['L2_MESHSET', 'local_id_internos', 'local_fac_internos']
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

    def run(self):
        M = self.mesh

        assert not self._loaded

        self.create_tags()
        # self.set_primal_level_l_meshsets(M)
        self.generate_dual_and_primal(M)
        self.get_elements(M)
        self.get_boundary_coarse_faces(M)
        # self.get_elements_2(M)
        self.loaded()

        # if not self._carregar:
        #     self.save_mesh(M)

    def generate_dual_and_primal(self, M):
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

        nx = int(Lx / lx)
        ny = int(Ly / ly)
        nz = int(Lz / lz)

        l1 = [cr1[0] * lx, cr1[1] * ly, cr1[2] * lz]
        l2 = [cr2[0] * lx, cr2[1] * ly, cr2[2] * lz]

        x1 = nx * lx
        y1 = ny * ly
        z1 = nz * lz

        L2_meshset = mb.create_meshset()
        mb.tag_set_data(self.tags['L2_MESHSET'], M.core.root_set, L2_meshset)

        lx2, ly2, lz2 = [], [], []
        # O valor 0.01 é adicionado para corrigir erros de ponto flutuante
        for i in range(int(Lx / l2[0])):    lx2.append(xmin + i * l2[0])
        for i in range(int(Ly / l2[1])):    ly2.append(ymin + i * l2[1])
        for i in range(int(Lz / l2[2])):    lz2.append(zmin + i * l2[2])
        lx2.append(Lx)
        ly2.append(Ly)
        lz2.append(Lz)

        lx1, ly1, lz1 = [], [], []
        for i in range(int(l2[0] / l1[0])):   lx1.append(i * l1[0])
        for i in range(int(l2[1] / l1[1])):   ly1.append(i * l1[1])
        for i in range(int(l2[2] / l1[2])):   lz1.append(i * l1[2])

        D_x = max(Lx - int(Lx / l1[0]) * l1[0], Lx - int(Lx / l2[0]) * l2[0])
        D_y = max(Ly - int(Ly / l1[1]) * l1[1], Ly - int(Ly / l2[1]) * l2[1])
        D_z = max(Lz - int(Lz / l1[2]) * l1[2], Lz - int(Lz / l2[2]) * l2[2])
        nD_x = int((D_x + 0.001) / l1[0])
        nD_y = int((D_y + 0.001) / l1[1])
        nD_z = int((D_z + 0.001) / l1[2])

        lxd1 = [xmin + dx0 / 100]
        for i in range(int(Lx / l1[0]) - 2 - nD_x):
            lxd1.append(l1[0] / 2 + (i + 1) * l1[0])
        lxd1.append(xmin + Lx - dx0 / 100)

        lyd1 = [ymin + dy0 / 100]
        for i in range(int(Ly / l1[1]) - 2 - nD_y):
            lyd1.append(l1[1] / 2 + (i + 1) * l1[1])
        lyd1.append(ymin + Ly - dy0 / 100)

        lzd1 = [zmin + dz0 / 100]

        for i in range(int(Lz / l1[2]) - 2 - nD_z):
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

        # add_parent_child(self, parent_meshset, child_meshset, exceptions = ()):
        ##-----------------------------------------------------------------
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
                    mb.add_parent_child(L2_meshset, l2_meshset)
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
                                mb.add_parent_child(l2_meshset, l1_meshset)
        # -------------------------------------------------------------------------------

        mv = M.core.root_set
        n_reord = 0
        interns = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags['D1']]),
                                                  np.array([0]))
        mb.tag_set_data(self.tags['reordered_id_1'], interns, np.arange(n_reord, n_reord + len(interns)))
        n_reord += len(interns)

        edges = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags['D1']]),
                                                np.array([1]))
        mb.tag_set_data(self.tags['reordered_id_1'], edges, np.arange(n_reord, n_reord + len(edges)))
        n_reord += len(edges)

        faces = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags['D1']]),
                                                np.array([2]))
        mb.tag_set_data(self.tags['reordered_id_1'], faces, np.arange(n_reord, n_reord + len(faces)))
        n_reord += len(faces)

        vertex = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags['D1']]),
                                                 np.array([3]))
        mb.tag_set_data(self.tags['reordered_id_1'], vertex, np.arange(n_reord, n_reord + len(vertex)))
        n_reord += len(vertex)

        local_id_int_tag = self.tags['local_id_internos']
        local_id_fac_tag = self.tags['local_fac_internos']
        ID_reordenado_tag = self.tags['reordered_id_1']
        mb.tag_set_data(local_id_int_tag, all_volumes,np.repeat(len(all_volumes)+1,len(all_volumes)))
        mb.tag_set_data(local_id_fac_tag, all_volumes,np.repeat(len(all_volumes)+1,len(all_volumes)))
        sgids = 0

        intern_adjs_by_dual=[]
        faces_adjs_by_dual=[]

        for i in range(len(lxd1)-1):
            x0=lxd1[i]
            x1=lxd1[i+1]
            box_x=np.array([[x0-0.01,ymin,zmin],[x1+0.01,ymax,zmax]])
            inds_vols_x = get_box(centroids, box_x)
            vols_x = all_volumes[inds_vols_x]
            x_centroids = centroids[inds_vols_x]
            map_x_centroids = dict(zip(range(len(x_centroids)), x_centroids))

            for j in range(len(lyd1)-1):
                y0=lyd1[j]
                y1=lyd1[j+1]
                box_y=np.array([[x0-0.01,y0-0.01,zmin],[x1+0.01,y1+0.01,zmax]])
                inds_vols_y = get_box(x_centroids, box_y)
                vols_y = vols_x[inds_vols_y]
                y_centroids = np.array([map_x_centroids[k] for k in inds_vols_y])
                map_y_centroids = dict(zip(range(len(y_centroids)), y_centroids))
                for k in range(len(lzd1)-1):
                    z0=lzd1[k]
                    z1=lzd1[k+1]
                    box_dual_1=np.array([[x0-0.01,y0-0.01,z0-0.01],[x1+0.01,y1+0.01,z1+0.01]])
                    inds_vols = get_box(y_centroids, box_dual_1)
                    vols = vols_y[inds_vols]
                    tipo=mb.tag_get_data(D1_tag,vols,flat=True)
                    inter=rng.Range(np.array(vols)[np.where(tipo==0)[0]])

                    mb.tag_set_data(local_id_int_tag,inter,range(len(inter)))
                    add_topology(inter,local_id_int_tag,intern_adjs_by_dual, mb, mtu, ID_reordenado_tag)


                    fac=rng.Range(np.array(vols)[np.where(tipo==1)[0]])
                    fac_centroids = np.array([map_y_centroids[k] for k in inds_vols])
                    # fac_centroids=np.array([M1.mtu.get_average_position([f]) for f in fac])

                    box_faces_x=np.array([[x0-lx/2,y0-ly/2,z0-lz/2],[x0+lx/2,y1+ly/2,z1+lz/2]])
                    box_faces_y=np.array([[x0-lx/2,y0-ly/2,z0-lz/2],[x1+lx/2,y0+ly/2,z1+lz/2]])
                    box_faces_z=np.array([[x0-lx/2,y0-ly/2,z0-lz/2],[x1+lx/2,y1+ly/2,z0+lz/2]])

                    inds_faces_x = get_box(fac_centroids, box_faces_x)
                    faces_x = fac[inds_faces_x]
                    # faces_x=get_box(fac, fac_centroids, box_faces_x, False)

                    inds_faces_y = get_box(fac_centroids, box_faces_y)
                    faces_y = fac[inds_faces_y]
                    # faces_y=get_box(fac, fac_centroids, box_faces_y, False)
                    f1=rng.unite(faces_x,faces_y)

                    inds_faces_z = get_box(fac_centroids, box_faces_z)
                    faces_z = fac[inds_faces_z]
                    # faces_z=get_box(fac, fac_centroids, box_faces_z, False)
                    f1=rng.unite(f1,faces_z)

                    if i==len(lxd1)-2:
                        box_faces_x2=np.array([[x1-lx/2,y0-ly/2,z0-lz/2],[x1+lx/2,y1+ly/2,z1+lz/2]])
                        inds_faces_x2 = get_box(fac_centroids, box_faces_x2)
                        faces_x2 = fac[inds_faces_x2]
                        # faces_x2=get_box(fac, fac_centroids, box_faces_x2, False)
                        f1=rng.unite(f1,faces_x2)

                    if j==len(lyd1)-2:
                        box_faces_y2=np.array([[x0-lx/2,y1-ly/2,z0-lz/2],[x1+lx/2,y1+ly/2,z1+lz/2]])
                        inds_faces_y2 = get_box(fac_centroids, box_faces_y2)
                        faces_y2 = fac[inds_faces_y2]
                        # faces_y2=get_box(fac, fac_centroids, box_faces_y2, False)
                        f1=rng.unite(f1,faces_y2)

                    if k==len(lzd1)-2:
                        box_faces_z2=np.array([[x0-lx/2,y0-ly/2,z1-lz/2],[x1+lx/2,y1+ly/2,z1+lz/2]])
                        inds_faces_z2 = get_box(fac_centroids, box_faces_z2)
                        faces_z2 = fac[inds_faces_z2]
                        # faces_z2=get_box(fac, fac_centroids, box_faces_z2, False)
                        f1=rng.unite(f1,faces_z2)

                    sgids+=len(f1)
                    mb.tag_set_data(local_id_fac_tag,f1,range(len(f1)))
                    add_topology(f1,local_id_fac_tag,faces_adjs_by_dual, mb, mtu, ID_reordenado_tag)

        self._data['intern_adjs_by_dual'] = np.array(intern_adjs_by_dual)
        self._data['faces_adjs_by_dual'] = np.array(faces_adjs_by_dual)

    def get_elements(self, M):
        assert not self._loaded

        mb = M.core.mb
        tags_fine = ['D', 'FINE_TO_PRIMAL_CLASSIC_']
        tags_coarse = ['PRIMAL_ID_']
        tag_mv = ['MV_']
        tag_reordered_id = ['reordered_id_']
        all_volumes = M.core.all_volumes
        dict_volumes = dict(zip(all_volumes, M.volumes.all))

        mvs = [M.core.root_set]

        for i in range(2):
            n = i + 1
            level = n
            name_tag_c = tags_coarse[0] + str(n)
            dual_fine_name = tags_fine[0] + str(n)
            primal_fine_name = tags_fine[1] + str(n)
            tag_reord_id = tag_reordered_id[0] + str(n)
            mv = mvs[i]
            n_reord = 0

            # coarse_volumes = mb.get_entities_by_type_and_tag(0, types.MBENTITYSET, np.array([self.tags[name_tag_c]]), np.array([None]))
            # self._coarse_volumes[level] = coarse_volumes

            interns = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]),
                                                      np.array([0]))
            if n not in [1]:
                mb.tag_set_data(self.tags[tag_reord_id], interns, np.arange(n_reord, len(interns)))
            n_reord += len(interns)
            interns = np.array([dict_volumes[k] for k in interns])
            # self._interns[level] = interns

            edges = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]),
                                                    np.array([1]))
            if n not in [1]:
                mb.tag_set_data(self.tags[tag_reord_id], edges, np.arange(n_reord, n_reord + len(edges)))
            n_reord += len(edges)
            edges = np.array([dict_volumes[k] for k in edges])
            # self._edges[level] = edges

            faces = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]),
                                                    np.array([2]))
            if n not in [1]:
                mb.tag_set_data(self.tags[tag_reord_id], faces, np.arange(n_reord, n_reord + len(faces)))
            n_reord += len(faces)
            faces = np.array([dict_volumes[k] for k in faces])
            # self._faces[level] = faces

            vertex = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]),
                                                     np.array([3]))
            if n not in [1]:
                mb.tag_set_data(self.tags[tag_reord_id], vertex, np.arange(n_reord, n_reord + len(vertex)))
            n_reord += len(vertex)
            vertexes = np.array([dict_volumes[k] for k in vertex])
            # self._vertex[level] = vertexes

            self._data[self.interns + str(level)] = interns
            self._data[self.faces + str(level)] = faces
            self._data[self.edges + str(level)] = edges
            self._data[self.vertex + str(level)] = vertexes

            coarse_volumes = []
            coarse_primal_ids = []
            if not self._carregar:
                for i, vert in enumerate(vertex):
                    primal_id = mb.tag_get_data(self.tags[primal_fine_name], vert, flat=True)[0]
                    coarse_volume = \
                    mb.get_entities_by_type_and_tag(0, types.MBENTITYSET, np.array([self.tags[name_tag_c]]),
                                                    np.array([primal_id]))[0]
                    # elements = mb.get_entities_by_handle(coarse_volume)
                    # ne = len(elements)
                    # mb.tag_set_data(self.tags[primal_fine_name], elements, np.repeat(i, ne))
                    # mb.tag_set_data(self.tags[name_tag_c], coarse_volume, i)
                    coarse_volumes.append(coarse_volume)
                    coarse_primal_ids.append(primal_id)
            else:
                for i, vert in enumerate(vertex):
                    primal_id = mb.tag_get_data(self.tags[primal_fine_name], vert, flat=True)[0]
                    coarse_volume = \
                    mb.get_entities_by_type_and_tag(0, types.MBENTITYSET, np.array([self.tags[name_tag_c]]),
                                                    np.array([primal_id]))[0]
                    # elements = mb.get_entities_by_handle(coarse_volume)
                    # ne = len(elements)
                    # mb.tag_set_data(self.tags[primal_fine_name], elements, np.repeat(i, ne))
                    # mb.tag_set_data(self.tags[name_tag_c], coarse_volume, i)
                    coarse_volumes.append(coarse_volume)
                    coarse_primal_ids.append(primal_id)

            # self._coarse_volumes[level] = np.array(coarse_volumes)
            self._data[self.coarse_volumes + str(level)] = np.array(coarse_volumes)
            dtype = [('elements', np.uint64), ('id', np.uint64)]
            structured_array = np.zeros(len(coarse_volumes), dtype=dtype)
            structured_array['elements'] = np.array(coarse_volumes)
            structured_array['id'] = np.array(coarse_primal_ids)
            # self._coarse_primal_id[level] = structured_array
            self._data[self.coarse_primal_id + str(level)] = np.array(coarse_primal_ids)

            nnn = tag_mv[0] + str(level)
            if not self._carregar:
                mv1 = mb.create_meshset()
                mb.add_entities(mv1, vertex)
                mb.tag_set_data(self.tags[nnn], M.core.root_set, mv1)
            else:
                mv1 = mb.tag_get_data(self.tags[nnn], M.core.root_set, flat=True)[0]

            # self.mvs[level] = mv1
            mvs.append(mv1)
            self._data[self.meshset_vertices + str(level)] = mv1

            fine_primal_id = mb.tag_get_data(self.tags[primal_fine_name], all_volumes, flat=True)
            # self._fine_primal_id[level] = fine_primal_id
            self._data[self.fine_primal_id + str(level)] = fine_primal_id

            fine_dual_id = mb.tag_get_data(self.tags[dual_fine_name], all_volumes, flat=True)
            # self._fine_dual_id[level] = fine_dual_id
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

        self._data[self.reordered_id + str(1)] = mb.tag_get_data(self.tags[self.reordered_id + str(1)], all_volumes, flat=True)
        self._data[self.reordered_id + str(2)] = mb.tag_get_data(self.tags[self.reordered_id + str(2)], all_volumes, flat=True)

    def get_boundary_coarse_faces(self, M):
        assert not self._loaded
        # meshsets_nv1 = self._coarse_volumes[1]
        # meshsets_nv2 = self._coarse_volumes[2]
        meshsets_nv1 = self._data[self.coarse_volumes + str(1)]
        meshsets_nv2 = self._data[self.coarse_volumes + str(2)]

        mb = M.core.mb
        mtu = M.core.mtu
        n_levels = 2

        name_tag_faces_boundary_meshsets = self.faces_boundary_meshset_level
        all_meshsets = [meshsets_nv1, meshsets_nv2]
        d_faces = dict(zip(M.core.all_faces, M.faces.all))

        from ....utils import pymoab_utils as utpy

        for i in range(n_levels):
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
            self._data[self.faces_boundary_meshset_level + str(i+1)] = np.array([d_faces[k] for k in faces_boundary])

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

    def save_mesh_dep0(self, M):

        M.state = 2
        np.save(direc.state_path, np.array([M.state]))
        np.save(direc.path_local_last_file_name, np.array([direc.names_outfiles_steps[2]]))
        M.core.print(file=direc.output_file + str(M.state))

    def save_mesh(self):
        M = self.mesh

        M.core.print(file=self.name_mesh, config_input='input_cards/print_settings.yml')





'''
def propriedade():
    doc = "The propriedade property."

    def fget(self):
        return self._propriedade

    def fset(self, value):
        self._propriedade = value

    def fdel(self):
        del self._propriedade

    return locals()

propriedade = property(**propriedade())
'''
