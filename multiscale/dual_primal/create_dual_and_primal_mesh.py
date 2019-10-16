import directories as direc
from utils.utils_old import getting_tag, get_box, Min_Max
# from pymoab import core, types, rng, topo_util
from pymoab import types, rng
import numpy as np

class DualPrimalMesh1:

    def __init__(self):
        self._loaded = False
        self.tags = dict()
        self.tags_to_infos = dict()
        self._coarse_volumes = dict()
        self._coarse_primal_id = dict()
        self._interns = dict()
        self._faces = dict()
        self._edges = dict()
        self._vertex = dict()
        self._fine_primal_id = dict()
        self._fine_dual_id = dict()
        self.mvs = dict()

    def create_tags(self, M):
        assert not self._loaded

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
            tipo = 'handle'
            entitie = 'root_set'
            t1 = types.MB_TYPE_HANDLE
            t2 = types.MB_TAG_MESH
            getting_tag(mb, name, n, t1, t2, True, entitie, tipo, self.tags, self.tags_to_infos)

        return 0

    def loaded(self):
        assert not self._loaded
        self._loaded = True

    def run(self, M):

        assert not self._loaded

        M.dualprimal = self
        self.create_tags(M)
        # self.set_primal_level_l_meshsets(M)
        self.generate_dual_and_primal(M)
        self.get_elements(M)
        self.loaded()

    def set_primal_level_l_meshsets(self, M):
        assert not self._loaded

        mb = M.core.mb
        l2_meshset = mb.create_meshset()
        mb.tag_set_data(self.tags['L2_MESHSET'], 0, l2_meshset)

        coarse_elements = M.coarse.elements
        all_volumes = np.array(M.core.all_volumes)



        import pdb; pdb.set_trace()

        for cvol in coarse_elements:
            volumes = all_volumes[cvol.volumes.global_id[:]]





        import pdb; pdb.set_trace()




        return 0

    def generate_dual_and_primal(self, M):
        assert not self._loaded

        def get_hs(M, coord_nodes):

            unis = np.array([np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])])
            nos0 = M.volumes.bridge_adjacencies(0, 2, 0)[0]
            n0 = coord_nodes[nos0[0]]
            hs = np.zeros(3)

            for i in range(1, len(nos0)):
                n1 = coord_nodes[nos0[i]]
                v = n0 - n1
                norma = np.linalg.norm(v)
                uni = np.absolute(v/norma)
                if np.allclose(uni, unis[0]):
                    hs[0] = norma
                elif np.allclose(uni, unis[1]):
                    hs[1] = norma
                elif np.allclose(uni, unis[2]):
                    hs[2] = norma

            return hs

        cr1 = direc.data_loaded['Crs']['Cr1']
        cr2 = direc.data_loaded['Crs']['Cr2']

        coord_nodes = M.data.centroids[direc.entities_lv0[0]]

        mb = M.core.mb

        Lx, Ly, Lz = coord_nodes.max(axis = 0)
        xmin, ymin, zmin = coord_nodes.min(axis = 0)
        xmax, ymax, zmax = Lx, Ly, Lz

        lx, ly, lz = get_hs(M, coord_nodes)
        dx0 = lx
        dy0 = ly
        dz0 = lz

        nx = int(Lx/lx)
        ny = int(Ly/ly)
        nz = int(Lz/lz)

        l1 = [cr1[0]*lx,cr1[1]*ly,cr1[2]*lz]
        l2 = [cr2[0]*lx,cr2[1]*ly,cr2[2]*lz]

        x1=nx*lx
        y1=ny*ly
        z1=nz*lz

        L2_meshset = mb.create_meshset()
        mb.tag_set_data(self.tags['L2_MESHSET'], 0, L2_meshset)

        lx2, ly2, lz2 = [], [], []
        # O valor 0.01 é adicionado para corrigir erros de ponto flutuante
        for i in range(int(Lx/l2[0])):    lx2.append(xmin+i*l2[0])
        for i in range(int(Ly/l2[1])):    ly2.append(ymin+i*l2[1])
        for i in range(int(Lz/l2[2])):    lz2.append(zmin+i*l2[2])
        lx2.append(Lx)
        ly2.append(Ly)
        lz2.append(Lz)

        lx1, ly1, lz1 = [], [], []
        for i in range(int(l2[0]/l1[0])):   lx1.append(i*l1[0])
        for i in range(int(l2[1]/l1[1])):   ly1.append(i*l1[1])
        for i in range(int(l2[2]/l1[2])):   lz1.append(i*l1[2])


        D_x=max(Lx-int(Lx/l1[0])*l1[0],Lx-int(Lx/l2[0])*l2[0])
        D_y=max(Ly-int(Ly/l1[1])*l1[1],Ly-int(Ly/l2[1])*l2[1])
        D_z=max(Lz-int(Lz/l1[2])*l1[2],Lz-int(Lz/l2[2])*l2[2])
        nD_x=int((D_x+0.001)/l1[0])
        nD_y=int((D_y+0.001)/l1[1])
        nD_z=int((D_z+0.001)/l1[2])

        lxd1=[xmin+dx0/100]
        for i in range(int(Lx/l1[0])-2-nD_x):
            lxd1.append(l1[0]/2+(i+1)*l1[0])
        lxd1.append(xmin+Lx-dx0/100)

        lyd1=[ymin+dy0/100]
        for i in range(int(Ly/l1[1])-2-nD_y):
            lyd1.append(l1[1]/2+(i+1)*l1[1])
        lyd1.append(ymin+Ly-dy0/100)

        lzd1=[zmin+dz0/100]

        for i in range(int(Lz/l1[2])-2-nD_z):
            lzd1.append(l1[2]/2+(i+1)*l1[2])
        lzd1.append(xmin+Lz-dz0/100)

        # print("definiu planos do nível 1")
        lxd2=[lxd1[0]]
        for i in range(1,int(len(lxd1)*l1[0]/l2[0])-1):
            lxd2.append(lxd1[int(i*l2[0]/l1[0]+0.0001)+1])
        lxd2.append(lxd1[-1])

        lyd2=[lyd1[0]]
        for i in range(1,int(len(lyd1)*l1[1]/l2[1])-1):
            lyd2.append(lyd1[int(i*l2[1]/l1[1]+0.00001)+1])
        lyd2.append(lyd1[-1])

        lzd2=[lzd1[0]]
        for i in range(1,int(len(lzd1)*l1[2]/l2[2])-1):
            lzd2.append(lzd1[int(i*l2[2]/l1[2]+0.00001)+1])
        lzd2.append(lzd1[-1])

        # print("definiu planos do nível 2")

        centroids = M.data.centroids[direc.entities_lv0[3]]
        all_volumes = np.array(M.core.all_volumes)
        dict_volumes = M.data.dict_elements[direc.entities_lv0[3]]

        D1_tag = self.tags['D1']
        D2_tag = self.tags['D2']
        primal_id_tag1 = self.tags['PRIMAL_ID_1']
        primal_id_tag2 = self.tags['PRIMAL_ID_2']
        fine_to_primal1_classic_tag = self.tags['FINE_TO_PRIMAL_CLASSIC_1']
        fine_to_primal2_classic_tag = self.tags['FINE_TO_PRIMAL_CLASSIC_2']

        nc1=0
        nc2=0

        # add_parent_child(self, parent_meshset, child_meshset, exceptions = ()):
        ##-----------------------------------------------------------------
        for i in range(len(lx2)-1):
            # t1=time.time()
            if i==len(lx2)-2:
                sx=D_x
            sy=0

            #################################################
            x0=lx2[i]
            x1=lx2[i+1]
            box_x=np.array([[x0-0.01,ymin,zmin],[x1+0.01,ymax,zmax]])
            inds_vols_x=get_box(centroids, box_x)
            vols_x = all_volumes[inds_vols_x]
            x_centroids=centroids[inds_vols_x]
            map_x_centroids = dict(zip(range(len(x_centroids)), x_centroids))
            ######################################

            for j in range(len(ly2)-1):
                if j==len(ly2)-2:
                    sy=D_y
                sz=0
                #########################
                y0=ly2[j]
                y1=ly2[j+1]
                box_y=np.array([[x0-0.01,y0-0.01,zmin],[x1+0.01,y1+0.01,zmax]])
                inds_vols_y=get_box(x_centroids, box_y)
                vols_y = vols_x[inds_vols_y]
                y_centroids=np.array([map_x_centroids[k] for k in inds_vols_y])
                map_y_centroids = dict(zip(range(len(y_centroids)), y_centroids))
                ###############
                for k in range(len(lz2)-1):
                    if k==len(lz2)-2:
                        sz=D_z
                    ########################################
                    z0=lz2[k]
                    z1=lz2[k+1]
                    # tb=time.time()
                    box_dual_1=np.array([[x0-0.01,y0-0.01,z0-0.01],[x1+0.01,y1+0.01,z1+0.01]])
                    inds_vols=get_box(y_centroids, box_dual_1)
                    vols = vols_y[inds_vols]
                    ####################
                    l2_meshset=mb.create_meshset()
                    cont=0
                    elem_por_L2=vols
                    mb.add_entities(l2_meshset,elem_por_L2)
                    ## alterar
                    # centroid_p2=np.array([self.mesh.mtu.get_average_position([np.uint64(v)]) for v in elem_por_L2])
                    inds_glob_p2 = np.array([dict_volumes[k] for k in elem_por_L2])
                    centroid_p2 = centroids[inds_glob_p2]
                    cx,cy,cz=centroid_p2[:,0],centroid_p2[:,1],centroid_p2[:,2]
                    posx=np.where(abs(cx-lxd2[i])<=l1[0]/1.9)[0]
                    posy=np.where(abs(cy-lyd2[j])<=l1[1]/1.9)[0]
                    posz=np.where(abs(cz-lzd2[k])<=l1[2]/1.9)[0]
                    f1a2v3=np.zeros(len(elem_por_L2),dtype=int)
                    f1a2v3[posx]+=1
                    f1a2v3[posy]+=1
                    f1a2v3[posz]+=1
                    mb.tag_set_data(D2_tag, elem_por_L2, f1a2v3)
                    mb.tag_set_data(fine_to_primal2_classic_tag, elem_por_L2, np.repeat(nc2,len(elem_por_L2)))
                    mb.add_parent_child(L2_meshset,l2_meshset)
                    sg=mb.get_entities_by_handle(l2_meshset)
                    # print(k, len(sg), time.time()-t1)
                    # t1=time.time()
                    mb.tag_set_data(primal_id_tag2, l2_meshset, nc2)
                    centroids_primal2=centroid_p2
                    nc2+=1
                    s1x=0
                    for m in range(len(lx1)):
                        a=int(l2[0]/l1[0])*i+m
                        if Lx-D_x==lx2[i]+lx1[m]+l1[0]:# and D_x==Lx-int(Lx/l1[0])*l1[0]:
                            s1x=D_x
                        s1y=0
                        for n in range(len(ly1)):
                            b=int(l2[1]/l1[1])*j+n
                            if Ly-D_y==ly2[j]+ly1[n]+l1[1]:# and D_y==Ly-int(Ly/l1[1])*l1[1]:
                                s1y=D_y
                            s1z=0

                            for o in range(len(lz1)):
                                c=int(l2[2]/l1[2])*k+o
                                if Lz-D_z==lz2[k]+lz1[o]+l1[2]:
                                    s1z=D_z
                                l1_meshset=mb.create_meshset()
                                box_primal1 = np.array([np.array([lx2[i]+lx1[m], ly2[j]+ly1[n], lz2[k]+lz1[o]]), np.array([lx2[i]+lx1[m]+l1[0]+s1x, ly2[j]+ly1[n]+l1[1]+s1y, lz2[k]+lz1[o]+l1[2]+s1z])])
                                # elem_por_L1 = get_box(elem_por_L2, centroids_primal2, box_primal1, False)
                                inds_elem_por_l1 = get_box(centroids_primal2, box_primal1)
                                elem_por_L1 = elem_por_L2[inds_elem_por_l1]
                                mb.add_entities(l1_meshset,elem_por_L1)
                                cont1=0
                                values_1=[]
                                for e in elem_por_L1:
                                    cont1+=1
                                    f1a2v3=0
                                    M_M=Min_Max(e, M.core)
                                    if (M_M[0]<lxd1[a] and M_M[1]>=lxd1[a]):
                                        f1a2v3+=1
                                    if (M_M[2]<lyd1[b] and M_M[3]>=lyd1[b]):
                                        f1a2v3+=1
                                    if (M_M[4]<lzd1[c] and M_M[5]>=lzd1[c]):
                                        f1a2v3+=1
                                    values_1.append(f1a2v3)
                                mb.tag_set_data(D1_tag, elem_por_L1,values_1)
                                mb.tag_set_data(fine_to_primal1_classic_tag, elem_por_L1, np.repeat(nc1,len(elem_por_L1)))
                                mb.tag_set_data(primal_id_tag1, l1_meshset, nc1)
                                nc1+=1
                                mb.add_parent_child(l2_meshset,l1_meshset)
        #-------------------------------------------------------------------------------

    def get_elements(self, M):

        mb = M.core.mb
        tags_fine = ['D', 'FINE_TO_PRIMAL_CLASSIC_']
        tags_coarse = ['PRIMAL_ID_']
        all_volumes = M.core.all_volumes
        dict_volumes = M.data.dict_elements[direc.entities_lv0[3]]

        mvs = [0]

        for i in range(2):
            n = i+1
            level = n
            name_tag_c = tags_coarse[0] + str(n)
            dual_fine_name = tags_fine[0] + str(n)
            primal_fine_name = tags_fine[1] + str(n)
            mv = mvs[i]

            # coarse_volumes = mb.get_entities_by_type_and_tag(0, types.MBENTITYSET, np.array([self.tags[name_tag_c]]), np.array([None]))
            # self._coarse_volumes[level] = coarse_volumes
            interns = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]), np.array([0]))
            interns = np.array([dict_volumes[k] for k in interns])
            self._interns[level] = interns

            edges = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]), np.array([1]))
            edges = np.array([dict_volumes[k] for k in edges])
            self._edges[level] = edges

            faces = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]), np.array([2]))
            faces = np.array([dict_volumes[k] for k in faces])
            self._faces[level] = faces

            vertex = mb.get_entities_by_type_and_tag(mv, types.MBHEX, np.array([self.tags[dual_fine_name]]), np.array([3]))
            vertexes = np.array([dict_volumes[k] for k in vertex])
            self._vertex[level] = vertexes

            coarse_volumes = []
            coarse_primal_ids = []
            for i, vert in enumerate(vertex):

                primal_id = mb.tag_get_data(self.tags[primal_fine_name], vert, flat=True)[0]
                coarse_volume = mb.get_entities_by_type_and_tag(0, types.MBENTITYSET, np.array([self.tags[name_tag_c]]), np.array([primal_id]))[0]
                elements = mb.get_entities_by_handle(coarse_volume)
                ne = len(elements)
                mb.tag_set_data(self.tags[primal_fine_name], elements, np.repeat(i, ne))
                mb.tag_set_data(self.tags[name_tag_c], coarse_volume, i)
                coarse_volumes.append(coarse_volume)
                coarse_primal_ids.append(i)

            self._coarse_volumes[level] = np.array(coarse_volumes)
            self._coarse_primal_id[level] = np.array(coarse_primal_ids)

            mv1 = mb.create_meshset()
            mb.add_entities(mv1, vertex)
            self.mvs[level] = mv1
            mvs.append(mv1)

            fine_primal_id = mb.tag_get_data(self.tags[primal_fine_name], all_volumes, flat=True)
            self._fine_primal_id[level] = fine_primal_id

            fine_dual_id = mb.tag_get_data(self.tags[dual_fine_name], all_volumes, flat=True)
            self._fine_dual_id = fine_dual_id
