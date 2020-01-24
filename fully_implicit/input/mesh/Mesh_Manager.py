import numpy as np
from math import pi, sqrt
from pymoab import core, types, rng, topo_util
import time
import os

class MeshManager:
    def __init__(self,mesh_file, dim=2):

        self.dimension = dim
        self.mb = core.Core()
        self.root_set = self.mb.get_root_set()
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        # import pdb; pdb.set_trace()
        self.mb.load_file(mesh_file)
        # self.mb.load_file('mesh/6x6.msh')

        self.physical_tag = self.mb.tag_get_handle("MATERIAL_SET")
        self.physical_sets = self.mb.get_entities_by_type_and_tag(
            0, types.MBENTITYSET, np.array(
            (self.physical_tag,)), np.array((None,)))

        self.dirichlet_tag = self.mb.tag_get_handle(
            "Dirichlet", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)

        self.neumann_tag = self.mb.tag_get_handle(
            "Neumann", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)

        #self.perm_tag = self.mb.tag_get_handle(
        #    "Permeability", 9, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)

        self.source_tag = self.mb.tag_get_handle(
            "Source term", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)


        self.all_volumes = self.mb.get_entities_by_dimension(0, self.dimension)
        self.all_nodes = self.mb.get_entities_by_dimension(0, 0)

        self.mtu.construct_aentities(self.all_nodes)
        self.all_faces = self.mb.get_entities_by_dimension(0, self.dimension-1) #ADJs=np.array([M1.mb.get_adjacencies(face, 3) for face in M1.all_faces])

        self.all_edges = self.mb.get_entities_by_dimension(0, 1)
        self.dirichlet_faces = set()
        self.neumann_faces = set()

        '''self.GLOBAL_ID_tag = self.mb.tag_get_handle(
            "Global_ID", 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, True)'''

        self.create_tags()
        #self.set_k_and_phi_structured_spe10()
        #self.set_information("PERM", self.all_volumes, 3)
        #self.get_boundary_faces()
        self.gravity = False
        self.gama = 10
        self.mi = 1
        t0=time.time()
        # print('set area')
        #self.get_kequiv_by_face_quad(self.all_faces)
        '''
        print('set area')
        for f in self.all_faces:
            self.set_area(f)'''
        # print("took",time.time()-t0)
        self.get_faces_boundary

    def create_tags(self):
        self.perm_tag = self.mb.tag_get_handle("PERM", 9, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.finos_tag = self.mb.tag_get_handle("finos", 1, types.MB_TYPE_HANDLE, types.MB_TAG_MESH, True)
        self.Swns_tag = self.mb.tag_get_handle("Swns", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.Swn1s_tag = self.mb.tag_get_handle("Swn1s", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.wells_dirichlet_tag = self.mb.tag_get_handle("WELLS_D", 1, types.MB_TYPE_HANDLE, types.MB_TAG_MESH, True)
        self.wells_neumann_tag = self.mb.tag_get_handle("WELLS_N", 1, types.MB_TYPE_HANDLE, types.MB_TAG_MESH, True)
        self.press_value_tag = self.mb.tag_get_handle("P", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.vazao_value_tag = self.mb.tag_get_handle("Q", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.all_faces_boundary_tag = self.mb.tag_get_handle("FACES_BOUNDARY", 1, types.MB_TYPE_HANDLE, types.MB_TAG_MESH, True)
        self.area_tag = self.mb.tag_get_handle("AREA", 3, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.GLOBAL_ID_tag = self.mb.tag_get_handle("G_ID_tag", 1, types.MB_TYPE_INTEGER, types.MB_TAG_SPARSE, True)
        self.ID_reordenado_tag = self.mb.tag_get_handle("ID_reord_tag", 1, types.MB_TYPE_INTEGER, types.MB_TAG_SPARSE, True)
        self.phi_tag = self.mb.tag_get_handle("PHI", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.k_eq_tag = self.mb.tag_get_handle("K_EQ", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.n_pressure_tag = self.mb.tag_get_handle("Pressure", 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)


    def create_vertices(self, coords):
        new_vertices = self.mb.create_vertices(coords)
        self.all_nodes.append(new_vertices)
        return new_vertices

    def create_element(self, poly_type, vertices):
        new_volume = self.mb.create_element(poly_type, vertices)
        self.all_volumes.append(new_volume)
        return new_volume

    def set_information(self, information_name, physicals_values,
                        dim_target, set_connect=False):
        information_tag = self.mb.tag_get_handle(information_name)
        for physical_value in physicals_values:
            for a_set in self.physical_sets:
                physical_group = self.mb.tag_get_data(self.physical_tag,
                                                      a_set, flat=True)

                if physical_group == physical:
                    group_elements = self.mb.get_entities_by_dimension(a_set, dim_target)

                    if information_name == 'Dirichlet':
                        # print('DIR GROUP', len(group_elements), group_elements)
                        self.dirichlet_faces = self.dirichlet_faces | set(
                                                    group_elements)

                    if information_name == 'Neumann':
                        # print('NEU GROUP', len(group_elements), group_elements)
                        self.neumann_faces = self.neumann_faces | set(
                                                  group_elements)

                    for element in group_elements:
                        self.mb.tag_set_data(information_tag, element, value)

                        if set_connect:
                            connectivities = self.mtu.get_bridge_adjacencies(
                                                                element, 0, 0)
                            self.mb.tag_set_data(
                                information_tag, connectivities,
                                np.repeat(value, len(connectivities)))

    def set_k(self):
        k = 1.0
        perm_tensor = [k, 0, 0,
                       0, k, 0,
                       0, 0, k]
        for v in self.all_volumes:
            self.mb.tag_set_data(self.perm_tag, v, perm_tensor)
            #v_tags=self.mb.tag_get_tags_on_entity(v)
            #print(self.mb.tag_get_data(v_tags[1],v,flat=True))

    def set_area(self, face):
        points = self.mtu.get_bridge_adjacencies(face, 2, 0)
        points = [self.mb.get_coords([vert]) for vert in points]
        if len(points) == 3:
            n1 = np.array(points[0] - points[1])
            n2 = np.array(points[0] - points[2])
            area = np.cross(n1, n2)/2.0

        #calculo da area para quadrilatero regular
        elif len(points) == 4:
            n = np.array([np.array(points[0] - points[1]), np.array(points[0] - points[2]), np.array(points[0] - points[3])])
            norms = np.array(list(map(np.linalg.norm, n)))
            ind_norm_max = np.where(norms == max(norms))[0]
            n = np.delete(n, ind_norm_max, axis = 0)
            area = np.cross(n[0], n[1])
        self.mb.tag_set_data(self.area_tag, face, area)

    def calc_area(self, face):
        points = self.mtu.get_bridge_adjacencies(face, 2, 0)
        points = [self.mb.get_coords([vert]) for vert in points]
        if len(points) == 3:
            n1 = np.array(points[0] - points[1])
            n2 = np.array(points[0] - points[2])
            area = np.cross(n1, n2)/2.0

        #calculo da area para quadrilatero regular
        elif len(points) == 4:
            n = np.array([np.array(points[0] - points[1]), np.array(points[0] - points[2]), np.array(points[0] - points[3])])
            norms = np.array(list(map(np.linalg.norm, n)))
            ind_norm_max = np.where(norms == max(norms))[0]
            n = np.delete(n, ind_norm_max, axis = 0)
            area = np.cross(n[0], n[1])
        return(area)

    def get_kequiv_by_face_quad(self, conj_faces):
        """
        retorna os valores de k equivalente para colocar na matriz
        a partir da face

        input:
            face: face do elemento
        output:
            kequiv: k equivalente
            elems: elementos vizinhos pela face
            s: termo fonte da gravidade
        """
        ADJs=np.array([self.mb.get_adjacencies(face, 3) for face in self.all_faces])

        centroids=np.array([self.mtu.get_average_position([v]) for v in self.all_volumes])

        ADJsv=np.array([self.mb.get_adjacencies(face, 0) for face in self.all_faces])

        ks=self.mb.tag_get_data(self.perm_tag, self.all_volumes)
        #vol_to_pos=dict(zip(M1.all_volumes,range(len(M1.all_volumes))))
        vol_to_pos=dict(zip(self.all_volumes,range(len(self.all_volumes))))
        cont=0
        K_eq=[]
        for f in self.all_faces:
            adjs=ADJs[cont]
            adjsv=ADJsv[cont]
            cont+=1
            if len(adjs)==2:
                v1=adjs[0]
                v2=adjs[1]
                k1 = ks[vol_to_pos[v1]].reshape([3, 3])
                k2 = ks[vol_to_pos[v2]].reshape([3, 3])
                centroid1 = centroids[vol_to_pos[v1]]
                centroid2 = centroids[vol_to_pos[v2]]
                direction = centroid2 - centroid1
                norm=np.linalg.norm(direction)
                uni = np.absolute(direction/norm)
                k1 = np.dot(np.dot(k1,uni), uni)
                k2 = np.dot(np.dot(k2,uni), uni)
                k_harm = (2*k1*k2)/(k1+k2)

                vertex_cent=np.array([self.mb.get_coords([np.uint64(a)]) for a in adjsv])
                dx=max(vertex_cent[:,0])-min(vertex_cent[:,0])
                dy=max(vertex_cent[:,1])-min(vertex_cent[:,1])
                dz=max(vertex_cent[:,2])-min(vertex_cent[:,2])
                if dx<0.001:
                    dx=1
                if dy<0.001:
                    dy=1
                if dz<0.001:
                    dz=1
                area=dx*dy*dz
                #area = self.mb.tag_get_data(self.area_tag, face, flat=True)[0]
                #s_gr = self.gama*keq*(centroid2[2]-centroid1[2])
                keq = k_harm*area/(self.mi*norm)

                K_eq.append(keq)
            else:
                K_eq.append(0.0)
        self.mb.tag_set_data(self.k_eq_tag, self.all_faces, K_eq)

    def set_k_and_phi_structured_spe10(self):
        ks = np.load('spe10_perms_and_phi.npz')['perms']
        phi = np.load('spe10_perms_and_phi.npz')['phi']

        nx = 60
        ny = 220
        nz = 85
        perms = []
        phis = []

        k = 1.0  #para converter a unidade de permeabilidade
        centroids=np.array([self.mtu.get_average_position([v]) for v in self.all_volumes])
        cont=0
        for v in self.all_volumes:
            centroid = centroids[cont]
            cont+=1
            ijk = np.array([centroid[0]//20.0, centroid[1]//10.0, centroid[2]//2.0])
            e = int(ijk[0] + ijk[1]*nx + ijk[2]*nx*ny)
            # perm = ks[e]*k
            # fi = phi[e]
            perms.append(ks[e]*k)
            phis.append(phi[e])

        self.mb.tag_set_data(self.perm_tag, self.all_volumes, perms)
        self.mb.tag_set_data(self.phi_tag, self.all_volumes, phis)

    def get_boundary_nodes(self):
        all_faces = self.dirichlet_faces | self.neumann_faces
        boundary_nodes = set()
        for face in all_faces:
            nodes = self.mtu.get_bridge_adjacencies(face, 2, 0)
            boundary_nodes.update(nodes)
        return boundary_nodes

    def get_faces_boundary(self):
        """
        cria os meshsets
        all_faces_set: todas as faces do dominio
        all_faces_boundary_set: todas as faces no contorno
        """
        all_faces_boundary_set = self.mb.create_meshset()

        for face in self.all_faces:
            size = len(self.mb.get_adjacencies(face, 3))
            self.set_area(face)
            if size < 2:
                self.mb.add_entities(all_faces_boundary_set, [face])

        self.mb.tag_set_data(self.all_faces_boundary_tag, 0, all_faces_boundary_set)
        self.all_faces_boundary=self.mb.get_entities_by_handle(all_faces_boundary_set)
    def get_non_boundary_volumes(self, dirichlet_nodes, neumann_nodes):
        volumes = self.all_volumes
        non_boundary_volumes = []
        for volume in volumes:
            volume_nodes = set(self.mtu.get_bridge_adjacencies(volume, 0, 0))
            if (volume_nodes.intersection(dirichlet_nodes | neumann_nodes)) == set():
                non_boundary_volumes.append(volume)
        return non_boundary_volumes

    def set_media_property(self, property_name, physicals_values,
                           dim_target=3, set_nodes=False):

        self.set_information(property_name, physicals_values,
                             dim_target, set_connect=set_nodes)

    def set_boundary_condition(self, boundary_condition, physicals_values,
                               dim_target=3, set_nodes=False):

        self.set_information(boundary_condition, physicals_values,
                             dim_target, set_connect=set_nodes)

    def get_tetra_volume(self, tet_nodes):
        vect_1 = tet_nodes[1] - tet_nodes[0]
        vect_2 = tet_nodes[2] - tet_nodes[0]
        vect_3 = tet_nodes[3] - tet_nodes[0]
        vol_eval = abs(np.dot(np.cross(vect_1, vect_2), vect_3))/6.0
        return vol_eval

    def get_boundary_faces(self):
        all_boundary_faces = self.mb.create_meshset()
        for face in self.all_faces:
            elems = self.mtu.get_bridge_adjacencies(face, 2, 3)
            if len(elems) < 2:
                self.mb.add_entities(all_boundary_faces, [face])

        self.mb.tag_set_data(self.all_faces_boundary_tag, 0, all_boundary_faces)
        self.all_boundary_faces=self.mb.get_entities_by_handle(all_boundary_faces)



    @staticmethod
    def imprima(self, text = None):
        m1 = self.mb.create_meshset()
        self.mb.add_entities(m1, self.all_nodes)
        m2 = self.mb.create_meshset()
        self.mb.add_entities(m2, self.all_edges)
        m3 = self.mb.create_meshset()
        self.mb.add_entities(m3, self.all_faces)

        if text == None:
            text = "output"
        extension = ".vtk"
        text1 = text + "-nodes" + extension
        text2 = text + "-edges" + extension
        text3 = text + "-faces" + extension

        self.mb.write_file(text1,[m1])
        self.mb.write_file(text2,[m2])
        self.mb.write_file(text3,[m3])

        if len(self.all_volumes)>0:
            m4 = self.mb.create_meshset()
            self.mb.add_entities(m4, self.all_volumes)
            text4 = text + "-volumes" + extension
            self.mb.write_file(text4,[m4])

        print(text, "Arquivos gerados")
