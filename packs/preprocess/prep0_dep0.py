from .. import directories as direc
from . import directories_impress as direc_impress
from ..utils.utils_old import get_box
import numpy as np
from .preprocess1 import set_saturation_regions


def set_permeability_and_phi_spe10(M):
    data_spe10 = np.load(direc.data_loaded['file_name_permeability'])
    ks = data_spe10['perms']
    phi = data_spe10['phi']
    phi = phi.flatten()

    nx = 60
    ny = 220
    nz = 85

    # centroids = M.data.centroids[direc.entities_lv0[3]]
    centroids = M.volumes.center(M.volumes.all)

    ijk0 = np.array([centroids[:, 0]//20.0, centroids[:, 1]//10.0, centroids[:, 2]//2.0])
    ee = ijk0[0] + ijk0[1]*nx + ijk0[2]*nx*ny
    ee = ee.astype(np.int32)

    M.data.variables[M.data.variables_impress['permeability']] = ks[ee] #permeabilidade
    M.data.variables[M.data.variables_impress['poro']] = phi[ee]  #porosidade

class Preprocess0:
    '''
    objetivo: setar k_harmonico, pretransmissibilidade, transmissibilidade, area, dist_centroid, u_normal nas faces
    permeabilidade, volume, centroide, porosidade nos volumes
    '''

    def __init__(self, M):
        self.mesh = M

        for key, item in direc.variables_impress.items():
            M.data.variables_impress[key] = item

        self.update_centroids_and_unormal(M)
        self.set_permeability_and_phi(M)
        self.set_area_hex_structured(M)
        self.set_k_harm_hex_structured(M)
        self.set_pretransmissibility(M)
        self.set_transmissibility_monophasic(M)
        # if direc.data_loaded['biphasic'] == True:
        #     set_saturation_regions(M)

        M.state = 0
        M.data.update_variables_to_mesh()
        M.data.save_info_data()
        M.data.export_variables_to_npz(direc.names_outfiles_variables_steps[0])
        M.core.print(file=direc.output_file+str(M.state))
        np.save(direc.state_path, np.array([M.state]))
        np.save(direc.path_local_last_file_name, np.array([direc.names_outfiles_steps[0]]))

    def set_area_hex_structured(self, M):

        def get_area(ind, normals, nodes_faces, coord_nodes):
            indice = np.where(normals[:,ind] == 1)[0][0]
            nos = nodes_faces[indice]
            normas = []
            vetores = []
            maior = 0
            ind_maior = 0
            for i in range(3):
                v = coord_nodes[nos[0]] - coord_nodes[nos[i+1]]
                vetores.append(v)
                normas.append(np.linalg.norm(v))
                if normas[i] > maior:
                    maior = normas[i]
                    ind_maior = i

            del normas[ind_maior]
            del vetores[ind_maior]

            area = normas[0] * normas[1]
            vetores = np.absolute(np.array(vetores))
            normas = np.array(normas)

            return area, vetores, normas

        unis = np.array([np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])])
        faces = M.faces.all
        n_faces = len(faces)
        normals = np.absolute(M.faces.normal[:])
        nodes_faces = M.faces.bridge_adjacencies(faces, 2, 0)
        coord_nodes = M.nodes.center(M.nodes.all)
        areas = []
        hs = np.zeros(3)
        for i in range(3):
            area, vetores, modulos = get_area(i, normals, nodes_faces, coord_nodes)
            unis2 = np.zeros([2,3])
            unis2[0] = vetores[0]/modulos[0]
            unis2[1] = vetores[1]/modulos[1]
            areas.append(area)
            for j in range(2):
                if np.allclose(unis2[j], unis[0]):
                    hs[0] = modulos[j]
                elif np.allclose(unis2[j], unis[1]):
                    hs[1] = modulos[j]
                elif np.allclose(unis2[j], unis[2]):
                    hs[2] = modulos[j]

        areas = np.array(areas)
        all_areas = np.dot(normals, areas)
        dist_cent = np.dot(normals, hs)
        volume = hs[0]*hs[1]*hs[2]
        n_volumes = M.data.len_entities[direc.entities_lv0[3]]

        dd = np.zeros([n_volumes, 3])
        for i in range(3):
            dd[:, i] = np.repeat(hs[i], n_volumes)

        M.data.variables[M.data.variables_impress['area']] = all_areas
        M.data.variables[M.data.variables_impress['dist_cent']] = dist_cent
        M.data.variables[M.data.variables_impress['volume']] = np.repeat(volume, n_volumes)
        M.data.variables[M.data.variables_impress['NODES']] = coord_nodes
        M.data.variables[M.data.variables_impress['hs']] = dd

    def set_permeability_and_phi(self, M):

        data_loaded = direc.data_loaded

        read = data_loaded[direc.names_data_loaded_lv0[0]]
        set_perm = data_loaded[direc.names_data_loaded_lv0[5]]
        set_phi = data_loaded[direc.names_data_loaded_lv0[6]]

        if read:
            set_permeability_and_phi_spe10(M)
        if set_perm:
            self.set_permeability_regions(M)
        if set_phi:
            self.set_phi_regions(M)

    def set_permeability_regions(self, M):

        centroids = M.data['centroid_volumes']

        for reg in direc.data_loaded[direc.names_data_loaded_lv0[2]]:
            d0 = direc.data_loaded[direc.names_data_loaded_lv0[2]][reg]
            tipo = d0[direc.names_data_loaded_lv2[0]]
            value = np.array([np.array(d0[direc.names_data_loaded_lv2[1]])])

            if tipo == direc.types_region_data_loaded[0]:
                # tamanho_variavel = M.data.info_data[M.data.variables_impress['permeability']][direc_impress.names_datas[0]]
                n_volumes = M.data.len_entities['volumes']
                # valor = valor.reshape([n, tamanho_variavel])
                M.data.variables[M.data.variables_impress['permeability']] = np.repeat(value, n_volumes, axis=0)

            elif tipo == direc.types_region_data_loaded[1]:
                p0 = d0[direc.names_data_loaded_lv2[2]]
                p1 = d0[direc.names_data_loaded_lv2[3]]
                points = np.array([np.array(p0), np.array(p1)])
                indices = get_box(centroids, points)
                n_volumes = len(indices)
                M.data.variables[M.data.variables_impress['permeability']][indices] = np.repeat(value, n_volumes, axis=0)

    def set_phi_regions(self, M):
        # TODO: atualizar essa funcao
        centroids = M.data['centroid_volumes']

        for reg in direc.data_loaded[direc.names_data_loaded_lv0[7]]:
            d0 = direc.data_loaded[direc.names_data_loaded_lv0[7]][reg]
            tipo = d0[direc.names_data_loaded_lv2[0]]
            value = d0[direc.names_data_loaded_lv2[1]]
            assert isinstance(value, float)

            if tipo == direc.types_region_data_loaded[0]:
                n_volumes = M.data.len_entities['volumes']
                M.data.variables[M.data.variables_impress['poro']] = np.repeat(value, n_volumes)

            elif tipo == direc.types_region_data_loaded[1]:
                p0 = d0[direc.names_data_loaded_lv2[2]]
                p1 = d0[direc.names_data_loaded_lv2[3]]
                points = np.array([np.array(p0), np.array(p1)])
                indices = get_box(centroids, points)
                n_volumes = len(indices)
                M.data.variables[M.data.variables_impress['poro']][indices] = np.repeat(value, n_volumes)

    def set_k_harm_hex_structured(self, M):
        '''
        considerando malha estruturada uniforme
        '''

        u_normal = M.data.variables[M.data.variables_impress['u_normal']]
        vols_viz_faces = M.data.elements_lv0[direc_impress.entities_lv0_0[1]]
        internal_faces = M.data.elements_lv0[direc_impress.entities_lv0_0[0]]
        boundary_faces = M.data.elements_lv0[direc_impress.entities_lv0_0[4]]
        centroids_volumes = M.data['centroid_volumes']
        ks = M.data.variables[M.data.variables_impress['permeability']].copy()
        # dist_cent = M.data.variables[M.data.variables_impress['dist_cent']]
        # areas = M.data.variables[M.data.variables_impress['area']]
        hs = M.data.variables[M.data.variables_impress['hs']]
        vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[2]]
        vols_viz_boundary_faces = M.data.elements_lv0[direc.entities_lv0_0[3]]

        k_harm_faces = np.zeros(len(internal_faces) + len(boundary_faces))
        pretransmissibility_faces = k_harm_faces.copy()

        u_normal_internal_faces = u_normal[internal_faces]
        ni = len(internal_faces)
        ks0 = ks[vols_viz_internal_faces[:, 0]]
        ks1 = ks[vols_viz_internal_faces[:, 1]]
        ks0 = ks0.reshape([ni, 3, 3]) * u_normal_internal_faces.reshape([ni, 1, 3])
        ks1 = ks1.reshape([ni, 3, 3]) * u_normal_internal_faces.reshape([ni, 1, 3])
        ks0 = ks0.sum(axis=2).sum(axis=1)
        ks1 = ks1.sum(axis=2).sum(axis=1)

        hi = np.zeros((ni, 2))
        hi[:, 0] = ((hs[vols_viz_internal_faces[:, 0]]*u_normal_internal_faces).sum(axis=1))/2
        hi[:, 1] = ((hs[vols_viz_internal_faces[:, 1]]*u_normal_internal_faces).sum(axis=1))/2
        k_harm_faces[internal_faces] = hi.sum(axis=1)/(hi[:, 0]/ks0 + hi[:, 1]/ks1)

        u_normal_b_faces = u_normal[boundary_faces]
        nb = len(boundary_faces)
        ks0 = ks[vols_viz_boundary_faces.flatten()]
        ks0 = ks0.reshape([nb, 3, 3]) * u_normal_b_faces.reshape([nb, 1, 3])
        ks0 = ks0.sum(axis=2).sum(axis=1)
        k_harm_faces[boundary_faces] = ks0

        # pretransmissibility_faces = (areas*k_harm_faces)/dist_cent

        M.data.variables[M.data.variables_impress['k_harm']] = k_harm_faces
        # M.data.variables[M.data.variables_impress['pretransmissibility']] = pretransmissibility_faces

    def set_transmissibility_monophasic(self, M):
        # mi_monofasic = 1.0
        mi_monofasic = direc.data_loaded['monophasic_data']['mi']
        pretransmissibility_faces = M.data.variables[M.data.variables_impress['pretransmissibility']]
        transmissibility = M.data.variables[M.data.variables_impress['transmissibility']]
        transmissibility = pretransmissibility_faces/mi_monofasic
        M.data.variables[M.data.variables_impress['transmissibility']] = transmissibility

    def update_centroids_and_unormal(self, M):

        M.data['centroid_volumes'] = M.volumes.center(M.volumes.all)
        M.data['centroid_faces'] = M.faces.center(M.faces.all)
        M.data['centroid_edges'] = M.edges.center(M.edges.all)
        M.data['centroid_nodes'] = M.nodes.center(M.nodes.all)
        M.data['u_normal'] = np.absolute(M.faces.normal[:])
        M.data['NODES'] = M.data['centroid_nodes'].copy()

    def set_pretransmissibility(self, M):
        areas = M.data['area']
        k_harm_faces = M.data['k_harm']
        dist_cent = M.data['dist_cent']
        pretransmissibility_faces = (areas*k_harm_faces)/dist_cent
        M.data[M.data.variables_impress['pretransmissibility']] = pretransmissibility_faces
        M.data.update_variables_to_mesh([M.data.variables_impress['pretransmissibility']])

    def conectivity_volumes(self, M):


        import networkx as nx

        vols_viz_faces = M.data.elements_lv0[direc.entities_lv0_0[3]]
        internal_faces = M.elements_lv0[direc.entities_lv0_0[0]]

        import pdb; pdb.set_trace()

        g = nx.Graph()
        for vols in vols_viz_faces:
            g.add_node(vols[0])
            g.add_node(vols[1])
            g.add_edge(vols[0], vols[1])

        M.graph_volumes = g

        return 0
