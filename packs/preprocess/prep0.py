from packs import directories as direc
from impress.preprocessor import directories as direc_impress
from packs.utils import get_box
import numpy as np


def set_permeability_and_phi_spe10(M):
    ks = np.load(direc.data_loaded['file_name_permeability'])['perms']
    phi = np.load(direc.data_loaded['file_name_permeability'])['phi']
    phi = phi.flatten()

    nx = 60
    ny = 220
    nz = 85
    perms = []
    phis = []

    k = 1.0  #para converter a unidade de permeabilidade
    centroids=M.data.centroids[direc.entities_lv0[3]]

    for centroid in centroids:
        ijk = np.array([centroid[0]//20.0, centroid[1]//10.0, centroid[2]//2.0])
        e = int(ijk[0] + ijk[1]*nx + ijk[2]*nx*ny)
        permeabilidade = ks[e]
        permeabilidade *= k
        perms.append(permeabilidade)
        phis.append(phi[e])

    perms = np.array(perms)
    phis = np.array(phis)

    M.data.variables[direc.variables_impress['permeability']] = perms #permeabilidade
    M.data.variables[direc.variables_impress['poro']] = phis  #porosidade

class Preprocess0:
    '''
    objetivo: setar k_harmonico, pretransmissibilidade, transmissibilidade, area, dist_centroid, u_normal nas faces
    permeabilidade, volume, centroide, porosidade nos volumes
    '''

    def __init__(self, M):

        for key, item in direc.variables_impress.items():
            M.data.variables_impress[key] = item

        self.set_permeability_and_phi(M)
        self.set_area_hex_structured(M)
        self.set_k_harm_and_pretransmissibility_hex_structured(M)
        self.set_transmissibility_monofasic(M)

        M.state = 0
        M.data.update_variables_to_mesh()
        M.data.save_info_data()
        M.data.export_variables_to_npz(direc.names_outfiles_variables_steps[0])
        M.core.print(text=direc.output_file+str(M.state))
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

        M.data.variables[direc.variables_impress['area']] = all_areas
        M.data.variables[direc.variables_impress['dist_cent']] = dist_cent
        M.data.variables[M.data.variables_impress['volume']] = np.repeat(volume, n_volumes)
        M.data.variables[M.data.variables_impress['NODES']] = coord_nodes
        M.data.variables[M.data.variables_impress['hs']] = dd

    def set_permeability_and_phi(self, M):
        data_loaded = direc.data_loaded

        read = data_loaded[direc.names_data_loaded_lv0[0]]
        if data_loaded[direc.names_data_loaded_lv0[2]]:
            if read:
                set_permeability_and_phi_spe10(M)
            else:
                self.set_permeability_regions(M)
                self.set_phi_regions(M)

    def set_permeability_regions(self, M):

        centroids = M.data.centroids[direc.entities_lv0[3]]
        n = len(centroids)

        for reg in direc.data_loaded[direc.names_data_loaded_lv0[2]]:
            d0 = direc.data_loaded[direc.names_data_loaded_lv0[2]][reg]
            tipo = d0[direc.names_data_loaded_lv2[0]]
            value = np.array(d0[direc.names_data_loaded_lv2[1]])

            if tipo == direc.types_region_data_loaded[0]:
                tamanho_variavel = M.data.info_data[direc.variables_impress['permeability']][direc_impress.names_datas[0]]
                # valor = valor.reshape([n, tamanho_variavel])
                for i in range(n):
                    M.data.variables[direc.variables_impress['permeability']][i] = value

            elif tipo == direc.types_region_data_loaded[1]:
                p0 = d0[direc.names_data_loaded_lv2[2]]
                p1 = d0[direc.names_data_loaded_lv2[3]]
                points = np.array([np.array(p0), np.array(p1)])
                indices = get_box(centroids, points)
                for i in indices:
                    M.data.variables[direc.variables_impress['permeability']][i] = value

    def set_phi_regions(self, M):
        # TODO: atualizar essa funcao
        centroids = M.data.centroids[direc.entities_lv0[3]]
        n = len(centroids)
        values = np.repeat(0.3, n)
        M.data.variables[direc.variables_impress['poro']] = values

    def set_k_harm_and_pretransmissibility_hex_structured(self, M):
        '''
        considerando malha estruturada uniforme
        '''
        normals = np.absolute(M.faces.normal[:])
        M.data.variables[direc.variables_impress['u_normal']] = normals
        vols_viz_faces = M.data.elements_lv0[direc_impress.entities_lv0_0[1]]
        internal_faces = M.data.elements_lv0[direc_impress.entities_lv0_0[0]]
        boundary_faces = M.data.elements_lv0[direc_impress.entities_lv0_0[4]]
        centroids_volumes = M.data.centroids[direc.entities_lv0[3]]
        ks = M.data.variables[direc.variables_impress['permeability']].copy()
        dist_cent = M.data.variables[direc.variables_impress['dist_cent']]
        areas = M.data.variables[direc.variables_impress['area']]
        k_harm_faces = M.data.variables[direc.variables_impress['k_harm']]
        hs = M.data.variables[M.data.variables_impress['hs']]
        pretransmissibility_faces = M.data.variables[M.data.variables_impress['pretransmissibility']]

        vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[2]]
        vols_viz_boundary_faces = M.data.elements_lv0[direc.entities_lv0_0[3]]
        shape_ks_0 = [vols_viz_internal_faces.shape[0], vols_viz_internal_faces.shape[1], 9]
        vols_viz_faces = M.data.elements_lv0[direc.entities_lv0_0[1]]

        ks_vols_viz_internal_faces = np.zeros(shape_ks_0)
        ks_vols_viz_internal_faces[:,0] = ks[vols_viz_internal_faces[:,0]]
        ks_vols_viz_internal_faces[:,1] = ks[vols_viz_internal_faces[:,1]]

        for i, f in enumerate(internal_faces):
            k0 = ks_vols_viz_internal_faces[i,0].reshape([3,3])
            k1 = ks_vols_viz_internal_faces[i,1].reshape([3,3])
            h0 = hs[vols_viz_internal_faces[i, 0]]
            h1 = hs[vols_viz_internal_faces[i, 1]]
            u_normal = normals[f]
            h0 = np.dot(h0, u_normal)
            h1 = np.dot(h1, u_normal)
            area = areas[f]
            dist = dist_cent[f]
            normal = normals[f]
            k00 = np.dot(np.dot(k0, normal), normal)
            k11 = np.dot(np.dot(k1, normal), normal)
            # k_harm = (2*k00*k11)/(k00+k11)
            k_harm = (h0+h1)/(h0/k00 + h1/k11)
            k_harm_faces[f] = k_harm
            pretransmissibility_faces[f] = (area*k_harm)/(dist)

        shape_ks_0 = [vols_viz_boundary_faces.shape[0], 9]
        ks_vols_viz_boundary_faces = np.zeros(shape_ks_0)
        ks_vols_viz_boundary_faces[:] = ks[vols_viz_boundary_faces[:]]

        for i, f in enumerate(boundary_faces):
            k0 = ks_vols_viz_boundary_faces[i].reshape([3,3])
            area = areas[f]
            dist = dist_cent[f]
            normal = normals[f]
            k00 = np.dot(np.dot(k0, normal), normal)
            k_harm = k00
            k_harm_faces[f] = k_harm
            pretransmissibility_faces[f] = (area*k_harm)/(dist)

    def set_transmissibility_monofasic(self, M):
        # mi_monofasic = 1.0
        mi_monofasic = direc.data_loaded['monophasic_data']['mi']
        pretransmissibility_faces = M.data.variables[M.data.variables_impress['pretransmissibility']]
        transmissibility = M.data.variables[M.data.variables_impress['transmissibility']]
        transmissibility = pretransmissibility_faces/mi_monofasic
        M.data.variables[M.data.variables_impress['transmissibility']] = transmissibility

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
