
import directories as direc
from impress.preprocessor import directories as direc_impress
from utils.utils_old import get_box
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

    def __init__(self, M):
        self.set_permeability_and_phi(M)
        self.set_area_hex_structured(M)
        self.set_k_harm_hex_structured(M)

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

        M.data.variables[direc.variables_impress['area']] = all_areas
        M.data.variables[direc.variables_impress['dist_cent']] = dist_cent

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

    def set_k_harm_hex_structured(self, M):
        '''
        considerando malha estruturada
        '''
        normals = np.absolute(M.faces.normal[:])
        vols_viz_faces = M.data.elements_lv0[direc_impress.entities_lv0_0[1]]
        internal_faces = M.data.elements_lv0[direc_impress.entities_lv0_0[0]]
        centroids_volumes = M.data.centroids[direc.entities_lv0[3]]
        ks = M.data.variables[direc.variables_impress['permeability']].copy()
        # ks = ks.reshape([len(ks), 3, 3])

        areas = M.data.variables[direc.variables_impress['area']]
        k_harm = np.zeros(len(areas))

        import pdb; pdb.set_trace()

        vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[2]]
        vols_viz_boundary_faces = M.data.elements_lv0[direc.entities_lv0_0[3]]

        import pdb; pdb.set_trace()
        shape_ks_0 = [vols_viz_internal_faces.shape[0], vols_viz_internal_faces.shape[1], 9]


        ks_vols_viz_internal_faces = np.zeros(shape_ks_0)
        ks_vols_viz_internal_faces[:,0] = ks[vols_viz_internal_faces[:,0]]
        ks_vols_viz_internal_faces[:,1] = ks[vols_viz_internal_faces[:,1]]

        for i, k01 in enumerate(ks_vols_viz_internal_faces):
            pass


        import pdb; pdb.set_trace()
