from .. import directories as direc
from . import directories_impress as direc_impress
from ..directories import data_loaded
from ..utils.utils_old import get_box
from math import pi
import numpy as np
from sympy.parsing.sympy_parser import parse_expr
from sympy import symbols, lambdify

# from .preprocess1 import set_saturation_regions


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

    # ijk0 = np.array([centroids[:, 0]//20.0, centroids[:, 1]//10.0, centroids[:, 2]//2.0])
    ijk0 = np.array([centroids[:, 0]//1.0, centroids[:, 1]//1.0, centroids[:, 2]//2.0])
    ee = ijk0[0] + ijk0[1]*nx + ijk0[2]*nx*ny
    ee = ee.astype(np.int32)

    M.data[M.data.variables_impress['permeability']] = ks[ee] #permeabilidade
    M.data[M.data.variables_impress['perm_z']] = ks[ee][:,-1]
    M.data[M.data.variables_impress['perm_x']] = ks[ee][:,0]

    M.data[M.data.variables_impress['poro']] = phi[ee]  #porosidade

class Preprocess0:
    '''
    objetivo: setar k_harmonico, pretransmissibilidade, transmissibilidade, area, dist_centroid, u_normal nas faces
    permeabilidade, volume, centroide, porosidade nos volumes
    '''

    def __init__(self, M, elements_lv0):

        self.elements_lv0 = elements_lv0

        for key, item in direc.variables_impress.items():
            M.data.variables_impress[key] = item

        self.run(M)
        M.data.export_to_npz()
        np.save("flying/permeability.npy",M.data[M.data.variables_impress['permeability']])
        M.data.update_variables_to_mesh()

    def set_area_hex_structured(self, M):
        """ Resolver Isso para malha diagonal amanha"""
        def get_area(ind, normals, nodes_faces, coord_nodes):

            if len(np.where(normals[:,ind] == 1)[0]) == 0:
                indice = np.where(normals[:,ind]!=0)[0][0]
            else: indice = np.where(normals[:,ind] == 1)[0][0]

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
                elif unis2[j][1] == 0 and unis[j][0]!=1 and unis[j][2]!=1: #caso inclinado em xz
                    hs[0] = modulos[j]

        areas = np.array(areas)
        all_areas = np.dot(normals**2, areas)
        dist_cent = np.dot(normals, hs)

        v = np.ones([3,1])
        dd = np.argwhere(np.round(normals@v, 1)!=1)
        if len(dd)>0:
            inclined_faces_normals = normals[dd.ravel()]
            xz_face = np.argwhere(inclined_faces_normals[:,1]==0).ravel()
            xz_face_angle = inclined_faces_normals[xz_face][0,2]/inclined_faces_normals[xz_face][0,0]
            theta = np.arctan(xz_face_angle)
        else: theta = pi/2

        volume = hs[0]*hs[1]*hs[2]*np.sin(theta)
        n_volumes = M.data.len_entities[direc.entities_lv0[3]]

        dd = np.zeros([n_volumes, 3])
        for i in range(3):
            dd[:, i] = np.repeat(hs[i], n_volumes)

        M.data[M.data.variables_impress['area']] = all_areas
        M.data[M.data.variables_impress['dist_cent']] = dist_cent
        M.data[M.data.variables_impress['volume']] = np.repeat(volume, n_volumes)
        M.data[M.data.variables_impress['NODES']] = coord_nodes
        M.data[M.data.variables_impress['hs']] = dd

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
            value = np.array(d0[direc.names_data_loaded_lv2[1]])
            value = np.array([float(v) for v in value.flatten()])
            value = value[np.newaxis,:]

            if tipo == direc.types_region_data_loaded[0]:
                # tamanho_variavel = M.data.info_data[M.data.variables_impress['permeability']][direc_impress.names_datas[0]]
                n_volumes = M.data.len_entities['volumes']
                # valor = valor.reshape([n, tamanho_variavel])
                M.data[M.data.variables_impress['permeability']] = np.repeat(value, n_volumes, axis=0)
            elif tipo == 'ring':
                ###########################################
                c0=d0['c0']
                r0=d0['r0']
                r1=d0['r1']
                x0=c0[0]
                y0=c0[1]
                z0=c0[2]

                x=centroids[:,0]
                y=centroids[:,1]
                z=centroids[:,2]
                cir_0=(x-x0)**2+(y-x0)**2+(z-z0)**2>r0**2
                cir_1=(x-x0)**2+(y-x0)**2+(z-z0)**2<r1**2
                indices=np.arange(len(centroids))[cir_0 & cir_1]
                n_volumes = len(indices)
                M.data[M.data.variables_impress['permeability']][indices] = np.repeat(value, n_volumes, axis=0)

                np.save("flying/permeability.npy",M.data[M.data.variables_impress['permeability']])

            elif tipo == 'cartesian_region':
                indices=np.repeat(False,len(centroids))
                x_inf=d0['x_inf']
                x_sup=d0['x_sup']
                y_inf=d0['y_inf']
                y_sup=d0['y_sup']
                z_inf=d0['z_inf']
                z_sup=d0['z_sup']
                xc=centroids[:,0]
                yc=centroids[:,1]
                zc=centroids[:,2]
                x, y, z = symbols('x y z')
                for i in range(len(x_inf)):
                    x_i=parse_expr(x_inf[i])
                    x_s=parse_expr(x_sup[i])
                    y_i=parse_expr(y_inf[i])
                    y_s=parse_expr(y_sup[i])
                    z_i=parse_expr(z_inf[i])
                    z_s=parse_expr(z_sup[i])
                    fx_i=lambdify([x,y,z],x_i)(xc,yc,zc)
                    fx_s=lambdify([x,y,z],x_s)(xc,yc,zc)
                    fy_i=lambdify([x,y,z],y_i)(xc,yc,zc)
                    fy_s=lambdify([x,y,z],y_s)(xc,yc,zc)
                    fz_i=lambdify([x,y,z],z_i)(xc,yc,zc)
                    fz_s=lambdify([x,y,z],z_s)(xc,yc,zc)
                    inds=(xc>fx_i) & (xc<fx_s) & (yc>fy_i) & (yc<fy_s) & (zc>fz_i) & (zc<fz_s)
                    indices=indices | inds
                indices=np.arange(len(centroids))[indices]
                n_volumes = len(indices)

                M.data[M.data.variables_impress['permeability']][indices] = np.repeat(value, n_volumes, axis=0)
                np.save("flying/permeability.npy",M.data[M.data.variables_impress['permeability']])
            elif tipo == 'polar_region':
                r, t = symbols('r t')
                indices=np.repeat(False,len(centroids))
                for i in range(len(d0['center'])):
                    center=d0['center'][i]
                    r_inf=parse_expr(d0['r_inf'][i])
                    r_sup=parse_expr(d0['r_sup'][i])
                    t_inf=parse_expr(d0['t_inf'][i])
                    t_sup=parse_expr(d0['t_sup'][i])
                    xc=centroids[:,0]-center[0]
                    yc=centroids[:,1]-center[1]
                    zc=centroids[:,2]-center[2]
                    rs=np.sqrt(xc*xc+yc*yc)
                    ts=np.arccos(xc/rs)
                    ts[yc<0]=2*3.141592653589793-ts[yc<0]
                    r_i=lambdify([r, t], r_inf)(rs,ts)
                    r_s=lambdify([r, t], r_sup)(rs,ts)
                    t_i=lambdify([r, t], t_inf)(rs,ts)
                    t_s=lambdify([r, t], t_sup)(rs,ts)
                    inds = (rs>r_i) & (rs<r_s) & (ts>t_i) & (ts<t_s)
                    indices=indices | inds
                indices=np.arange(len(centroids))[indices]
                n_volumes = len(indices)
                M.data[M.data.variables_impress['permeability']][indices] = np.repeat(value, n_volumes, axis=0)
                np.save("flying/permeability.npy",M.data[M.data.variables_impress['permeability']])

            elif tipo == direc.types_region_data_loaded[1]:
                p0 = d0[direc.names_data_loaded_lv2[2]]
                p1 = d0[direc.names_data_loaded_lv2[3]]
                points = np.array([np.array(p0), np.array(p1)])
                indices = get_box(centroids, points)
                n_volumes = len(indices)
                M.data[M.data.variables_impress['permeability']][indices] = np.repeat(value, n_volumes, axis=0)
                np.save("flying/permeability.npy",M.data[M.data.variables_impress['permeability']])


        #######################
        ## deletar
        permx = M.data[M.data.variables_impress['permeability']][:, 0]
        M.data['verif_po'][:] = permx
        #######################

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
                M.data[M.data.variables_impress['poro']] = np.repeat(value, n_volumes)
                M.data['poro'] = np.repeat(value, n_volumes)

            elif tipo == direc.types_region_data_loaded[1]:
                p0 = d0[direc.names_data_loaded_lv2[2]]
                p1 = d0[direc.names_data_loaded_lv2[3]]
                points = np.array([np.array(p0), np.array(p1)])
                indices = get_box(centroids, points)
                n_volumes = len(indices)
                M.data[M.data.variables_impress['poro']][indices] = np.repeat(value, n_volumes)

    def set_k_harm_hex_structured(self, M):
        '''
        considerando malha estruturada uniforme
        '''

        u_normal = M.data[M.data.variables_impress['u_normal']]
        vols_viz_faces = self.elements_lv0['neig_faces']
        internal_faces = self.elements_lv0['internal_faces']
        #internal_faces = internal_faces_or[M.faces.center[internal_faces_or][:,2]!=0] #just for now
        boundary_faces = self.elements_lv0['boundary_faces']
        centroids_volumes = M.data['centroid_volumes']
        ks = M.data[M.data.variables_impress['permeability']].copy()
        # dist_cent = M.data.variables[M.data.variables_impress['dist_cent']]
        # areas = M.data.variables[M.data.variables_impress['area']]
        hs = M.data[M.data.variables_impress['hs']]
        vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        vols_viz_boundary_faces = self.elements_lv0['neig_boundary_faces']

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

        # keq_internal_faces = ((ks0/hi[:,0])*(ks1/hi[:,1]))/((ks0/hi[:,0]) + (ks1/hi[:,1]))

        #Make unitary dimentions#Atemçao
        # hi=np.ones_like(hi)

        k_harm_faces[internal_faces] = hi.sum(axis=1)/(hi[:, 0]/ks0 + hi[:, 1]/ks1)

        u_normal_b_faces = u_normal[boundary_faces]
        nb = len(boundary_faces)
        ks0 = ks[vols_viz_boundary_faces]
        ks0 = ks0.reshape([nb, 3, 3]) * u_normal_b_faces.reshape([nb, 1, 3])
        ks0 = ks0.sum(axis=2).sum(axis=1)
        k_harm_faces[boundary_faces] = ks0
        # k_harm_faces[boundary_faces] = 0.0

        M.data[M.data.variables_impress['k_harm']] = k_harm_faces

    def set_transmissibility_monophasic(self, M):
        # mi_monofasic = 1.0
        mi_monofasic = direc.data_loaded['monophasic_data']['mi']
        pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
        transmissibility = M.data[M.data.variables_impress['transmissibility']].copy()
        transmissibility = pretransmissibility_faces/mi_monofasic

        M.data[M.data.variables_impress['transmissibility']] = transmissibility

    def update_centroids_and_unormal(self, M):

        M.data['centroid_volumes'] = M.volumes.center(M.volumes.all)
        M.data['centroid_faces'] = M.faces.center(M.faces.all)
        M.data['centroid_edges'] = M.edges.center(M.edges.all)
        M.data['centroid_nodes'] = M.nodes.center(M.nodes.all)
        M.data['u_normal'] = np.absolute(M.faces.normal[:])
        M.data['NODES'] = M.data['centroid_nodes'].copy()

    def set_pretransmissibility(self, M):
        areas = M.data['area'].copy()
        areas=np.ones_like(areas)
        k_harm_faces = M.data['k_harm']
        dist_cent = M.data['dist_cent'].copy()
        # dist_cent = np.ones_like(dist_cent)
        pretransmissibility_faces = (areas*k_harm_faces)/dist_cent
        M.data[M.data.variables_impress['pretransmissibility']] = pretransmissibility_faces
        # M.data.update_variables_to_mesh([M.data.variables_impress['pretransmissibility']])

    def initial_gama(self, M):
        gama_mono = data_loaded['monophasic_data']['gama']
        M.data['gama'] = np.repeat(gama_mono, len(M.data['gama']))

    def adj_matrix_volumes_volumes(self, M):

        vols = self.elements_lv0['volumes']
        n = len(vols)
        vizs = M.volumes.bridge_adjacencies(vols, 2, 3)
        adj_matrix = np.full((n, n), False, dtype=bool)
        lines = []
        cols = []
        for v, viz in zip(vols, vizs):
            n2 = len(viz)
            lines.append(np.repeat(v, n2).astype(int))
            cols.append(viz)

        lines = np.concatenate(lines)
        cols = np.concatenate(cols)
        data = np.full(len(cols), True, dtype=bool)
        adj_matrix[lines, cols] = data
        self.elements_lv0['adj_matrix_volumes_volumes'] = adj_matrix

    def run(self, M):
        self.update_centroids_and_unormal(M)
        self.set_permeability_and_phi(M)
        self.set_area_hex_structured(M)
        self.set_k_harm_hex_structured(M)
        self.set_pretransmissibility(M)
        self.set_transmissibility_monophasic(M)
        self.initial_gama(M)
        self.adj_matrix_volumes_volumes(M)
