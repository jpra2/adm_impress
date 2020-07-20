import pdb
import numpy as np

class TpfaPreprocess:
    '''
        Preprocessamento para tpfa
    '''

    def get_equivalent_permeability_faces_from_diagonal_permeability(self, volumes, faces, volumes_adj_internal_faces, volumes_adj_boundary_faces, unit_normal_vector, permeability, internal_faces, block_dimension):

        '''
            volumes: volumes
            faces: faces
            volumes_adj_faces: volumes vizinhos de face
            unit_normal_vector: vetor normal unitario das faces
            permeability: permeabilidade dos volumes
            internal_faces: faces_internas
            block_dimension: dimensões dos volumes

            output:
                eq_perm: permeabilidade equivalente nas faces

        '''

        boundary_faces = np.setdiff1d(faces, internal_faces)
        k_harm_faces = np.empty(len(internal_faces) + len(boundary_faces))
        u_normal = unit_normal_vector

        u_normal_internal_faces = np.absolute(unit_normal_vector[internal_faces])
        v0 = volumes_adj_internal_faces
        ni = len(internal_faces)

        ks0 = permeability[v0[:, 0]]
        ks1 = permeability[v0[:, 1]]

        ks0 = ks0 * u_normal_internal_faces.reshape([ni, 1, 3])
        ks1 = ks1 * u_normal_internal_faces.reshape([ni, 1, 3])
        ks0 = ks0.sum(axis=2).sum(axis=1)
        ks1 = ks1.sum(axis=2).sum(axis=1)
        hi = self.get_h_internal_faces(ni, block_dimension, v0, u_normal_internal_faces)
        k_harm_faces[internal_faces] = hi.sum(axis=1)/(hi[:, 0]/ks0 + hi[:, 1]/ks1)

        vols_viz_boundary_faces = volumes_adj_boundary_faces

        u_normal_b_faces = np.absolute(u_normal[boundary_faces])
        nb = len(boundary_faces)
        ks0 = permeability[vols_viz_boundary_faces]
        ks0 = ks0.reshape([nb, 3, 3]) * u_normal_b_faces.reshape([nb, 1, 3])
        ks0 = ks0.sum(axis=2).sum(axis=1)
        k_harm_faces[boundary_faces] = ks0

        return k_harm_faces

    def get_h_internal_faces(self, n_internal_faces, block_dimension, volumes_adj_internal_faces, u_normal_internal_faces):

        ni = n_internal_faces
        vols_viz_internal_faces = volumes_adj_internal_faces
        hs = block_dimension

        hi = np.zeros((ni, 2))
        hi[:, 0] = ((hs[vols_viz_internal_faces[:, 0]]*u_normal_internal_faces).sum(axis=1))/2
        hi[:, 1] = ((hs[vols_viz_internal_faces[:, 1]]*u_normal_internal_faces).sum(axis=1))/2

        return hi

    def get_areas_faces(self, faces, points_faces):

        areas = []

        for face in faces:
            points = points_faces[face]
            area = self.get_area_from_points(points, np.mean(points, axis=0))
            areas.append(area)

        return np.array(areas)

    def get_area_from_points(self, points, centroid):

        sorted_points = self.get_sorted_points(points, centroid)
        area = self.get_area_from_sorted_points(sorted_points, centroid)
        return area

    def get_sorted_points(self, points, centroid):

        lim = 1e-20

        dty = [('point', np.int), ('angle', np.float)]

        p0 = points[0]
        v0 = centroid - p0
        modv0 = np.linalg.norm(v0)

        reference = (np.cross(v0, centroid - points[1]))
        reference = reference/np.linalg.norm(reference)

        str_array = np.zeros(len(points), dtype=dty)
        str_array['point'] = np.arange(len(points))

        for i, point in enumerate(points[1:]):
            v1 = point - centroid
            modv1 = np.linalg.norm(v1)
            intern_product = np.dot(v1, v0)
            cross = np.cross(v0, v1)
            cross_norm = np.linalg.norm(cross)
            if cross_norm == 0:
                str_array['angle'][i+1] = np.pi
                continue
            ref2 = cross/cross_norm
            coss = intern_product/(modv1*modv0)
            angle = np.arccos(coss)

            if coss == 0:
                if np.allclose(ref2, reference):
                    pass
                else:
                    angle += np.pi
            elif np.allclose(ref2, reference):
                pass
            else:
                angle = 2*np.pi - angle

            str_array['angle'][i+1] = angle

        sorted_points = np.sort(str_array, order='angle')
        sorted_p = points[sorted_points['point']]
        return sorted_p

    def get_area_from_sorted_points(self, sorted_points, centroid):

        areas = []
        n = len(sorted_points)

        for i in range(n):
            tri_points = np.array([centroid, sorted_points[i], sorted_points[i-1]])
            area_tri = self.get_area_from_triangle(tri_points)
            areas.append(area_tri)

        return np.sum(areas)

    def get_area_from_triangle(self, points):

        if len(points) != 3:
            raise ValueError(f'triangle not defined: {len(points)} data')

        v0 = points[1] - points[0]
        v1 = points[2] - points[0]

        area_tri = np.linalg.norm(np.cross(v0, v1))/2

        return area_tri

    def get_u_normal_internal_faces(self, volumes_adj_internal_faces, abs_u_normal_internal_faces, centroid_volumes):

        delta_x = centroid_volumes[volumes_adj_internal_faces[:,1]] - centroid_volumes[volumes_adj_internal_faces[:,0]]
        normals = np.linalg.norm(delta_x, axis=1)
        normals = normals.reshape(len(normals), 1)
        u_normal_direction = delta_x/normals
        return u_normal_direction

    def get_k_volumes_internal_faces_direction(self, volumes_adj_internal_faces, permeability, abs_u_normal_internal_faces):
        perms_adjs = permeability[volumes_adj_internal_faces]
        perms_resp = np.empty((len(perms_adjs), 2))
        abs_u = abs_u_normal_internal_faces

        for i in range(perms_adjs.shape[0]):
            perms_resp[i,0] = np.dot(np.dot(perms_adjs[i, 0], abs_u[i]), abs_u[i])
            perms_resp[i,1] = np.dot(np.dot(perms_adjs[i, 1], abs_u[i]), abs_u[i])

        return perms_resp
