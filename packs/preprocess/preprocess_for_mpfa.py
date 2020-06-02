import pdb
from packs.data_class.data_manager import DataManager
from packs.directories import only_mesh_name
import numpy as np

class PreprocessForMpfa(DataManager):

    def __init__(self, data_name='PreprocessForMpfa', load=False):
        data_name = data_name + '_' + only_mesh_name + '.npz'
        super().__init__(data_name=data_name, load=load)
        self.type_int = np.dtype(int)
        self.dt_for_points_subfaces = [('centroid_of_face', self.type_int), ('midpoint_of_edge_1', self.type_int), ('centroid_of_node', self.type_int), ('midpoint_of_edge_2', self.type_int)]
        self.name_cent_volume = 'centroid_of_volume'
        self.name_cent_face = 'centroid_of_face'
        self.name_cent_node = 'centroid_of_node'

    def create_midpoints_edges(self, edges, nodes, centroid_nodes, nodes_edges):
        '''
            input:
                edges: id das edges
                nodes: id dos nodes
                centroid_nodes: ponto dos nodes
                nodes_edges: nodes das edges
            output:
                midpoints_edges
        '''
        midpoints_edges = np.zeros((len(edges), 3))
        cents_nodes_edges = centroid_nodes[nodes_edges]
        midpoints_edges = np.mean(cents_nodes_edges, axis=1)
        return midpoints_edges

    def create_subfaces(self, faces, edges, edges_faces, nodes_faces, centroid_faces, centroid_nodes, midpoints_edges, set_ids=False):
        '''
            define os pontos que pertencem as subfaces
            input:
                faces: faces_ids
                edges: edges ids
                edges_faces: edges das faces
                nodes_faces: nodes das faces
                centroid_faces: centroid das faces
                centroid_nodes: centroid dos nodes
            output:
                subfaces

                todas_normais: normais das subfaces
                todos_pontos: pontos das subfaces
                local_ids_subfaces: local id das subfaces
        '''

        todas_normais = []
        todos_pontos = []
        local_ids_subfaces = []

        for face in faces:
            cent_face = centroid_faces[face]
            edges = edges_faces[face]
            mid_edges = midpoints_edges[edges]
            nodes = nodes_faces[face]
            cent_nodes = centroid_nodes[nodes]

            normais, points_areas, local_id_subfaces = self.func1(face, cent_face, edges, mid_edges, nodes, cent_nodes, set_ids)
            todos_pontos.append(points_areas)
            todas_normais.append(normais)
            local_ids_subfaces.append(local_id_subfaces)

        todas_normais = np.array(todas_normais)
        todos_pontos = np.array(todos_pontos)
        local_ids_subfaces = np.array(local_ids_subfaces)

        return todas_normais, todos_pontos, local_ids_subfaces

    def func1(self, face, cent_face, edges, midpoint_edges, nodes, cent_nodes, set_ids):
        '''
        input:
            cent_face: centroid da face
            edges: edges das faces
            midpoint_edges: midpoint das edges
            nodes: nodes das faces
            cent_nodes: centroid dos nodes
        output:
            normais: normais das subfaces
            points_areas: vertices das subfaces
            local_id_subfaces: local id das subfaces
        '''

        mid_edges = midpoint_edges
        vectors_midpoints = mid_edges - cent_face
        vectors_nodes = cent_nodes - cent_face
        norms_vectors_midpoins = np.linalg.norm(vectors_midpoints, axis=1)
        norms_vectors_nodes = np.linalg.norm(vectors_nodes, axis=1)
        unit_sentido_test = np.cross(vectors_nodes[0], vectors_midpoints[0])
        unit_sentido_test = unit_sentido_test/np.linalg.norm(unit_sentido_test)

        normais = []
        points_areas = []
        all_points2 = []

        for i, vec1 in enumerate(vectors_nodes):
            norm1 = norms_vectors_nodes[i]
            angles = []
            all_norms = []
            for j, vec2 in enumerate(vectors_midpoints):
                norm2 = norms_vectors_midpoins[j]
                angle, normal = self.func2(vec1, vec2, norm1, norm2, unit_sentido_test)
                angles.append(angle)
                all_norms.append(normal)

            angles = np.array(angles)
            ids = np.arange(len(angles))

            imin = ids[angles == angles.min()][0]
            imax = ids[angles == angles.max()][0]

            midp1 = mid_edges[imin]
            midp2 = mid_edges[imax]

            points = np.array([cent_face, midp1, cent_nodes[i], midp2])
            if set_ids == True:
                points2 = np.array([face, edges[imin], nodes[i], edges[imax]])
                all_points2.append(points2)

            normais.append((all_norms[imax] + all_norms[imin])/2)
            points_areas.append(points)

        normais = np.array(normais)
        points_areas = np.array(points_areas)
        local_id_subfaces = np.arange(len(normais))
        all_points2 = np.array(all_points2)

        if set_ids == True:
            dt = self.dt_for_points_subfaces
            structured_array = np.zeros(len(all_points2), dtype=dt)
            structured_array['centroid_of_face'] = all_points2[:,0]
            structured_array['midpoint_of_edge_1'] = all_points2[:,1]
            structured_array['centroid_of_node'] = all_points2[:,2]
            structured_array['midpoint_of_edge_2'] = all_points2[:,3]
            points_areas = structured_array

        return normais, points_areas, local_id_subfaces

    def func2(self, vec1, vec2, norm1, norm2, unit_sentido_test):
        '''
        output:
            angle: angle between vec1, vec2
            normal: normal between vec1, vec2
        '''

        angle = (180/np.pi)*np.arccos(np.dot(vec1, vec2)/(norm1 + norm2))
        normal = np.cross(vec1, vec2)
        unit2 = normal/np.linalg.norm(normal)
        if np.allclose(unit_sentido_test, unit2):
            pass
        else:
            angle = 360 - angle
            normal *= -1

        return angle, normal

    def get_points_of_subfaces_from_structured_array(self, structured_array, centroid_faces, centroid_nodes, midpoints_edges):

        arr = structured_array
        cent_faces = centroid_faces[arr['centroid_of_face']]
        mipd1 = midpoints_edges[arr['midpoint_of_edge_1']]
        cent_nodes = centroid_nodes[arr['centroid_of_node']]
        midp2 = midpoints_edges[arr['midpoint_of_edge_2']]
        poo = np.hstack([cent_faces, mipd1, cent_nodes, midp2])
        poo = poo.reshape((poo.shape[0], 4, 3))
        return poo

    def get_points_of_subfaces_from_complete_structured_array(self, comp_structured_array, centroid_faces, centroid_nodes, midpoints_edges):

        points = []

        for arr in comp_structured_array:
            points.append(
                self.get_points_of_subfaces_from_structured_array(
                    arr,
                    centroid_faces,
                    centroid_nodes,
                    midpoints_edges
                )
            )

        points = np.array(points)
        return points

    def create_all_subvolumes_mpfao(self, volumes, faces, faces_volumes, faces_nodes, nodes_volumes):
        '''
            input:
                volumes: volumes
                faces: faces
                centroid_volumes: centroid volumes
                centroid_faces: centroid faces
                faces_volumes: faces dos volumes
                faces_nodes: faces dos nodes
                nodes_volumes: nodes dos volumes
            output:
                subvolumes
        '''

        ## tamanho de todos subvolumes
        gids_subvolumes = []

        # ## gids faces of subvolumes
        # gids_faces_of_subvolumes = []

        # ## tamanho das faces dos subvolumes
        # primeiro_volume = []

        # tamanho de todos subvolumes
        points_of_subvolumes = []

        # # tamanho de todas as faces dos subvolumes
        # points_of_faces_of_subvolumes = []

        # ## tamanho de todas as faces dos subvolumes
        # normals_of_faces_of_subvolumes = []

        # tamanho de volumes
        from_volumes_to_gids_subvolumes = []

        # ## tamanho de todos subvolumes
        # from_subvolumes_to_gids_faces_of_subvolumes = []

        n = 0

        for volume in volumes:

            nodes_volume = nodes_volumes[volume]
            faces_nodes_volume = faces_nodes[nodes_volume]
            faces_volume = faces_volumes[volume]

            local_ids_subvolumes, points_subvolumes = self.create_subvolumes_from_volume_mpfao(
                                                        volume, nodes_volume,
                                                        faces_nodes_volume, faces_volume
                                                      )
            n_subvolumes = local_ids_subvolumes.shape[0]
            from_volumes_to_gids_subvolumes.append(np.arange(n, n+n_subvolumes))
            [points_of_subvolumes.append(pt) for pt in points_subvolumes]
            n += n_subvolumes

        gids_subvolumes = np.arange(n)
        from_volumes_to_gids_subvolumes = np.array(from_volumes_to_gids_subvolumes)
        points_of_subvolumes = np.array(points_of_subvolumes)

        return gids_subvolumes, from_volumes_to_gids_subvolumes, points_of_subvolumes

    def create_subvolumes_from_volume_mpfao(self, volume, nodes_volume, faces_nodes_volume, faces_volume):

        points_subvolumes = []

        for i, node in enumerate(nodes_volume):
            faces_node = faces_nodes_volume[i]

            # points_subvolume = self.create_subvolume_mpfao(volume, node, faces_node, faces_volume)
            # points_subvolumes.append(points_subvolume)
            points_subvolumes.append(
                self.create_subvolume_mpfao(volume, node, faces_node, faces_volume)
            )

        points_subvolumes = np.array(points_subvolumes)
        local_ids_subvolumes = np.arange(len(points_subvolumes))

        return local_ids_subvolumes, points_subvolumes

    def create_subvolume_mpfao(self, volume, node, faces_node, faces_volume):

        dt = [(self.name_cent_volume, self.type_int), (self.name_cent_node, self.type_int)]
        points = [np.array([volume]), np.array([node])]

        faces_intersect = np.intersect1d(faces_node, faces_volume)
        for i, f1 in enumerate(faces_intersect):
            dt.append((self.name_cent_face + '_' + str(i), self.type_int))
            points.append(np.array([f1]))

        points_subvolume = np.zeros(1, dtype=dt)
        for i, pt in enumerate(points):
            points_subvolume[0][i] = pt

        return points_subvolume

    def get_points_from_st_subvolume_mpfao(self, st_points_subvolume, centroid_faces, centroid_volumes, centroid_nodes):

        

        pass
