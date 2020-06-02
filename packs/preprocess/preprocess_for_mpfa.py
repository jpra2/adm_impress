import pdb
from packs.data_class.data_manager import DataManager
from packs.directories import only_mesh_name
import numpy as np

class PreprocessForMpfa(DataManager):

    def __init__(self, data_name='PreprocessForMpfa', load=False):
        data_name = data_name + '_' + only_mesh_name + '.npz'
        super().__init__(data_name=data_name, load=load)

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

    def create_subfaces(self, faces, edges, edges_faces, nodes_faces, centroid_faces, centroid_nodes, midpoints_edges):
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

            normais, points_areas, local_id_subfaces = self.func1(cent_face, edges, mid_edges, nodes, cent_nodes)

            todos_pontos.append(points_areas)
            todas_normais.append(normais)
            local_ids_subfaces.append(local_id_subfaces)

        todas_normais = np.array(todas_normais)
        todos_pontos = np.array(todos_pontos)
        local_ids_subfaces = np.array(local_ids_subfaces)

        return todas_normais, todos_pontos, local_ids_subfaces

    def func1(self, cent_face, edges, midpoint_edges, nodes, cent_nodes):
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

            normais.append((all_norms[imax] + all_norms[imin])/2)
            points_areas.append(points)

        normais = np.array(normais)
        points_areas = np.array(points_areas)
        local_id_subfaces = np.arange(len(normais))

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
