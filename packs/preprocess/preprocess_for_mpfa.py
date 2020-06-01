import pdb
from packs.data_class.data_manager import DataManager
from packs.directories import only_mesh_name
import numpy as np

class PreprocessForMpfa(DataManager):

    def __init__(self, data_impress, data_name='PreprocessForMpfa', load=False):
        data_name = data_name + '_' + only_mesh_name + '.npz'
        super().__init__(data_name=data_name, load=load)
        self.centroid_faces = data_impress['centroid_faces']
        self.centroid_nodes = data_impress['centroid_nodes']

    def init_data(self, M, elements_lv0, data_impress):
        n_faces = len(elements_lv0['faces'])
        n_internal_faces = len(elements_lv0['internal_faces'])
        n_volumes = len(elements_lv0['volumes'])
        n_nodes = len(elements_lv0['nodes'])
        n_edges = len(elements_lv0['edges'])

        self.create_midpoints_edges(M, elements_lv0, data_impress)
        self.create_subfaces(M, elements_lv0, data_impress)

    def create_midpoints_edges(self, M, elements_lv0, data_impress):
        self['midpoints_edges'] = np.zeros((len(elements_lv0['edges']), 3))

        nodes_edges = M.edges.bridge_adjacencies(elements_lv0['edges'], 0, 0)
        cents_nodes_edges = data_impress['centroid_nodes'][nodes_edges]
        self['midpoints_edges'][:] = np.mean(cents_nodes_edges, axis=1)

    def create_subfaces(self, M, elements_lv0, data_impress):
        '''
            define os vertices que pertencem as subfaces
        '''
        edges_faces = elements_lv0['faces_edge_edges']
        centroid_faces = data_impress['centroid_faces']
        centroid_nodes = data_impress['centroid_nodes']
        nodes_faces = elements_lv0['faces_node_nodes']
        midpoints_edges = self['midpoints_edges']
        edges = elements_lv0['edges']

        todas_normais = []
        todos_pontos = []
        local_ids_subfaces = []

        for face in elements_lv0['faces']:
            cent_face = centroid_faces[face]
            edges = edges_faces[face]
            mid_edges = midpoints_edges[edges]
            nodes = nodes_faces[face]
            cent_nodes = centroid_nodes[nodes]

            vectors_midpoints = mid_edges - cent_face
            vectors_nodes = cent_nodes - cent_face
            norms_vectors_midpoins = np.linalg.norm(vectors_midpoints, axis=1)
            norms_vectors_nodes = np.linalg.norm(vectors_nodes, axis=1)
            unit_sentido_test = np.cross(vectors_nodes[0], vectors_midpoints[0])
            unit_sentido_test = unit_sentido_test/np.linalg.norm(unit_sentido_test)

            # pdb.set_trace()

            normais = []
            points_areas = []

            for i, vec1 in enumerate(vectors_nodes):
                norm1 = norms_vectors_nodes[i]
                angles = []
                all_norms = []
                for j, vec2 in enumerate(vectors_midpoints):
                    norm2 = norms_vectors_midpoins[j]
                    angle = (180/np.pi)*np.arccos(np.dot(vec1, vec2)/(norm1 + norm2))
                    normal = np.cross(vec1, vec2)
                    unit2 = normal/np.linalg.norm(normal)
                    if np.allclose(unit_sentido_test, unit2):
                        pass
                    else:
                        angle = 360 - angle
                        normal *= -1
                    angles.append(angle)
                    all_norms.append(normal)

                angles = np.array(angles)
                ids = np.arange(len(angles))

                imin = ids[angles == angles.min()][0]
                imax = ids[angles == angles.max()][0]

                midp1 = mid_edges[imin]
                midp2 = mid_edges[imax]

                ######################
                # pontos da subface
                # face = id da face cujo centroide == cent_face
                # id_midp1 = id do edge cujo midpoint == midp1
                # id_midp2 = id do edge cujo midpoint == midp2
                # id_node = id do node cujo centroid_nodes == cent_nodes[i]
                points = np.array([cent_face, midp1, cent_nodes[i], midp2])
                # id_midp1 = edges[imin]
                # id_node = nodes[i]
                # id_midp2 = edges[imax]
                # points = np.array([face, id_midp1, id_node, id_midp2])
                # points3 = np.array([centroid_faces[face], midpoints_edges[edges[imin]], centroid_nodes[nodes[i]], midpoints_edges[edges[imax]]])
                # points = points3
                ################################

                normais.append((all_norms[imax] + all_norms[imin])/2)
                points_areas.append(points)

            normais = np.array(normais)
            points_areas = np.array(points_areas)
            local_ids_subfaces.append(np.arange(len(normais)))

            todos_pontos.append(points_areas)
            todas_normais.append(normais)

        todas_normais = np.array(todas_normais)
        todos_pontos = np.array(todos_pontos)
        local_ids_subfaces = np.array(local_ids_subfaces)

        self['normals_subfaces'] = todas_normais
        self['points_subfaces'] = todos_pontos
        self['local_ids_subfaces'] = local_ids_subfaces
