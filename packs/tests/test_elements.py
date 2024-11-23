
from packs.data_class import Elements
import pdb
import numpy as np
import unittest

class TestElements(unittest.TestCase):

    def setUp(self):

        from packs.load.preprocessor0 import M
        self.M = M
        self.elements = Elements()
        self.elements.set_mesh_elements(
            M.volumes.all,
            M.faces.all,
            M.edges.all,
            M.nodes.all,
            M.faces.internal,
            M.nodes.internal
        )

        self.elements.create_adj_matrix_volumes_to_faces(elements.volumes, elements.faces, M.volumes.bridge_adjacencies(elements.volumes, 3, 2))
        self.elements.create_adj_matrix_faces_to_edges(elements.faces, elements.edges, M.faces.bridge_adjacencies(elements.faces, 2, 1))
        self.elements.create_adj_matrix_edges_to_nodes(elements.edges, elements.nodes, M.edges.bridge_adjacencies(elements.edges, 1, 0))

    def test_v_to_f(self):
        v_to_f = elements.volumes_to_faces(elements.volumes)
        t1 = self.M.volumes.bridge_adjacencies(elements.volumes, 3, 2)
        assert np.allclose(v_to_f, t1)

    def test_f_to_e(self):
        f_to_e = elements.faces_to_edges(elements.faces)
        t1 = self.M.faces.bridge_adjacencies(elements.faces, 2, 1)
        assert np.allclose(f_to_e, t1)

    def test_e_to_n(self):
        e_to_n = elements.edges_to_nodes(elements.edges)
        t1 = self.M.edges.bridge_adjacencies(elements.edges, 1, 0)
        assert np.allclose(e_to_n, t1)






# elements = Elements()
# elements.set_mesh_elements(
#     M.volumes.all,
#     M.faces.all,
#     M.edges.all,
#     M.nodes.all,
#     M.faces.internal,
#     M.nodes.internal
# )
#
#
# elements.create_adj_matrix_volumes_to_faces(elements.volumes, elements.faces, M.volumes.bridge_adjacencies(elements.volumes, 3, 2))
# elements.create_adj_matrix_faces_to_edges(elements.faces, elements.edges, M.faces.bridge_adjacencies(elements.faces, 2, 1))
# elements.create_adj_matrix_edges_to_nodes(elements.edges, elements.nodes, M.edges.bridge_adjacencies(elements.edges, 1, 0))
#
# v_to_f = elements.volumes_to_faces(elements.volumes)
# t1 = M.volumes.bridge_adjacencies(elements.volumes, 3, 2)
# assert np.allclose(v_to_f, t1)
#
# f_to_e = elements.faces_to_edges(elements.faces)
# t1 = M.faces.bridge_adjacencies(elements.faces, 2, 1)
# assert np.allclose(f_to_e, t1)
#
# e_to_n = elements.edges_to_nodes(elements.edges)
# t1 = M.edges.bridge_adjacencies(elements.edges, 1, 0)
# assert np.allclose(e_to_n, t1)
