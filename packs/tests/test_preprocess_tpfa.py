from packs.preprocess import TpfaPreprocess
from packs.load.preprocessor0 import M
from packs.data_class import Elements
import pdb
import numpy as np
from packs.data_class import GeometricData

elements = Elements()
elements.set_mesh_elements(
    volumes=M.volumes.all,
    faces=M.faces.all,
    edges=M.edges.all,
    nodes=M.nodes.all,
    boundary_faces=M.faces.boundary,
    boundary_nodes=M.nodes.boundary
)

elements.create_adj_matrix_volumes_to_faces(elements.volumes, elements.faces, M.volumes.bridge_adjacencies(elements.volumes, 3, 2))
elements.create_adj_matrix_faces_to_edges(elements.faces, elements.edges, M.faces.bridge_adjacencies(elements.faces, 2, 1))
elements.create_adj_matrix_edges_to_nodes(elements.edges, elements.nodes, M.edges.bridge_adjacencies(elements.edges, 1, 0))

elements.insert_information('v0_int_faces', M.faces.bridge_adjacencies(elements.internal_faces, 2, 3), tipo='array')
elements.insert_information('v0_bound_faces', M.faces.bridge_adjacencies(elements.boundary_faces, 2, 3), tipo='array')


n_vols = len(elements.volumes)
geom = GeometricData()
geom['unit_normal_vector'] = M.faces.normal
geom['block_dimension'] = np.repeat(1.0, n_vols*3).reshape(n_vols, 3)
geom['centroid_volumes'] = M.volumes.center[:]
geom['centroid_faces'] = M.faces.center[:]
geom['centroid_edges'] = M.edges.center[:]
geom['centroid_nodes'] = M.nodes.center[:]

M.permeability[:] = np.array([1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0])
perms = M.permeability[:]
perms = perms.reshape(len(perms), 3, 3)

preprocess = TpfaPreprocess()

keq = preprocess.get_equivalent_permeability_faces_from_diagonal_permeability(
    elements.volumes,
    elements.faces,
    elements.get_value('v0_int_faces'),
    elements.get_value('v0_bound_faces'),
    geom['unit_normal_vector'],
    perms,
    elements.internal_faces,
    geom['block_dimension']
)

nodes_of_faces = elements.faces_to_nodes(elements.faces)
points_faces = geom['centroid_nodes'][nodes_of_faces]
geom['areas'] = preprocess.get_areas_faces(elements.faces, points_faces)




pdb.set_trace()
