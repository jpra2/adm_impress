from packs.tpfa.biphasic import TpfaBiphasicCons
from packs.tpfa.monophasic import TpfaMonophasic
from packs.data_class import GeometricData, Elements
from packs.preprocess import TpfaPreprocess
from packs.running.initial_mesh_properties import initial_mesh
import numpy as np
import pdb


M, _, data_impress, wells = initial_mesh()
perms = M.permeability[:]
perms = perms.reshape(len(perms), 3, 3)
saturation = M.saturation[:]

elements = Elements()
elements.set_mesh_elements(
    M.volumes.all,
    M.faces.all,
    M.edges.all,
    M.nodes.all,
    M.faces.boundary,
    M.nodes.boundary
)
elements.create_adj_matrix_volumes_to_faces(elements.volumes, elements.faces, M.volumes.bridge_adjacencies(elements.volumes, 3, 2))
elements.create_adj_matrix_faces_to_edges(elements.faces, elements.edges, M.faces.bridge_adjacencies(elements.faces, 2, 1))
elements.create_adj_matrix_edges_to_nodes(elements.edges, elements.nodes, M.edges.bridge_adjacencies(elements.edges, 1, 0))
elements.insert_value('volumes_adj_internal_faces', elements.faces_to_volumes(elements.internal_faces), 'array')
n_vols = len(elements.volumes)

preprocess = TpfaPreprocess()
geom = GeometricData()
geom['abs_u_normal_faces'] = M.faces.normal
geom['block_dimension'] = np.repeat(1.0, n_vols*3).reshape(n_vols, 3)
geom['volume'] = np.prod(geom['block_dimension'], axis=1)
geom['centroid_volumes'] = M.volumes.center[:]
geom['centroid_faces'] = M.faces.center[:]
geom['centroid_edges'] = M.edges.center[:]
geom['centroid_nodes'] = M.nodes.center[:]
nodes_of_faces = elements.faces_to_nodes(elements.faces)
points_faces = geom['centroid_nodes'][nodes_of_faces]
geom['areas'] = preprocess.get_areas_faces(elements.faces, points_faces)
geom['unit_direction_internal_faces'] = preprocess.get_u_normal_internal_faces(elements.get_value('volumes_adj_internal_faces'), geom['abs_u_normal_faces'][elements.internal_faces], geom['centroid_volumes'])

pdb.set_trace()

monophasic = TpfaMonophasic()
biphasic = TpfaBiphasicCons()
###
##initialize
biphasic.set_current_saturation(saturation)
###


pdb.set_trace()
