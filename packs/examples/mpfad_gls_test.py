import os
from packs.manager.meshmanager import CreateMeshProperties, MeshProperty, create_initial_mesh_properties, load_mesh_properties
from packs.manager.mesh_data import MeshData
import packs.defpaths as defpaths
import numpy as np
import pandas as pd
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import CalculateGlsWeight2D
from packs.manager.generic_data import VerticesWeight
from packs.utils import calculate_face_properties

mesh_path = os.path.join(defpaths.mesh, '2d_unstructured.msh')
mesh_name = 'mpfad_gls_test'
weight_name = 'weights_gls_test'

##################################################
# mesh_properties = create_initial_mesh_properties(mesh_path, mesh_name)
# mesh_properties.export_data()
# import pdb; pdb.set_trace()
###################################################

######################################################
# mesh_properties: MeshProperty = load_mesh_properties(mesh_name)

# mesh_properties.update_data(
#     {
#       'nodes_centroids': mesh_properties.nodes_centroids[:, 0:2],
#       'faces_centroids': mesh_properties.faces_centroids[:, 0:2]
#       }
# )

# datas_to_rename = {
#     'faces_adj_by_nodes': 'faces_of_nodes',
#     'faces_adj_by_edges': 'adjacencies',
#     'nodes_adj_by_nodes': 'nodes_of_nodes',
#     'edges_adj_by_nodes': 'edges_of_nodes'
# }

# mesh_properties.rename_data(datas_to_rename)

# mesh_properties.export_data()
# import pdb; pdb.set_trace()
#################################################

###################################################
# mesh_properties = load_mesh_properties(mesh_name)

# norma, unitary_normal_edges = calculate_face_properties.create_unitary_normal_edges_xy_plane(
#     mesh_properties.nodes_of_edges,
#     mesh_properties.nodes_centroids,
#     mesh_properties.adjacencies,
#     mesh_properties.faces_centroids,
#     mesh_properties.bool_boundary_edges
# )

# mesh_properties.insert_data({
#     'unitary_normal_edges': unitary_normal_edges
# })

# mesh_properties.export_data()
# import pdb; pdb.set_trace()
######################################

########################################
# mesh_properties = load_mesh_properties(mesh_name)
# permeability = np.zeros((len(mesh_properties.faces), 2, 2))
# permeability[:, [0, 1], [0, 1]] = 1
# mesh_properties.insert_data({'permeability': permeability})
# mesh_properties.export_data()
# import pdb; pdb.set_trace()
########################################



#####################################################
mesh_properties = load_mesh_properties(mesh_name)

weight = VerticesWeight()
weight.insert_name(name=weight_name)

calculate_weight = CalculateGlsWeight2D(**mesh_properties.get_all_data())
nodes_weights = calculate_weight.get_weights_internal_nodes()
import pdb; pdb.set_trace()
##################################333


print(mesh_properties)