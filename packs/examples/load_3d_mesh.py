from packs.manager.meshmanager import create_initial_3D_mesh_prperties, MeshProperty, load_mesh_properties
from packs import defpaths
import os
import numpy as np
import scipy.io as io

mesh_properties_name = 'exemplo_estruturado_3d'
mesh_name = 'cube_structured.msh'
mesh_path = os.path.join(defpaths.mesh, mesh_name)

# mesh_properties = create_initial_3D_mesh_prperties(mesh_path, mesh_properties_name)
# mesh_properties.export_data()
# import pdb; pdb.set_trace()
# mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)


def calculate_volumes_structured(mesh_properties: MeshProperty):
    
    volume_of_volumes = np.zeros(len(mesh_properties.volumes))
    volumes_centroids = np.zeros((len(mesh_properties.volumes), 3))

    for volume in mesh_properties.volumes:
        nodes_of_volume = mesh_properties.nodes_of_volumes[volume]
        vcentroids = mesh_properties.nodes_centroids[nodes_of_volume]
        minp = vcentroids.min(axis=0)
        maxp = vcentroids.max(axis=0)
        diagonal = np.linalg.norm(maxp - minp)
        vol = (diagonal**3)/(3*np.sqrt(3))
        centroid = np.mean(vcentroids, axis=0)
        volume_of_volumes[volume] = vol
        volumes_centroids[volume] = centroid
    
    return volume_of_volumes, volumes_centroids

# # ######################################
# mesh_properties = load_mesh_properties(mesh_properties_name)
# vol_volumes, volumes_centroids = calculate_volumes_structured(mesh_properties)
# mesh_properties.insert_data({
#     'vol_volumes': vol_volumes,
#     'volumes_centroids': volumes_centroids

# })
# mesh_properties.export_data()
# import pdb; pdb.set_trace()
# ######################################



def update_properties(mesh_properties: MeshProperty):

    internal_faces = mesh_properties.faces[~mesh_properties.bool_boundary_faces]
    adj_internal_faces = mesh_properties.volumes_adj_by_faces[internal_faces]

    internal_faces_remap = np.arange(len(internal_faces))

    

mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)



mesh_properties.remove_datas(
    [
        'bool_boundary_faces'
    ]
)

import pdb; pdb.set_trace()

matlab_data_path = os.path.join(defpaths.flying, 'variables_p1.mat')
io.savemat(matlab_data_path, mesh_properties.get_all_data())

import pdb; pdb.set_trace()
print(mesh_properties)
