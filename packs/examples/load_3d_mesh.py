from packs.manager.meshmanager import create_initial_3D_mesh_prperties, MeshProperty, load_mesh_properties
from packs import defpaths
import os
import numpy as np
import scipy.io as io

mesh_properties_name = 'exemplo_estruturado_3d'
mesh_name = 'cube_structured.msh'
mesh_path = os.path.join(defpaths.mesh, mesh_name)

def define_p1(mesh_properties: MeshProperty):

    n = len(mesh_properties.volumes)
    permeability = np.zeros((n , 3, 3))
    dia = np.arange(3)
    permeability[:, dia, dia] = 1

    mesh_properties.insert_data({
        'volumes_perm': permeability
    })


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

def calculate_areas(mesh_properties):
    for face in mesh_properties.faces:
        pass

def update_properties(mesh_properties: MeshProperty):

    internal_faces = mesh_properties.faces[~mesh_properties.bool_boundary_faces]
    adj_internal_faces = mesh_properties.volumes_adj_by_faces[internal_faces].copy()
    internal_faces_remap = np.arange(len(internal_faces))

    h_dist = np.zeros((len(internal_faces), 2))
    centroids_internal_faces = np.zeros((len(internal_faces), 3))

    for i, face in enumerate(internal_faces):
        nodes_of_face = mesh_properties.nodes_of_faces[face]
        vcentroids = mesh_properties.nodes_centroids[nodes_of_face]
        centroid_face = np.mean(vcentroids, axis=0)
        centroids_internal_faces[i,:] = centroid_face
        volumes_adj = adj_internal_faces[i]
        c1 = mesh_properties.volumes_centroids[volumes_adj[1]]
        c0 = mesh_properties.volumes_centroids[volumes_adj[0]]
        h_dist[i, 0] = np.linalg.norm(centroid_face - c0)
        h_dist[i, 1] = np.linalg.norm(centroid_face - c1)
    
    mesh_properties.remove_datas(
        [
            'bool_boundary_faces',
            'volumes_adj_by_faces',
            'faces',
            'edges',
            'volumes_adj_by_nodes',
            'nodes_of_faces',
            'nodes_centroids',
            'nodes_of_volumes'
        ]
    )



    mesh_properties.insert_data({
        'internal_faces': internal_faces_remap,
        'adjacencies': adj_internal_faces,
        'centroids_internal_faces': centroids_internal_faces,
        'h_dist': h_dist
    })


def first_func():
    mesh_properties = create_initial_3D_mesh_prperties(mesh_path, mesh_properties_name)
    mesh_properties.export_data()

def second_func():

    mesh_properties = load_mesh_properties(mesh_properties_name)
    vol_volumes, volumes_centroids = calculate_volumes_structured(mesh_properties)
    mesh_properties.insert_data({
        'vol_volumes': vol_volumes,
        'volumes_centroids': volumes_centroids

    })
    mesh_properties.export_data()

def thr_func():
    
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    define_p1(mesh_properties)
    update_properties(mesh_properties)
    mesh_properties.export_data()

def save_mat_file():
    
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    matlab_data_path = os.path.join(defpaths.flying, 'variables_p1.mat')
    io.savemat(matlab_data_path, mesh_properties.get_all_data())


first_func()
second_func()
thr_func()
save_mat_file()


mesh_properties = load_mesh_properties(mesh_properties_name)
import pdb; pdb.set_trace()
print(mesh_properties)
