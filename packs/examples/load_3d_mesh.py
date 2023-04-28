from packs.manager.meshmanager import create_initial_3D_mesh_prperties, MeshProperty, load_mesh_properties
from packs import defpaths
import os
import numpy as np
import scipy.io as io
from packs.utils.utils_old import get_box

# mesh_properties_name = 'exemplo_20x1x1'
# mesh_name = '20x1x1.h5m'

# mesh_properties_name = 'exemplo_30x30x1'
# mesh_name = '30x30x1.h5m'

mesh_properties_name = 'exemplo_21x21x1'
mesh_name = '21x21x1.h5m'


mesh_path = os.path.join(defpaths.mesh, mesh_name)

def define_p1(mesh_properties: MeshProperty):

    n = len(mesh_properties.volumes)
    permeability = np.zeros((n , 3, 3))
    dia = np.arange(3)
    permeability[:, dia, dia] = 1

    xyzmin = mesh_properties.volumes_centroids.min(axis=0)
    xyzmax = mesh_properties.volumes_centroids.max(axis=0)

    xp1 = mesh_properties.volumes[mesh_properties.volumes_centroids[:,0] <= xyzmin[0] + 0.01]
    xp0 = mesh_properties.volumes[mesh_properties.volumes_centroids[:,0] >= xyzmax[0] - 0.01]

    values_p1 = np.repeat(1, len(xp1))
    values_p0 = np.repeat(0, len(xp0))

    saturations1 = np.repeat(1, len(xp1))

    volumes_pressure_defined = np.concatenate([xp1, xp0])
    pressure_defined_values = np.concatenate([values_p1, values_p0])

    volumes_saturation_defined = xp1
    saturation_defined_values = saturations1

    porosity = np.repeat(0.2, n)

    bool_internal_faces = ~mesh_properties.bool_boundary_faces
    adj_internal_faces = mesh_properties.volumes_adj_by_faces[bool_internal_faces].copy()
    upwind = np.full(adj_internal_faces.shape, False, dtype=bool)
    test = np.isin(adj_internal_faces, xp1)
    test1 = (test[:,0] == True) & (test[:,1] == True)
    upwind[test1, 0] = True
    test2 = ~test1
    upwind[test2,0] = True
    assert upwind.sum() == upwind.shape[0] 
    
    mesh_properties.insert_data({
        'volumes_perm': permeability,
        'porosity': porosity,
        'volumes_pressure_defined': volumes_pressure_defined + 1,
        'pressure_defined_values': pressure_defined_values,
        'volumes_saturation_defined': volumes_saturation_defined + 1,
        'saturation_defined_values': saturation_defined_values,
        'upwind': upwind
    })

def define_p2(mesh_properties: MeshProperty):
    """
        malha 20x1x1.h5m
    """

    n = len(mesh_properties.volumes)
    permeability = np.zeros((n , 3, 3))
    dia = np.arange(3)
    permeability[:, dia, dia] = 1

    limites = np.array([
        [0, 0, 0],
        [10, 1, 1]
    ])

    indexes = get_box(mesh_properties.volumes_centroids, limites)
    for index in indexes:
        permeability[index, dia, dia] = 0.01

    xyzmin = mesh_properties.volumes_centroids.min(axis=0)
    xyzmax = mesh_properties.volumes_centroids.max(axis=0)

    xp1 = mesh_properties.volumes[mesh_properties.volumes_centroids[:,0] <= xyzmin[0] + 0.01]
    xp0 = mesh_properties.volumes[mesh_properties.volumes_centroids[:,0] >= xyzmax[0] - 0.01]

    values_p1 = np.repeat(100, len(xp1))
    values_p0 = np.repeat(0, len(xp0))

    saturations1 = np.repeat(1, len(xp1))

    volumes_pressure_defined = np.concatenate([xp1, xp0])
    pressure_defined_values = np.concatenate([values_p1, values_p0])

    volumes_saturation_defined = xp1
    saturation_defined_values = saturations1

    injectors = xp1
    producers = xp0

    porosity = np.repeat(0.2, n)

    bool_internal_faces = ~mesh_properties.bool_boundary_faces
    adj_internal_faces = mesh_properties.volumes_adj_by_faces[bool_internal_faces].copy()
    upwind = np.full(adj_internal_faces.shape, False, dtype=bool)
    test = np.isin(adj_internal_faces, xp1)
    test1 = (test[:,0] == True) & (test[:,1] == True)
    upwind[test1, 0] = True
    test2 = ~test1
    upwind[test2,0] = True
    assert upwind.sum() == upwind.shape[0] 
    
    mesh_properties.insert_data({
        'volumes_perm': permeability,
        'porosity': porosity,
        'volumes_pressure_defined': volumes_pressure_defined + 1,
        'pressure_defined_values': pressure_defined_values,
        'volumes_saturation_defined': volumes_saturation_defined + 1,
        'saturation_defined_values': saturation_defined_values,
        'upwind': upwind,
        'injectors': injectors+1,
        'producers': producers+1
    })

def define_p3(mesh_properties: MeshProperty):
    """
        malha 30x30x1.h5m
    """

    n = len(mesh_properties.volumes)
    permeability = np.zeros((n , 3, 3))
    dia = np.arange(3)
    permeability[:, dia, dia] = 1

    # limites = np.array([
    #     [9, 9, 0],
    #     [19, 19, 1]
    # ])

    # indexes = get_box(mesh_properties.volumes_centroids, limites)
    # for index in indexes:
    #     permeability[index, dia, dia] = 0.01

    limites_xp1 = np.array([
        [0, 0, 0],
        [1, 1, 1]
    ])
    
    indexes_xp1 = get_box(mesh_properties.volumes_centroids, limites_xp1)
    xp1 = mesh_properties.volumes[indexes_xp1]

    limites_xp0 = np.array([
        [29, 29, 0],
        [30, 30, 1]
    ])
    indexes_xp0 = get_box(mesh_properties.volumes_centroids, limites_xp0)
    xp0 = mesh_properties.volumes[indexes_xp0]

    values_p1 = np.repeat(100, len(xp1))
    values_p0 = np.repeat(0, len(xp0))

    saturations1 = np.repeat(1, len(xp1))

    volumes_pressure_defined = np.concatenate([xp1, xp0])
    pressure_defined_values = np.concatenate([values_p1, values_p0])

    volumes_saturation_defined = xp1
    saturation_defined_values = saturations1

    injectors = xp1
    producers = xp0

    porosity = np.repeat(0.2, n)

    bool_internal_faces = ~mesh_properties.bool_boundary_faces
    adj_internal_faces = mesh_properties.volumes_adj_by_faces[bool_internal_faces].copy()
    upwind = np.full(adj_internal_faces.shape, False, dtype=bool)
    test = np.isin(adj_internal_faces, xp1)
    test1 = (test[:,0] == True) & (test[:,1] == True)
    upwind[test1, 0] = True
    test2 = ~test1
    upwind[test2,0] = True
    assert upwind.sum() == upwind.shape[0] 
    
    mesh_properties.insert_data({
        'volumes_perm': permeability,
        'porosity': porosity,
        'volumes_pressure_defined': volumes_pressure_defined + 1,
        'pressure_defined_values': pressure_defined_values,
        'volumes_saturation_defined': volumes_saturation_defined + 1,
        'saturation_defined_values': saturation_defined_values,
        'upwind': upwind,
        'injectors': injectors+1,
        'producers': producers+1
    })

def define_pp2(mesh_properties: MeshProperty):
    """
        malha 21x21x1.h5m
    """

    n = len(mesh_properties.volumes)
    permeability = np.zeros((n , 3, 3))
    dia = np.arange(3)
    permeability[:, dia, dia] = 1

    limites = np.array([
        [7, 7, 0],
        [14, 14, 1]
    ])

    indexes = get_box(mesh_properties.volumes_centroids, limites)
    for index in indexes:
        permeability[index, dia, dia] = 0.001

    limites_xp1 = np.array([
        [0, 0, 0],
        [1, 1, 1]
    ])
    
    indexes_xp1 = get_box(mesh_properties.volumes_centroids, limites_xp1)
    xp1 = mesh_properties.volumes[indexes_xp1]

    limites_xp0 = np.array([
        [20, 20, 0],
        [21, 21, 1]
    ])
    indexes_xp0 = get_box(mesh_properties.volumes_centroids, limites_xp0)
    xp0 = mesh_properties.volumes[indexes_xp0]

    values_p1 = np.repeat(100, len(xp1))
    values_p0 = np.repeat(0, len(xp0))

    saturations1 = np.repeat(0.8, len(xp1))

    volumes_pressure_defined = np.concatenate([xp1, xp0])
    pressure_defined_values = np.concatenate([values_p1, values_p0])

    volumes_saturation_defined = xp1
    saturation_defined_values = saturations1

    injectors = xp1
    producers = xp0

    porosity = np.repeat(0.2, n)

    bool_internal_faces = ~mesh_properties.bool_boundary_faces
    adj_internal_faces = mesh_properties.volumes_adj_by_faces[bool_internal_faces].copy()
    upwind = np.full(adj_internal_faces.shape, False, dtype=bool)
    test = np.isin(adj_internal_faces, xp1)
    test1 = (test[:,0] == True) & (test[:,1] == True)
    upwind[test1, 0] = True
    test2 = ~test1
    upwind[test2,0] = True
    assert upwind.sum() == upwind.shape[0] 
    
    mesh_properties.insert_data({
        'volumes_perm': permeability,
        'porosity': porosity,
        'volumes_pressure_defined': volumes_pressure_defined + 1,
        'pressure_defined_values': pressure_defined_values,
        'volumes_saturation_defined': volumes_saturation_defined + 1,
        'saturation_defined_values': saturation_defined_values,
        'upwind': upwind,
        'injectors': injectors+1,
        'producers': producers+1
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

def calculate_areas(mesh_properties: MeshProperty):
    bool_internal_faces = ~mesh_properties.bool_boundary_faces
    areas = np.zeros(bool_internal_faces.sum())
    unitary_vector_internal = np.zeros((bool_internal_faces.sum(), 3))
    for i, face in enumerate(mesh_properties.faces[bool_internal_faces]):
        faces_nodes = mesh_properties.nodes_of_faces[face]
        nodes_points = mesh_properties.nodes_centroids[faces_nodes]
        centroid = np.mean(nodes_points, axis=0)
        minp = nodes_points.min(axis=0)
        maxp = nodes_points.max(axis=0)
        diag = np.linalg.norm(maxp - minp)
        area = (diag/np.sqrt(2))**2
        areas[i] = area
        v1 = centroid - nodes_points[0]
        v2 = centroid - nodes_points[1]
        normal = np.cross(v1, v2)
        unormal = normal/np.linalg.norm(normal)
        volumes_adj = mesh_properties.volumes_adj_by_faces[face]
        vtest = mesh_properties.volumes_centroids[volumes_adj[1]] - centroid
        test = np.dot(unormal, vtest)
        if test > 0:
            pass
        else:
            unormal = -1*unormal
        unitary_vector_internal[i,:] = unormal

    
    mesh_properties.insert_data({
        'internal_areas': areas,
        'unitary_vector_internal': unitary_vector_internal
    })

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

    # ifaces = mesh_properties.faces[~mesh_properties.bool_boundary_faces]
    # bool_ifaces_per_volume = np.isin(mesh_properties.faces_of_volumes, ifaces)
    
    mesh_properties.remove_datas(
        [
            'bool_boundary_faces',
            'volumes_adj_by_faces',
            'faces',
            'edges',
            'volumes_adj_by_nodes',
            'nodes_of_faces',
            'nodes_centroids',
            'nodes_of_volumes',
            'nodes',
            'faces_of_volumes'
        ]
    )

    mesh_properties.update_data({
        'volumes': mesh_properties.volumes + 1
    })

    internal_faces_remap = internal_faces_remap + 1
    adj_internal_faces = adj_internal_faces + 1
    
    mesh_properties.insert_data({
        'internal_faces': internal_faces_remap,
        'adjacencies': adj_internal_faces,
        'centroids_internal_faces': centroids_internal_faces,
        'h_dist': h_dist
        # 'bool_ifaces_per_volume': bool_ifaces_per_volume
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
    calculate_areas(mesh_properties)
    mesh_properties.export_data()

def thr_func():
    
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    # define_p1(mesh_properties)
    # define_p2(mesh_properties)
    # define_p3(mesh_properties)
    define_pp2(mesh_properties)
    update_properties(mesh_properties)
    mesh_properties.export_data()

def save_mat_file():
    
    mesh_properties: MeshProperty = load_mesh_properties(mesh_properties_name)
    # matlab_data_path = os.path.join(defpaths.flying, 'variables_p1.mat')
    matlab_data_path = os.path.join(defpaths.flying, 'variables_p2.mat')
    io.savemat(matlab_data_path, mesh_properties.get_all_data())


first_func()
second_func()
thr_func()
save_mat_file()


mesh_properties = load_mesh_properties(mesh_properties_name)
import pdb; pdb.set_trace()
print(mesh_properties)
