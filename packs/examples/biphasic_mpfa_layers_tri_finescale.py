from packs.biphasic.relative_perm.brooks_and_corey import BrooksAndCorey
from packs.biphasic.mobility import BiphasicMobility
from packs.biphasic.unstructured.mobility_mesh_elements import direct_edges_mobility
from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths
from packs.manager import MeshProperty, MeshData, BoundaryConditions, SimulationData
from packs.multiscale.unstructured.test.test_cross import set_weights_nodes, set_fine_transmissibility
# from packs.multiscale.unstructured.test.test_brazil import define_faces_in_losangle, set_permeability
from packs.mpfa_methods.flux_calculation.lsds_method import LsdsFluxCalculation
from packs.mpfa_methods.weight_interpolation.gls_weight_2d import get_gls_nodes_weights

import os
import numpy as np
from typing import Tuple
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

def get_properties():

    fine_mesh_properties_name = 'finescale_layers'
    fine_mesh_path = defpaths.finescale_mesh_layer

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name)

    return fine_properties, fine_mesh_path

def get_nodes_org_from_faces_of_nodes_object(fp: MeshProperty):
    n_faces_of_nodes = fp['n_faces_of_nodes']

    all_nodes_org = []
    faces_of_nodes_org = []
    n_nodes = np.arange(1, n_faces_of_nodes.max()+1)

    for i in n_nodes:
        test = n_faces_of_nodes == i
        v4 = fp.nodes[test]
        all_nodes_org.append(v4)
    
        ft = fp['faces_of_nodes'][v4].copy()
        ft2 = np.concatenate(ft)
        ft3 = ft2.reshape((ft.shape[0], i))
        faces_of_nodes_org.append(ft3)
    
    all_nodes_org = np.array(all_nodes_org, dtype='O')
    faces_of_nodes_org = np.array(faces_of_nodes_org, dtype='O')

    return all_nodes_org, faces_of_nodes_org, n_nodes

def get_internal_nodes_org_from_faces_of_nodes_object(fp: MeshProperty):
    n_faces_of_nodes = fp['n_faces_of_nodes']
    internal_nodes = fp.internal_nodes

    all_nodes_org = []
    faces_of_nodes_org = []
    n_nodes = np.arange(2, n_faces_of_nodes.max()+1)
    new_n_nodes = []

    for i in n_nodes:
        test = n_faces_of_nodes == i
        v4 = fp.nodes[test]
        v4 = np.intersect1d(v4, internal_nodes)
        if v4.shape[0] == 0:
            continue
        all_nodes_org.append(v4)
    
        ft = fp['faces_of_nodes'][v4].copy()
        ft2 = np.concatenate(ft)
        ft3 = ft2.reshape((ft.shape[0], i))
        faces_of_nodes_org.append(ft3)
        new_n_nodes.append(i)
    
    all_nodes_org = np.array(all_nodes_org, dtype='O')
    faces_of_nodes_org = np.array(faces_of_nodes_org, dtype='O')
    new_n_nodes = np.array(new_n_nodes)
    
    return all_nodes_org, faces_of_nodes_org, new_n_nodes

def set_boundary_conditions(fine_properties: MeshProperty):
    bc = BoundaryConditions()

    nodes_centroids = fine_properties['nodes_centroids']
    faces = fine_properties['faces']
    faces_centroids = fine_properties['faces_centroids']

    xmin, ymin = nodes_centroids.min(axis=0)
    xmax, ymax = nodes_centroids.max(axis=0)

    c_p1 = np.array([xmin, ymin])
    c_p0 = np.array([xmax, ymax])

    dists = np.linalg.norm(faces_centroids - c_p1, axis=1)
    face_p1 = faces[dists <= dists.min()][0]
    dists[:] = np.linalg.norm(faces_centroids - c_p0, axis=1)
    face_p0 = faces[dists <= dists.min()][0]
    
    faces_pressure = np.array([face_p0])
    pressure_presc = np.array([0.0])

    faces_neumann = np.array([face_p1])
    neummann_presc_faces = np.array([1.0])

    bc.set_boundary('dirichlet_volumes', faces_pressure, pressure_presc)

    bc.set_boundary('neumann_volumes', faces_neumann, neummann_presc_faces)

    walls_edges = fine_properties.boundary_edges

    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    fine_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    bc.set_boundary('water_saturation_volumes', np.array([face_p1]), np.array([1.0]))
    bc.set_boundary('water_saturation_edges', np.array([]), np.array([]))

    bc.set_boundary('injectors', np.array([face_p1]), np.array([True]))
    bc.set_boundary('producers', np.array([face_p0]), np.array([True]))

    bc.update_zero_bcs()

    return bc


def run5():

    
    dt = 0.00005
    max_vpi = 1.3
    loop = 0
    max_loop = np.inf
    load = False
    loop_intervals = 10

    cumulative_oil = 0.0
    cumulative_water = 0.0
    vpi = 0.0

    relative_perm = BrooksAndCorey()
    biphasic_mobility = BiphasicMobility()
    lsds = LsdsFluxCalculation()
    simulation_data = SimulationData('biphasic_fine_mesh')

    fp, fine_mesh_path = get_properties()

    nodes_org, faces_of_nodes_org, n_nodes_org = fp.get_internal_nodes_org_from_faces_of_nodes_object()






    import pdb; pdb.set_trace()




    pass

