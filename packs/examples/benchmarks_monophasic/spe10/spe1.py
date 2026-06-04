from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths
from packs.examples.benchmarks_monophasic.spe10.others import set_permeability
from packs.manager import MeshData, MeshProperty, BoundaryConditions
from packs.examples.benchmarks_monophasic.cross.test_cross_2 import set_weights_nodes, set_fine_transmissibility_without_bc_v2, set_fine_transmissibility_v2
from packs.multiscale.msrsb import create_msrsb_structure as cms
from packs import defpaths
from packs.multiscale.transmissibility_correction.algorithimic_monotone import AlgorithimicMonotone

import os
import numpy as np
import anndata as ad
import scipy.sparse as sp
import shutil

def load_data(layer, nCr):
    rel_path = os.path.join(defpaths.remove_folder, f'layer_{layer}_Cr{nCr}')
    
    T_bc = sp.load_npz(os.path.join(rel_path, 'T_bc.npz'))
    b = np.load(os.path.join(rel_path, 'b.npy'))
    OP = sp.load_npz(os.path.join(rel_path, 'OP.npz'))
    OR = sp.load_npz(os.path.join(rel_path, 'OR.npz'))
    resp = {
        'T': T_bc,
        'b': b,
        'OP': OP,
        'OR': OR
    }
    return resp
    
    

def save_data(T_bc, b_bc, OP, OR, layer, nCr):
    rel_path = os.path.join(defpaths.remove_folder, f'layer_{layer}_Cr{nCr}')
    if os.path.exists(rel_path):
        shutil.rmtree(rel_path)
        
    os.makedirs(rel_path, exist_ok=True)
    sp.save_npz(os.path.join(rel_path, 'T_bc.npz'), T_bc)
    np.save(os.path.join(rel_path, 'b.npy'), b_bc)
    sp.save_npz(os.path.join(rel_path, 'OP.npz'), OP)
    sp.save_npz(os.path.join(rel_path, 'OR.npz'), OR)

def set_boundary_conditions(fine_properties: MeshProperty):
    bc = BoundaryConditions()

    nodes_centroids = fine_properties['nodes_centroids']
    faces = fine_properties['faces']
    faces_centroids = fine_properties['faces_centroids']

    xmin, ymin = nodes_centroids.min(axis=0)
    xmax, ymax = nodes_centroids.max(axis=0)

    # c_p1 = np.array([xmin, ymin])
    # c_p0 = np.array([xmax, ymax])

    c_p1 = np.array([xmin, ymax])
    c_p0 = np.array([xmax, ymin])

    dists = np.linalg.norm(faces_centroids - c_p1, axis=1)
    face_p1 = faces[dists <= dists.min()][0]
    dists[:] = np.linalg.norm(faces_centroids - c_p0, axis=1)
    face_p0 = faces[dists <= dists.min()][0]
    
    faces_pressure = np.array([face_p0, face_p1])
    pressure_presc = np.array([101, 1.0])

    # faces_neumann = np.array([face_p1])
    # neummann_presc_faces = np.array([2.0])

    # areas = fine_properties['areas']
    # edges_dim = fine_properties.edges_dim
    # area_face_p1 = areas[face_p1]

    # import pdb; pdb.set_trace()

    bc.set_boundary('dirichlet_volumes', faces_pressure, pressure_presc)

    # bc.set_boundary('neumann_volumes', faces_neumann, neummann_presc_faces)

    walls_edges = fine_properties['edges'][fine_properties['bool_boundary_edges']]

    edges_values = np.repeat(0.0, walls_edges.shape[0])
    bc.set_boundary('neumann_edges', walls_edges, edges_values)

    fine_properties.insert_or_update_data({
        'neumann_edges': bc['neumann_edges']['id'],
        'neumann_edges_value': bc['neumann_edges']['value']
    })

    # bc.set_boundary('water_saturation_volumes', np.array([face_p1]), np.array([1.0]))
    # bc.set_boundary('water_saturation_edges', np.array([]), np.array([]))

    bc.set_boundary('injectors', np.array([face_p1]), np.array([True]))
    bc.set_boundary('producers', np.array([face_p0]), np.array([True]))

    bc.update_zero_bcs()

    return bc

def plot_perm_layer(fp: MeshProperty, fine_mesh_path: str):
    
    permx = fp['permeability'][:,0,0]
    
    mesh_data = MeshData(mesh_path=fine_mesh_path)
    mesh_data.create_tag('permx')
    mesh_data.insert_tag_data('permx', permx, 'faces')
    mesh_data.export_all_elements_type_to_vtk('teste1_spe', 'faces')

def get_properties():
    
    fine_mesh_path = os.path.join('mesh', 'spe10', 'layer_spe10_perturbada.msh')
    fine_mesh_properties_name = 'spe10_perturbada' 
    fine_mesh_path_v4 = fine_mesh_path

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name, mesh_name_v4=fine_mesh_path_v4)

    return fine_properties, fine_mesh_path


def run4(layer=36, nCr=30):
    
    disjointed = False
    n_levels_adj = 3
    op_tolerance = 0.05
    op_max_it = 10

    fp, fine_mesh_path = get_properties()
    set_permeability(fp, layer=layer)
    bc = set_boundary_conditions(fp)
    set_weights_nodes(fp, update=False)
    transm = set_fine_transmissibility_without_bc_v2(fp)
    T = transm['T_tpfa']
    T_complete = transm['transmissibility_without_bc']
    ag = AlgorithimicMonotone()
    T_for_OP = ag.get_monotone_matrix(T_complete)
    
    resp = set_fine_transmissibility_v2(fp, bc)
    T_bc = resp['transmissibility']
    b_bc = resp['source']
    
    nparts = int(T.shape[0]/nCr)
    primal = cms.create_partition(T, nparts, disjointed=disjointed, nvols_mean=nCr)
    all_regions = cms.create_support_region_and_boundary_v3(T, primal, n_levels_adj=n_levels_adj, ext='_1')
    
    
    OR = cms.get_OR_finite_volume(primal)
    OP, op_iterations, emax = cms.get_msrsb_prolongation_operator_v3(
        all_regions['ind_support'],
        all_regions['ptr_support'],
        T_for_OP,
        OR.transpose(copy=True),
        tol_op=op_tolerance,
        maxit=op_max_it
    )
    
    save_data(T_bc, b_bc, OP, OR, layer, nCr)
    
    
    
    
    
    
    
    
    #  OR1, OP1 = get_OP_and_OR_v3(
    #     'OP1',
    #     primal1,
    #     regions1['ind_support'],
    #     regions1['ptr_support'],
    #     T,
    #     op_tolerance,
    #     matrices_path,
    #     my_params['op_iterations_path'],
    #     my_params['operator_times_path']
    # ) 
    
    
    
    
    
    
    
    
    
    