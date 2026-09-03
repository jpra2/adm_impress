from packs.mpfa_methods.mesh_preprocess import MpfaPreprocess, preprocess_mesh
from packs import defpaths
from packs.examples.benchmarks_monophasic.spe10.others import set_permeability
# from packs.manager import MeshData, MeshProperty, BoundaryConditions
from packs.manager import MeshProperty, BoundaryConditions
from packs.examples.benchmarks_monophasic.cross.test_cross_2 import set_weights_nodes, set_fine_transmissibility_without_bc_v2, set_fine_transmissibility_v2
from packs.multiscale.msrsb import create_msrsb_structure as cms
from packs import defpaths
from packs.multiscale.transmissibility_correction.algorithimic_monotone import AlgorithimicMonotone
from packs.multiscale.transmissibility_correction.enhanced import Enhanced
from packs.utils.multiscale_methods import print_fine_interfaces_coarse_mesh_2d
from packs.multiscale.preconditioner import multilevel_class as mcl

import os
import numpy as np
import anndata as ad
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve, factorized
import shutil
from pyamg.krylov import fgmres
import time

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
    from packs.manager.mesh_data import MeshData
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

def get_params():
    return {
        'eps': 5e-3,
        'n_levels_adj': 3,
        'tol_op': 0.05,
        'op_max_it': 10,
        'maxiter': int(1e4),
        'tol': 1e-8,
        'restart': 100
    }

def identify_filtered_interfaces(T_filtered: sp.csc_matrix, adjacencies: np.ndarray):
        
    Tdata = sp.find(T_filtered)
    
    interfaces = np.full(adjacencies.shape[0], False, dtype=bool)
    
    
    lines_f = Tdata[0]
    cols_f = Tdata[1]
    v1 = lines_f == cols_f
    v2 = ~v1
    lines_f = lines_f[v2]
    cols_f = cols_f[v2]
    
    for i, j in zip(lines_f, cols_f):
        test1 = (adjacencies[:, 0] == i) & (adjacencies[:, 1] == j)
        test2 = (adjacencies[:, 0] == j) & (adjacencies[:, 1] == i)
        test = test1 | test2
        interfaces [:] = interfaces | test
    
    
    interfaces[:] = ~interfaces
    
    filtered_interfaces = np.arange(interfaces.shape[0])[interfaces]
    
    return filtered_interfaces
    
        
        
def mount_kmatrix(fp: MeshProperty):
    
    
    
    adjacencies = fp['adjacencies']
    internal_edges = fp.internal_edges
    permeability = fp['permeability']
    unitary_normal = fp['unitary_normal_edges']
    faces_centroids = fp['faces_centroids']
    edges_centroids = fp.edges_centroids
    area = fp['areas']
    edges_dim = fp.edges_dim
    faces = fp['faces']
    
    
    n_internal = unitary_normal[internal_edges]

    adj_internal = adjacencies[internal_edges]
    elems_L = adj_internal[:, 0]  # Elementos à esquerda
    elems_R = adj_internal[:, 1]  # Elementos à direita

    # --- 2. Extrair a submatriz K_2x2 para os elementos vizinhos ---
    # Fatiamos os eixos de K para pegar apenas as dimensões x e y (:2, :2)
    K_L = permeability[elems_L, :2, :2]  
    K_R = permeability[elems_R, :2, :2]  
    
    # --- 3. Executar a operação n^T * K * n vetorizada ---
    # kn_L e kn_R terão dimensão (len(nii),) contendo o escalar de cada face interna
    kn_L = np.einsum('id,ide,ie->i', n_internal, K_L, n_internal)
    kn_R = np.einsum('id,ide,ie->i', n_internal, K_R, n_internal)
    
    # --- 1. Dados geométricos das faces internas (Filtrados por nii) ---
    # Substitua pelas variáveis reais da sua malha:
    d_L = np.linalg.norm(faces_centroids[adjacencies[internal_edges, 0]] - edges_centroids[internal_edges], axis=1)          # Dimensão: (len(nii),) - Distância centro_L -> face
    d_R = np.linalg.norm(faces_centroids[adjacencies[internal_edges, 0]] - edges_centroids[internal_edges], axis=1)         # Dimensão: (len(nii),) - Distância centro_R -> face
    area_f = edges_dim[internal_edges]    # Dimensão: (len(nii),) - Comprimento/Área da face
    
    mu = 1.0
    
    numerador_k = d_L + d_R
    denominador_k = (d_L / kn_L) + (d_R / kn_R)
    k_interface = numerador_k / denominador_k
    # k_interface = k_interface * area_f
    # k_interface = k_interface / (d_L + d_R)
    # transmissibilidade = (k_interface * area_f) / (mu * (d_L + d_R))
    
    biedges = internal_edges
    
    l2 = np.concatenate([adjacencies[biedges, 0], adjacencies[biedges, 0], adjacencies[biedges, 1], adjacencies[biedges, 1]])
    c2 = np.concatenate([adjacencies[biedges, 0], adjacencies[biedges, 1], adjacencies[biedges, 1], adjacencies[biedges, 0]])
    d2 = np.concatenate([k_interface,   -k_interface,  k_interface,  -k_interface])
    T_tpfa = sp.csc_matrix((d2,(l2,c2)), shape=(faces.shape[0],faces.shape[0]))
    
    return {'T_matrix': T_tpfa}
    
    

    
    
    
    
    
    

def run4(layer=36, nCr=81):
    
    disjointed = False
    n_levels_adj = 3
    my_params = get_params()

    fp, fine_mesh_path = get_properties()
    set_permeability(fp, layer=layer)
    bc = set_boundary_conditions(fp)
    set_weights_nodes(fp, update=False)
    transm = set_fine_transmissibility_without_bc_v2(fp)
    T: sp.csc_matrix = transm['T_tpfa']
    
    TPFA_mask: sp.csr_matrix = T.copy()
    TPFA_mask.data[:] = 1
    TPFA_mask = TPFA_mask.astype(bool)
    
    
    T_complete = transm['transmissibility_without_bc']
    
    ## interfaces with n^t K n
    resp5 = mount_kmatrix(fp) 
    T_matrix = resp5['T_matrix']
    
    T2 = T.copy()
    T2.data[:] = 1.0
    
    T3: sp.csc_matrix = T2.multiply(T_complete)
    T3.setdiag(0)
    soma = np.array(T3.sum(axis=1)).flatten()
    T3.setdiag(-soma)    
    
    
    # ag = AlgorithimicMonotone()
    # T_for_OP = ag.get_monotone_matrix(T_complete)

    enhanced = Enhanced()
    T_for_OP = enhanced.get_enhanced_matrix(T_complete)
    
    resp = set_fine_transmissibility_v2(fp, bc)
    T_bc = resp['transmissibility']
    b_bc = resp['source']
    
    nparts = int(T.shape[0]/nCr)
    
    primal = cms.create_partition(T_matrix, nparts, disjointed=disjointed, nvols_mean=nCr)
    all_regions = cms.create_support_region_and_boundary_v3(T_matrix, primal, n_levels_adj=n_levels_adj, ext='_1')
    
    OR = cms.get_OR_finite_volume(primal)
    OP, op_iterations, emax = cms.get_msrsb_prolongation_operator_v3(
        all_regions['ind_support'],
        all_regions['ptr_support'],
        T_for_OP,
        OR.transpose(copy=True),
        **my_params
    )
    
    
    T_strong = cms.define_strong_coupled(T_matrix, **my_params)
    primal_strong = cms.create_partition(T_strong, nparts=nparts, disjointed=True, nvols_mean=nCr)
    all_regions_strong = cms.create_support_region_and_boundary_v3(T_strong, primal_strong, n_levels_adj=n_levels_adj, ext='_1')
    fp.insert_or_update_data({'primal_id_level1': primal_strong})
    
    ## plotar interfaces da primal
    # print_fine_interfaces_coarse_mesh_2d(
    #     fp,
    #     fine_mesh_path,
    #     1,
    #     'edges_primal_intarfaces'
    # )
    
    OR_strong = cms.get_OR_finite_volume(primal_strong)
    OP_strong = OR_strong.transpose(copy=True)
    Msupport = sp.csr_matrix((np.repeat(1.0, all_regions_strong['ind_support'].shape[0]), all_regions_strong['ind_support'], all_regions_strong['ptr_support']), shape=OP_strong.transpose().shape).transpose()

    OR_strong = cms.get_OR_finite_volume(primal_strong)
    OP_strong, op_iterations, emax = cms.get_msrsb_prolongation_operator_v3(
        all_regions_strong['ind_support'],
        all_regions_strong['ptr_support'],
        T_for_OP,
        OR_strong.transpose(copy=True),
        **my_params
    )
    
    
    M2_strong = cms.define_strong_coupled_v2(T_complete, **my_params)
    M2_strong_mask = M2_strong.copy()
    M2_strong_mask.data[:] = 1
    M2_strong_mask = TPFA_mask.multiply(M2_strong_mask)
    
    primal_strong_M2 = cms.create_partition(M2_strong_mask, nparts=nparts, disjointed=True, nvols_mean=nCr)
    all_regions_strong_M2 = cms.create_support_region_and_boundary_v3(M2_strong_mask, primal_strong_M2, n_levels_adj=n_levels_adj, ext='_1')
    
    OR_strong_M2 = cms.get_OR_finite_volume(primal_strong_M2)
    OP_strong_M2 = OR_strong_M2.transpose(copy=True)
    Msupport_M2 = sp.csr_matrix((np.repeat(1.0, all_regions_strong_M2['ind_support'].shape[0]), all_regions_strong_M2['ind_support'], all_regions_strong_M2['ptr_support']), shape=OP_strong_M2.transpose().shape).transpose()
    
    S_M2 = cms.get_S_matrix(M2_strong)
    
    OP_strong_M2, op_iterations, emax = cms.get_msrsb_prolongation_operator_v3(
        all_regions_strong_M2['ind_support'],
        all_regions_strong_M2['ptr_support'],
        S_M2,
        OR_strong_M2.transpose(copy=True),
        **my_params
    )
    
    M3_strong = cms.define_strong_coupled_v3(T_complete, **my_params)
    M3_strong_mask = M3_strong.copy()
    M3_strong_mask.data[:] = 1
    M3_strong_mask = TPFA_mask.multiply(M3_strong_mask)
    
    primal_strong_M3 = cms.create_partition(M3_strong_mask, nparts=nparts, disjointed=True, nvols_mean=nCr)
    all_regions_strong_M3 = cms.create_support_region_and_boundary_v3(M3_strong_mask, primal_strong_M3, n_levels_adj=n_levels_adj, ext='_1')
    
    OR_strong_M3 = cms.get_OR_finite_volume(primal_strong_M3)
    OP_strong_M3 = OR_strong_M3.transpose(copy=True)
    Msupport_M3 = sp.csr_matrix((np.repeat(1.0, all_regions_strong_M3['ind_support'].shape[0]), all_regions_strong_M3['ind_support'], all_regions_strong_M3['ptr_support']), shape=OP_strong_M3.transpose().shape).transpose()
    
    S_M3 = cms.get_S_matrix(M3_strong)
    
    OP_strong_M3, op_iterations, emax = cms.get_msrsb_prolongation_operator_v3(
        all_regions_strong_M3['ind_support'],
        all_regions_strong_M3['ptr_support'],
        S_M3,
        OR_strong_M3.transpose(copy=True),
        **my_params
    )
    
    M4_strong = cms.define_strong_coupled_v2(T_for_OP, **my_params)
    M4_strong_mask = M4_strong.copy()
    M4_strong_mask.data[:] = 1
    M4_strong_mask = TPFA_mask.multiply(M4_strong_mask)
    
    primal_strong_M4 = cms.create_partition(M4_strong_mask, nparts=nparts, disjointed=True, nvols_mean=nCr)
    all_regions_strong_M4 = cms.create_support_region_and_boundary_v3(M4_strong_mask, primal_strong_M4, n_levels_adj=n_levels_adj, ext='_1')
    
    OR_strong_M4 = cms.get_OR_finite_volume(primal_strong_M4)
    OP_strong_M4 = OR_strong_M4.transpose(copy=True)
    Msupport_M4 = sp.csr_matrix((np.repeat(1.0, all_regions_strong_M4['ind_support'].shape[0]), all_regions_strong_M4['ind_support'], all_regions_strong_M4['ptr_support']), shape=OP_strong_M4.transpose().shape).transpose()
    
    S_M4 = cms.get_S_matrix(M4_strong)
    
    OP_strong_M4, op_iterations, emax = cms.get_msrsb_prolongation_operator_v3(
        all_regions_strong_M4['ind_support'],
        all_regions_strong_M4['ptr_support'],
        S_M4,
        OR_strong_M4.transpose(copy=True),
        **my_params
    )
    
    M5_strong = cms.compute_strength_matrix_symetric(T_complete, **my_params)
    M5_strong_mask = M5_strong.copy()
    M5_strong_mask.data[:] = 1
    M5_strong_mask = TPFA_mask.multiply(M5_strong_mask)
    
    primal_strong_M5 = cms.create_partition(M5_strong_mask, nparts=nparts, disjointed=True, nvols_mean=nCr)
    all_regions_strong_M5 = cms.create_support_region_and_boundary_v3(M5_strong_mask, primal_strong_M5, n_levels_adj=n_levels_adj, ext='_1')
    
    OR_strong_M5 = cms.get_OR_finite_volume(primal_strong_M5)
    OP_strong_M5 = OR_strong_M5.transpose(copy=True)
    Msupport_M5 = sp.csr_matrix((np.repeat(1.0, all_regions_strong_M5['ind_support'].shape[0]), all_regions_strong_M5['ind_support'], all_regions_strong_M5['ptr_support']), shape=OP_strong_M5.transpose().shape).transpose()
    
    S_M5 = cms.get_S_matrix(M5_strong)
    
    OP_strong_M5, op_iterations, emax = cms.get_msrsb_prolongation_operator_v3(
        all_regions_strong_M5['ind_support'],
        all_regions_strong_M5['ptr_support'],
        S_M5,
        OR_strong_M5.transpose(copy=True),
        **my_params
    )
    
    
    
    norm_b_bc = np.linalg.norm(b_bc)

    definitions = {
        'A': T_bc,
        'b': b_bc,
        'restart': 10,
        'tol': 1e-8,
        'x0': np.zeros_like(b_bc),
        'maxiter': int(50)
    }
    
    xref = spsolve(T_bc, b_bc)   
        

    finescale_smoother = mcl.Ilu0Precond(definitions['A'])
    # finescale_smoother_ilu1 = mcl.Ilu1Precond(definitions['A'])
    
    residuals_fs = []
    definitions.update({
        'M': finescale_smoother,
        'residuals': residuals_fs
    })
    method = 'Finescale Ilu0'
    
    def callback_func(xk):
        residuo = definitions.get('residuals')
        it = len(residuo)
        res = residuo[-1]
        print(f"Method: {method}, Iteracao: {it}, Residuo: {res:.2f}")
    
    t1 = time.perf_counter()
    x3, info = fgmres(**definitions, callback=callback_func)
    dt3 = time.perf_counter() - t1
    
    

    residual_op = []    
    precond = mcl.MultiScaleSmoother(T_bc, OP, OP.transpose(copy=True), finescale_smoother)
    definitions.update({
        'M': precond,
        'residuals': residual_op
    })
    method = 'MsRSB Method'
    t1 = time.perf_counter()
    x1, info = fgmres(**definitions, callback=callback_func)
    dt1 = time.perf_counter() - t1
    
    

    residual_op_strong = []
    precond = mcl.MultiScaleSmoother(T_bc, OP_strong, OP_strong.transpose(copy=True), finescale_smoother)
    definitions.update({
        'M': precond,
        'residuals': residual_op_strong
    })
    method = 'f-MsRSB Method'
    t1 = time.perf_counter()
    x2, info = fgmres(**definitions, callback=callback_func)
    dt2 = time.perf_counter() - t1
    
    residual_op_M2 = []    
    precond = mcl.MultiScaleSmoother(T_bc, OP_strong_M2, OP_strong_M2.transpose(copy=True), finescale_smoother)
    definitions.update({
        'M': precond,
        'residuals': residual_op_M2
    })
    method = 'MsRSB Method_M2'
    t1 = time.perf_counter()
    xM2, info = fgmres(**definitions, callback=callback_func)
    dtM2 = time.perf_counter() - t1
    
    residual_op_M3 = []    
    precond = mcl.MultiScaleSmoother(T_bc, OP_strong_M3, OP_strong_M3.transpose(copy=True), finescale_smoother)
    definitions.update({
        'M': precond,
        'residuals': residual_op_M3
    })
    method = 'MsRSB Method_M3'
    t1 = time.perf_counter()
    xM3, info = fgmres(**definitions, callback=callback_func)
    dtM3 = time.perf_counter() - t1
    
    residual_op_M4 = []    
    precond = mcl.MultiScaleSmoother(T_bc, OP_strong_M4, OP_strong_M4.transpose(copy=True), finescale_smoother)
    definitions.update({
        'M': precond,
        'residuals': residual_op_M4
    })
    method = 'MsRSB Method_M4'
    t1 = time.perf_counter()
    xM4, info = fgmres(**definitions, callback=callback_func)
    dtM4 = time.perf_counter() - t1
    
    residual_op_M5 = []    
    precond = mcl.MultiScaleSmoother(T_bc, OP_strong_M5, OP_strong_M5.transpose(copy=True), finescale_smoother)
    definitions.update({
        'M': precond,
        'residuals': residual_op_M5
    })
    method = 'MsRSB Method_M5'
    t1 = time.perf_counter()
    xM5, info = fgmres(**definitions, callback=callback_func)
    dtM5 = time.perf_counter() - t1
    
    print('#'*30)
    print(f"""
          OP1: MsRSB tradicional com o metodo enhanced
          OP_strong: f-MsRSB tradicional usando como filtro o n^t K n
          OP_M2: f-MsRSB tradicional
          OP_M3: usando o peso do f-MsRSB calculado em cada face
          OP_M4: Enhanced + f-MsRSB tradicional
          OP_M5: f-MsRSB usando o acumulado dos pesos (teste)
    """)
    print('#'*30)

    print('#'*30)
    print(f"""
        Comparacao precondicionador:
        Tempo:
            OP: {dt1:.2f}
            OP_strong: {dt2:.2f}
            OP_M2: {dtM2:.2f}
            OP_M3: {dtM3:.2f}
            OP_M4: {dtM4:.2f}
            OP_M5: {dtM5:.2f}

        Iteracoes:
            OP: {len(residual_op)}
            OP_strong: {len(residual_op_strong)}
            OP_M2: {len(residual_op_M2)}
            OP_M3: {len(residual_op_M3)}
            OP_M4: {len(residual_op_M4)}
            OP_M5: {len(residual_op_M5)}
            
        Residuo:
            OP: {residual_op[-1]:.2f}
            OP_strong: {residual_op_strong[-1]:.2f}
            OP_M2: {residual_op_M2[-1]:.2f}
            OP_M3: {residual_op_M3[-1]:.2f}
            OP_M4: {residual_op_M4[-1]:.2f}
            OP_M5: {residual_op_M5[-1]:.2f}
    """)
    print('#'*30)
    
    xapp1, err1, err1l2, res1 = get_app_solution(T_bc, b_bc, OP, OP.transpose(copy=True), xref)
    xapp2, err2, err2l2, res2 = get_app_solution(T_bc, b_bc, OP_strong, OP_strong.transpose(copy=True), xref)
    xapp3, err3, err3l2, res3 = get_app_solution(T_bc, b_bc, OP_strong_M2, OP_strong_M2.transpose(copy=True), xref)
    xapp4, err4, err4l2, res4 = get_app_solution(T_bc, b_bc, OP_strong_M3, OP_strong_M3.transpose(copy=True), xref)
    xapp5, err5, err5l2, res5 = get_app_solution(T_bc, b_bc, OP_strong_M4, OP_strong_M4.transpose(copy=True), xref)
    xapp6, err6, err6l2, res6 = get_app_solution(T_bc, b_bc, OP_strong_M5, OP_strong_M5.transpose(copy=True), xref)
    
    print('#'*30)
    print(f"""
            Comparacao solucao_aproximada:
            Inf norm:
                OP: {np.abs(err1).max():.2f}
                OP_strong: {np.abs(err2).max():.2f}
                OP_M2: {np.abs(err3).max():.2f}
                OP_M3: {np.abs(err4).max():.2f}
                OP_M4: {np.abs(err5).max():.2f}
                OP_M5: {np.abs(err6).max():.2f}
    
            2 norm:
                OP: {err1l2:.2f}
                OP_strong: {err2l2:.2f}
                OP_M2: {err3l2:.2f}
                OP_M3: {err4l2:.2f}
                OP_M4: {err5l2:.2f}
                OP_M5: {err6l2:.2f}
                
            Residuo:
                OP: {res1:.2f}
                OP_strong: {res2:.2f}
                OP_M2: {res3:.2f}
                OP_M3: {res4:.2f}
                OP_M4: {res5:.2f}
                OP_M5: {res6:.2f}
        """)
    print('#'*30)
    
    

    import ipdb; ipdb.set_trace()




    
    n = T_strong.shape[0]
    diag1 = T_strong.diagonal()
    D1 = sp.spdiags(1/diag1, 0, n, n).tocsc()
    
    diag2 = T_for_OP.diagonal()
    D2 = sp.spdiags(1/diag2, 0, n, n).tocsc()
    
    alpha = 0
    omega = 2/3
    for i in range(8):
        if i > 5:
            alpha = 0.7
            
        new_OP = OP_strong - omega*((1-alpha)*(D1*T_strong) + alpha*(D2*T_for_OP))*OP_strong
        new_OP = Msupport.multiply(new_OP)
        soma = np.array(new_OP.sum(axis=1)).flatten()
        soma[:] = 1/soma
        new_OP.data *= soma[new_OP.indices]
        emax = np.absolute((OP_strong - new_OP).data).max()
        print(i)
        print(emax)
        OP_strong = new_OP
    
    xf = spsolve(T_bc, b_bc)
    
    LU1 = factorized(OR @ T_bc @ OP)
    x1_app = OP @ LU1(OR @ b_bc)
    
    LU2 = factorized(OR_strong @ T_bc @ OP_strong)
    x2_app = OP_strong @ LU2(OR_strong @ b_bc)
    
    error1 = np.abs(xf - x1_app)
    error2 = np.abs(xf - x2_app)
    
    l2_error1 = np.linalg.norm(error1)
    l2_error2 = np.linalg.norm(error2)
    
    
    
    import pdb; pdb.set_trace()
    
    
    
    
    from packs.manager.mesh_data import MeshData
    mesh_data = MeshData(dim=3, mesh_path=fine_mesh_path)
    mesh_data.create_tag('permx')
    mesh_data.insert_tag_data('permx', fp['permeability'][:,0,0], elements_type='faces')
    
    filtered_interfaces = identify_filtered_interfaces(T_strong, fp['adjacencies'])
    
    mesh_data.export_all_elements_type_to_vtk('spe_perms', element_type='faces')
    mesh_data.export_only_the_elements('filtered_interfaces', element_type='edges', elements_array=filtered_interfaces)
    
    
    
    
    
    
    # save_data(T_bc, b_bc, OP, OR, layer, nCr)
    
    
    
    
    
    
    
    
    
    
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
    
def get_app_solution(A, b, OP, OR, xref):
    Ac = OR @ A @ OP
    sol = OP @ spsolve(Ac, OR @ b)
    err = xref - sol
    errl2 = np.linalg.norm(err)
    res = np.linalg.norm(b - A @ sol)
    return sol, err, errl2, res

def run5(layer=36, nCr=81):
    
    disjointed = False
    n_levels_adj = 3
    my_params = get_params()

    fp, fine_mesh_path = get_properties()
    set_permeability(fp, layer=layer)
    bc = set_boundary_conditions(fp)
    set_weights_nodes(fp, update=False)
    transm = set_fine_transmissibility_without_bc_v2(fp)
    T: sp.csc_matrix = transm['T_tpfa']

    T_complete = transm['transmissibility_without_bc']
    
    
    ag = AlgorithimicMonotone()
    T_for_OP = ag.get_monotone_matrix(T_complete)
    
    resp = set_fine_transmissibility_v2(fp, bc)
    T_bc = resp['transmissibility']
    b_bc = resp['source']
    
    nparts = int(T.shape[0]/nCr)
    
    ## primal criada com adjacencia tpfa
    primal = cms.create_partition(T, nparts)
    all_regions = cms.create_support_region_and_boundary_v3(T_complete, primal, n_levels_adj=n_levels_adj, ext='_1')
    fp.insert_or_update_data({'primal_id_level1': primal})
    
    OR = cms.get_OR_finite_volume(primal)
    OP, op_iterations, emax = cms.get_msrsb_prolongation_operator_v3(
        all_regions['ind_support'],
        all_regions['ptr_support'],
        T_for_OP,
        OR.transpose(copy=True),
        **my_params
    )
    
    
    # # plotar interfaces da primal
    # print_fine_interfaces_coarse_mesh_2d(
    #     fp,
    #     fine_mesh_path,
    #     1,
    #     'edges_primal_MsRSB'
    # )
    
    
    
    T_strong = cms.define_strong_coupled(T_complete, **my_params)
    primal_strong = cms.create_partition(T_strong, nparts=nparts, disjointed=True, nvols_mean=nCr)
    all_regions_strong = cms.create_support_region_and_boundary_v3(T_strong, primal_strong, n_levels_adj=n_levels_adj, ext='_1')
    fp.insert_or_update_data({'primal_id_level1': primal_strong})
    
    # # plotar interfaces da primal
    # print_fine_interfaces_coarse_mesh_2d(
    #     fp,
    #     fine_mesh_path,
    #     1,
    #     'edges_primal_f-MsRSB'
    # )
    
    
    OR_strong = cms.get_OR_finite_volume(primal_strong)
    OP_s, op_iterations_s, emax_s = cms.get_msrsb_prolongation_operator_v3(
        all_regions_strong['ind_support'],
        all_regions_strong['ptr_support'],
        T_strong,
        OR_strong.transpose(copy=True),
        **my_params
    )
    
    xf = spsolve(T_bc, b_bc)
    
    M_MsRSB = mcl.MultiScaleIlu0Smoother(T_bc, OP, OP.transpose(copy=True))
    M_fMsRSB = mcl.MultiScaleIlu0Smoother(T_bc, OP_s, OP_s.transpose(copy=True))
    
    r_MsRSB = []
    r_fMsRSB = []
    
    def meu_callback_Ms(rk):
        itM = len(r_MsRSB)
        msg = f"Norma do residuo: {r_MsRSB[-1]:.2e}, Iteracao: {itM}"
        print(msg)
        
    def meu_callback_fMs(rk):
        itfM = len(r_fMsRSB)
        msg = f"Norma do residuo: {r_fMsRSB[-1]:.2e}, Iteracao: {itfM}"
        print(msg)
    
    
    
    t1 = time.perf_counter()
    x_fMsRSB, exitcode = fgmres(T_bc, b_bc, residuals=r_fMsRSB, M=M_fMsRSB, maxiter=my_params['maxiter'], restart=my_params['restart'], tol=my_params['tol'], callback=meu_callback_fMs)
    dt_fMsRSB = time.perf_counter() - t1
    
    t1 = time.perf_counter()
    x_MsRSB, exitcode = fgmres(T_bc, b_bc, residuals=r_MsRSB, M=M_MsRSB, maxiter=my_params['maxiter'], restart=my_params['restart'], tol=my_params['tol'], callback=meu_callback_Ms)
    dt_MsRSB = time.perf_counter() - t1
    
    
    print(f'MsRSB: {dt_MsRSB}')
    print(f'f-MsRSB: {dt_fMsRSB}')
    
    print(f'error MsRSB: {np.linalg.norm(xf - x_MsRSB)}')
    print(f'error f-MsRSB: {np.linalg.norm(xf - x_fMsRSB)}')
    
    print(f'Len res MsRSB: {len(r_MsRSB)}')
    print(f'Len res fMsRSB: {len(r_fMsRSB)}')
    
    
    
    
    
    
    
    
    
    import pdb; pdb.set_trace()
    
    
    
    
    from packs.manager.mesh_data import MeshData
    mesh_data = MeshData(dim=3, mesh_path=fine_mesh_path)
    mesh_data.create_tag('permx')
    mesh_data.insert_tag_data('permx', fp['permeability'][:,0,0], elements_type='faces')
    
    filtered_interfaces = identify_filtered_interfaces(T_strong, fp['adjacencies'])
    
    mesh_data.export_all_elements_type_to_vtk('spe_perms', element_type='faces')
    mesh_data.export_only_the_elements('filtered_interfaces', element_type='edges', elements_array=filtered_interfaces)   
    
    
# teste1: malha 15x15 com a barreira nas 5 camadas do meio
# teste2: spe10 camada 16
# objetivo: filtrar as mesmas faces do f-MsRSB
    
    
    
    
    
    