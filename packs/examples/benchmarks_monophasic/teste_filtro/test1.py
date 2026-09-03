from packs.mpfa_methods.mesh_preprocess import preprocess_mesh
from packs import defpaths
from packs.manager import BoundaryConditions
from packs.examples.benchmarks_monophasic.cross.test_cross_2 import set_weights_nodes, set_fine_transmissibility_v2, set_fine_transmissibility_without_bc_v2
from packs.manager.meshmanager2 import MeshProperty
from packs.multiscale.msrsb import create_msrsb_structure as cms
from packs.examples.benchmarks_monophasic.spe10.spe1 import identify_filtered_interfaces
from packs.manager.meshmanager2 import MeshProperty

import os
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from pathlib import Path
import gmsh
import copy



def get_R(theta):
    
    R = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])
    
    return R
 
    
def set_permeability(fp:MeshProperty):
    
    k1 = np.array([
        [1000, 0],
        [0, 1]
    ])
    
    k2 = np.array([
            [100, 0],
            [0, 10]
    ])*1
    
    theta = np.deg2rad(45)
    R = get_R(theta)
    
    k_rot = R @ k1 @ R.T
    # k_rot = k1*100
    # k2_rot = R @ k2 @ R.T
    
    nfaces = len(fp.faces)
    perm = np.zeros((nfaces, 2, 2))
    perm[:] = k1
    
    xcentroids = fp.faces_centroids[:,0]
    
    test1 = xcentroids>3
    
    perm[test1] = k_rot
    # perm[test1] = k2
    
    
    
    
    
    fp.insert_or_update_data({'permeability': perm})
    fp.export_data()
    

def set_permeability_bar(fp:MeshProperty):
    
    k1 = np.array([
        [1e5, 0],
        [0, 1e5]
    ])
    
    k2 = np.array([
            [1, 0],
            [0, 1]
    ])
    
    
    nfaces = len(fp.faces)
    perm = np.zeros((nfaces, 2, 2))
    perm[:] = k1
    
    xcentroids = fp.faces_centroids[:,0]
    
    test1 = xcentroids<6
    test2 = xcentroids>4
    test3 = test1 & test2
    
    perm[test3] = k2
    
    fp.insert_or_update_data({'permeability': perm})
    fp.export_data()
    
    



def set_boundary_conditions(fine_properties: MeshProperty):
    bc = BoundaryConditions()

    nodes_centroids = fine_properties['nodes_centroids']
    
    
    faces = fine_properties['faces']
    faces_centroids = fine_properties['faces_centroids']
    xcentroids = faces_centroids[:,0]
    
    xmin, ymin = nodes_centroids.min(axis=0)
    xmax, ymax = nodes_centroids.max(axis=0)

    # c_p1 = np.array([xmin, ymin])
    # c_p0 = np.array([xmax, ymax])

    c_p1 = np.array([xmin, ymax])
    c_p0 = np.array([xmax, ymin])
    
    t1 = xcentroids < 1
    t2 = xcentroids > 9
    
    dists = np.linalg.norm(faces_centroids - c_p1, axis=1)
    # face_p1 = faces[dists <= dists.min()][0]
    face_p1 = faces[t1]
    pressure_p1 = np.repeat(1e5, t1.sum())
    # dists[:] = np.linalg.norm(faces_centroids - c_p0, axis=1)
    # face_p0 = faces[dists <= dists.min()][0]
    face_p0 = faces[t2]
    pressure_p0 = np.repeat(1, t2.sum())
    
    faces_pressure = np.concatenate([face_p0, face_p1])
    pressure_presc = np.concatenate([pressure_p0, pressure_p1])

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

    # bc.set_boundary('injectors', np.array([face_p1]), np.array([True]))
    # bc.set_boundary('producers', np.array([face_p0]), np.array([True]))
    
    bc.set_boundary('injectors', face_p1, np.repeat(True, t1.sum()))
    bc.set_boundary('producers', face_p0, np.repeat(True, t2.sum()))

    bc.update_zero_bcs()

    return bc

def set_boundary_conditions_pert(fine_properties: MeshProperty):
    bc = BoundaryConditions()

    nodes_centroids = fine_properties['nodes_centroids']
    
    
    faces = fine_properties['faces']
    faces_centroids = fine_properties['faces_centroids']
    xcentroids = faces_centroids[:,0]
    
    xmin, ymin = nodes_centroids.min(axis=0)
    xmax, ymax = nodes_centroids.max(axis=0)

    # c_p1 = np.array([xmin, ymin])
    # c_p0 = np.array([xmax, ymax])

    c_p1 = np.array([xmin, ymax])
    c_p0 = np.array([xmax, ymin])
    
    t1 = xcentroids < 1
    t2 = xcentroids > 9
    
    dists = np.linalg.norm(faces_centroids - c_p1, axis=1)
    # face_p1 = faces[dists <= dists.min()][0]
    face_p1 = faces[t1]
    pressure_p1 = np.repeat(1e5, t1.sum())
    # dists[:] = np.linalg.norm(faces_centroids - c_p0, axis=1)
    # face_p0 = faces[dists <= dists.min()][0]
    face_p0 = faces[t2]
    pressure_p0 = np.repeat(1, t2.sum())
    
    faces_pressure = np.concatenate([face_p0, face_p1])
    pressure_presc = np.concatenate([pressure_p0, pressure_p1])

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

    # bc.set_boundary('injectors', np.array([face_p1]), np.array([True]))
    # bc.set_boundary('producers', np.array([face_p0]), np.array([True]))
    
    bc.set_boundary('injectors', face_p1, np.repeat(True, t1.sum()))
    bc.set_boundary('producers', face_p0, np.repeat(True, t2.sum()))

    bc.update_zero_bcs()

    return bc


def get_properties():
    
    fine_mesh_path = os.path.join('mesh', 'quadrado_estruturado10x10.msh')
    fine_mesh_properties_name = 'quadrado_estruturado10x10' 
    fine_mesh_path_v4 = fine_mesh_path

    fine_properties = preprocess_mesh(fine_mesh_path, fine_mesh_properties_name, mesh_name_v4=fine_mesh_path_v4)

    return fine_properties, fine_mesh_path

def plot_filtered_interfaces(mesh_path, fp: MeshProperty, T_strong: sp.csr_matrix, ext=''):
    from packs.manager.mesh_data import MeshData
    mesh_data = MeshData(dim=2, mesh_path=mesh_path)
    mesh_data.create_tag('permx')
    mesh_data.insert_tag_data('permx', fp['permeability'][:,0,0], elements_type='faces')
    mesh_data.create_tag('fid', data_type='int')
    mesh_data.insert_tag_data('fid', fp.faces, elements_type='faces')
    
    filtered_interfaces = identify_filtered_interfaces(T_strong, fp['adjacencies'])
    filtered_interfaces = np.setdiff1d(filtered_interfaces, fp.boundary_edges)
    
    mesh_data.export_all_elements_type_to_vtk('spe_perms', element_type='faces')
    mesh_data.export_only_the_elements('filtered_interfaces'+ext, element_type='edges', elements_array=filtered_interfaces)   
    


def test_edge(x, y, fp: MeshProperty, T_withoutbc, debug=False):
    edges_centroids = fp.edges_centroids
    selected_edge = fp.edges[
        ((np.abs(edges_centroids[:, 0] - x)) < 1e-10) &
        ((np.abs(edges_centroids[:, 1] - y)) < 1e-10) 
    ]
    
    nodes_selected_edge = fp.nodes_of_edges[selected_edge[0]]
    faces_of_nodes_edge = fp.faces_of_nodes[nodes_selected_edge]
    kl = fp.adjacencies[selected_edge[0]]
    
    nodes_weight = fp.nodes_weights
    xi_params = fp.xi_params
    
    xi_params_edge = xi_params[selected_edge]
    
    node_weight_A = nodes_weight[nodes_weight['node_id'] == nodes_selected_edge[0]]
    node_weight_B = nodes_weight[nodes_weight['node_id'] == nodes_selected_edge[1]]
    
    fA = node_weight_A['face_id']
    fB = node_weight_B['face_id']
    
    transm_K = T_withoutbc[kl[0]]
    transm_L = T_withoutbc[kl[1]]
    
    if debug == True:
        print(transm_K)
        print('#'*30)
        print(transm_L)
        print('#'*30)
        print(node_weight_A)
        print('#'*30)
        print(node_weight_B)
        print('#'*30)
        print(fA)
        print('#'*30)
        print(fB)
        print('#'*30)
        print(kl)
        import ipdb; ipdb.set_trace()

def get_params():
    return {
        'eps': 1e-3,
        'theta': 0.25,
        'gamma': 0.9
    }

def perturb_internal_nodes(points, boundary_nodes, scale_factor=0.0005):
    
    """
    Perturba os nós internos usando uma distribuição gaussiana.
    scale_factor: fração do menor lado do bounding box usada como desvio padrão.
    """
    # 1. Identificar fronteira e nós internos
    
    all_nodes = np.arange(len(points))
    internal_nodes = np.setdiff1d(all_nodes, boundary_nodes)
    
    if len(internal_nodes) == 0:
        print("Nenhum nó interno encontrado para perturbar.")
        return points
        
    # 2. Calcular escala da perturbação baseada no tamanho da malha
    bbox = np.ptp(points, axis=0) # Diferença entre max e min de cada eixo
    scale = np.min(bbox) * scale_factor 
    
    # 3. Gerar perturbação gaussiana
    perturbation = np.random.normal(0, scale, size=(len(internal_nodes), 2))
    
    # 4. Limitar a perturbação máxima (clip) para evitar inversão de elementos
    max_perturbation = scale * 2.5
    perturbation = np.clip(perturbation, -max_perturbation, max_perturbation)
    
    # 5. Aplicar perturbação
    points_perturbed = points.copy()
    
    # points_perturbed[internal_nodes] += perturbation
    
    dx = np.random.uniform(low=-scale_factor, high=scale_factor, size=internal_nodes.shape[0])
    dy = np.random.uniform(low=-scale_factor, high=scale_factor, size=internal_nodes.shape[0])
    points_perturbed[internal_nodes,0] += dx
    points_perturbed[internal_nodes,1] += dy
    
    print(f"Perturbados {len(internal_nodes)} nós internos. Escala usada: {scale:.4f}")
    return points_perturbed

def ordenar_antihorario(pontos):
    """
    Ordena pontos 2D em ordem anti-horária em torno do centroide.
    
    Parâmetros:
        pontos: array numpy de shape (N, 2) com as coordenadas (x, y)
    
    Retorna:
        pontos_ordenados: array (N, 2) ordenado anti-horariamente
        indices: array com os índices da ordenação original
    """
    pontos = np.asarray(pontos)
    
    # 1. Calcular o centroide (média dos pontos)
    centroide = np.mean(pontos, axis=0)
    
    # 2. Calcular os ângulos polares de cada ponto em relação ao centroide
    # arctan2 retorna valores em [-pi, pi], perfeito para ordenação angular
    angulos = np.arctan2(pontos[:, 1] - centroide[1], 
                         pontos[:, 0] - centroide[0])
    
    # 3. Obter os índices que ordenam os ângulos (anti-horário = crescente)
    indices = np.argsort(angulos)[::-1]
    
    return pontos[indices], indices

def create_pertubed_mesh(fp: MeshProperty, new_mesh_name: str):
    from pymoab import core, types, rng, topo_util
    
    nodes_cetroids = fp.nodes_centroids
    boundary_nodes = fp.boundary_nodes
    z = np.repeat(0.0, len(nodes_cetroids))
    
    nodes_of_faces = fp['nodes_of_faces']
    
    new_nodes_of_faces = []
    for nodes in nodes_of_faces:
        new_nodes_of_faces.append(nodes[ordenar_antihorario(nodes_cetroids[nodes])[1]])
    new_nodes_of_faces = np.array(new_nodes_of_faces)
    
    perturbed_coords = perturb_internal_nodes(nodes_cetroids, boundary_nodes)
    perturbed_coords = np.column_stack([perturbed_coords, z])
    
    
    
    
    mb = core.Core()
    root_set = mb.get_root_set()
    mtu = topo_util.MeshTopoUtil(mb)
    verts = mb.create_vertices(perturbed_coords.flatten())
    
    quads = [verts[quad] for quad in new_nodes_of_faces]
    # quads = [verts[new_nodes_of_faces[0]]]
    
    # import ipdb; ipdb.set_trace()
    
    
    mb.create_elements(types.MBQUAD, quads)
    
    mname = new_mesh_name+'.vtk'
    export_mesh_name = str(Path(defpaths.results) / mname)
    
    mb.write_file(export_mesh_name)
    # mb.write_file(os.path.join('mesh', 'test_mesh.vtk'))
    

def perturbar_nos_internos(arquivo_entrada, arquivo_saida, fator_perturbacao=0.015, seed=42):
    gmsh.initialize()
    gmsh.model.add("malha_perturbada")
    gmsh.open(arquivo_entrada)
    
    # Obter todos os nós
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_tags = np.array(node_tags, dtype=int)
    node_coords = np.array(node_coords).reshape(-1, 3)
    
    print(f"Total de nós: {len(node_tags)}")
    
    # Identificar nós de fronteira (em 2D, dimensão 1)
    boundary_nodes = set()
    for dim, tag in gmsh.model.getEntities(1):
        line_nodes, _, _ = gmsh.model.mesh.getNodes(dim, tag)
        boundary_nodes.update(line_nodes)
    
    # Identificar nós internos
    all_nodes = set(node_tags)
    internal_nodes = np.array(list(all_nodes - boundary_nodes), dtype=int)
    
    print(f"Nós de fronteira: {len(boundary_nodes)}")
    print(f"Nós internos: {len(internal_nodes)}")
    
    # Calcular escala da perturbação
    bbox = np.ptp(node_coords, axis=0)
    tamanho_caracteristico = np.min(bbox[:2])
    escala = tamanho_caracteristico * fator_perturbacao
    
    # Gerar perturbação
    np.random.seed(seed)
    perturbacao = np.random.normal(0, escala, size=(len(internal_nodes), 3))
    max_perturb = escala * 2.5
    perturbacao = np.clip(perturbacao, -max_perturb, max_perturb)
    
    # Aplicar perturbação
    node_coords_perturbadas = node_coords.copy()
    tag_to_idx = {tag: i for i, tag in enumerate(node_tags)}
    
    for i, tag in enumerate(internal_nodes):
        idx = tag_to_idx[tag]
        node_coords_perturbadas[idx] += perturbacao[i]
        node_coords_perturbadas[idx, 2] = 0.0  # Manter z=0 para 2D
    
    # ==========================================
    # CORREÇÃO: Usar setNode (singular) em loop
    # ==========================================
    print("Atualizando coordenadas dos nós...")
    for i, tag in enumerate(node_tags):
        coord = node_coords_perturbadas[i]
        # setNode(tag, coord, parametricCoord)
        gmsh.model.mesh.setNode(int(tag), coord, [])
    
    # Salvar a malha
    gmsh.write(arquivo_saida)
    print(f"Malha salva em: {arquivo_saida}")
    
    gmsh.finalize()
    


def identificar_nos_fronteira_topologica(node_tags, cells_dict):
    """
    Identifica nós de fronteira pela topologia da malha (robusto para qualquer malha 2D).
    Uma aresta é de fronteira se aparece em apenas 1 elemento.
    """
    all_edges = []
    
    for cell_type, connectivity in cells_dict.items():
        # Extrair arestas baseado no tipo de elemento
        if cell_type == "triangle":
            # 3 arestas por triângulo
            edges = connectivity[:, [[0, 1], [1, 2], [2, 0]]]
        elif cell_type == "quad":
            # 4 arestas por quadrilátero
            edges = connectivity[:, [[0, 1], [1, 2], [2, 3], [3, 0]]]
        else:
            continue
        
        # Ordenar os nós de cada aresta (para que [1,2] e [2,1] sejam iguais)
        edges_sorted = np.sort(edges, axis=2)
        edges_reshaped = edges_sorted.reshape(-1, 2)
        all_edges.append(edges_reshaped)
    
    if not all_edges:
        return set()
    
    # Empilhar todas as arestas
    all_edges = np.vstack(all_edges)
    
    # Contar ocorrências de cada aresta
    unique_edges, counts = np.unique(all_edges, axis=0, return_counts=True)
    
    # Arestas que aparecem apenas 1 vez são de fronteira
    boundary_edges = unique_edges[counts == 1]
    
    # Nós de fronteira são todos os nós que compõem essas arestas
    boundary_nodes = set(np.unique(boundary_edges))
    
    print(f"Arestas de fronteira encontradas: {len(boundary_edges)}")
    print(f"Nós de fronteira identificados: {len(boundary_nodes)}")
    
    return boundary_nodes


def perturbar_nos_internos_2d(arquivo_entrada, arquivo_saida, fator_perturbacao=0.013, seed=42):
    """
    Lê uma malha 2D, perturba APENAS nós internos (fronteira preservada) e salva.
    """
    gmsh.initialize()
    gmsh.model.add("malha_perturbada")
    gmsh.open(arquivo_entrada)
    
    # ==========================================
    # 1. EXTRAIR NÓS E ELEMENTOS
    # ==========================================
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_tags = np.array(node_tags, dtype=int)
    node_coords = np.array(node_coords).reshape(-1, 3)
    
    # Mapeamento tag -> índice
    tag_to_idx = {int(tag): i for i, tag in enumerate(node_tags)}
    
    print(f"Total de nós: {len(node_tags)}")
    
    # ==========================================
    # 2. EXTRAIR CONECTIVIDADE DOS ELEMENTOS
    # ==========================================
    cells_dict = {}
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, -1)
    
    # Mapeamento de tipos gmsh para nomes
    type_map = {
        2: "triangle",    # Triângulo (3 nós)
        3: "quad"         # Quadrilátero (4 nós)
    }
    
    for etype, enodes in zip(elem_types, elem_node_tags):
        if etype in type_map:
            cell_type = type_map[etype]
            n_nodes_per_elem = 3 if etype == 2 else 4
            connectivity = np.array(enodes, dtype=int).reshape(-1, n_nodes_per_elem)
            
            # Converter tags do gmsh para índices 0-based
            conn_indices = np.array([[tag_to_idx[int(t)] for t in elem] 
                                     for elem in connectivity])
            
            if cell_type in cells_dict:
                cells_dict[cell_type] = np.vstack([cells_dict[cell_type], conn_indices])
            else:
                cells_dict[cell_type] = conn_indices
    
    print(f"Tipos de elementos encontrados: {list(cells_dict.keys())}")
    for ctype, conn in cells_dict.items():
        print(f"  {ctype}: {len(conn)} elementos")
    
    # ==========================================
    # 3. IDENTIFICAR NÓS DE FRONTEIRA (TOPOLOGIA)
    # ==========================================
    boundary_nodes = identificar_nos_fronteira_topologica(node_tags, cells_dict)
    
    # Identificar nós internos
    all_nodes = set(range(len(node_tags)))
    internal_nodes = np.array(list(all_nodes - boundary_nodes), dtype=int)
    
    print(f"Nós internos: {len(internal_nodes)}")
    
    if len(internal_nodes) == 0:
        print("Nenhum nó interno encontrado!")
        gmsh.finalize()
        return
    
    # ==========================================
    # 4. CALCULAR ESCALA E APLICAR PERTURBAÇÃO
    # ==========================================
    bbox = np.ptp(node_coords, axis=0)
    tamanho_caracteristico = np.min(bbox[:2])  # Apenas x e y
    escala = tamanho_caracteristico * fator_perturbacao
    
    print(f"Escala de perturbação: {escala:.6f}")
    
    np.random.seed(seed)
    perturbacao = np.random.normal(0, escala, size=(len(internal_nodes), 3))
    
    # Limitar perturbação máxima
    max_perturb = escala * 2.5
    perturbacao = np.clip(perturbacao, -max_perturb, max_perturb)
    
    # Aplicar perturbação APENAS nos nós internos
    node_coords_perturbadas = node_coords.copy()
    for i, idx in enumerate(internal_nodes):
        node_coords_perturbadas[idx] += perturbacao[i]
        node_coords_perturbadas[idx, 2] = 0.0  # Manter z=0 para 2D
    
    # ==========================================
    # 5. ATUALIZAR NÓS NO MODELO GMSH
    # ==========================================
    print("Atualizando coordenadas dos nós...")
    for i, tag in enumerate(node_tags):
        coord = node_coords_perturbadas[i]
        gmsh.model.mesh.setNode(int(tag), coord, [])
    
    # ==========================================
    # 6. SALVAR E FINALIZAR
    # ==========================================
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
    gmsh.write(arquivo_saida)
    print(f"\nMalha salva em: {arquivo_saida}")
    
    gmsh.finalize()







def run5():
    
    my_params = get_params()

    fp, fine_mesh_path = get_properties()
    set_permeability_bar(fp)
    bc = set_boundary_conditions(fp)
    set_weights_nodes(fp)
    resp = set_fine_transmissibility_v2(fp, bc)
    T_bc = resp['transmissibility']
    b_bc = resp['source']
    solution = spsolve(T_bc, b_bc)
    resp2 = set_fine_transmissibility_without_bc_v2(fp)
    T_complete = resp2['transmissibility_without_bc']
    
    T_strong1 = cms.define_strong_coupled(T_complete, **my_params)
    T_strong2 = cms.define_strong_coupled_v3(T_complete, **my_params)
    T_strong3 = cms.compute_strength_matrix_symetric(T_complete, **my_params)
    
    f1 = identify_filtered_interfaces(T_strong1, fp['adjacencies'])
    f2 = identify_filtered_interfaces(T_strong2, fp['adjacencies'])
    f3 = identify_filtered_interfaces(T_strong3, fp['adjacencies'])
    
    # plot_filtered_interfaces(fine_mesh_path, fp, T_strong1, ext='_default_strong')
    # plot_filtered_interfaces(fine_mesh_path, fp, T_strong2, ext='_max_row')
    # plot_filtered_interfaces(fine_mesh_path, fp, T_strong3, ext='_accum')
    
    # ## edge na interface
    # test_edge(3, 2.5, fp, T_complete, debug=True)
    # ## edge na regiao x > media
    # test_edge(4, 2.5, fp, T_complete, debug=True)
    # ## edge na regiao x < media
    # test_edge(2, 2.5, fp, T_complete, debug=True)
    
    
    
    
    
    
    new_mesh_name = Path(defpaths.mesh) / 'perturbed_mesh_10x10.msh'
    perturbar_nos_internos_2d(fine_mesh_path, str(new_mesh_name))
    mesh_properties_pert_name = 'quadrado_estruturado10x10_pert'
    fp_pert = preprocess_mesh(str(new_mesh_name), mesh_properties_pert_name)
    fp_pert.insert_or_update_data({
        'permeability': fp['permeability'].copy()
    })
    bc_pert = set_boundary_conditions_pert(fp_pert)
    set_weights_nodes(fp_pert)
    resp_pert = set_fine_transmissibility_v2(fp_pert, bc_pert)
    T_bc_pert = resp_pert['transmissibility']
    b_bc_pert = resp_pert['source']
    solution_pert = spsolve(T_bc_pert, b_bc_pert)
    resp2_pert = set_fine_transmissibility_without_bc_v2(fp_pert)
    T_complete_pert = resp2_pert['transmissibility_without_bc']
    
    T_strong1_pert = cms.define_strong_coupled(T_complete_pert, **my_params)
    T_strong2_pert = cms.define_strong_coupled_v3(T_complete_pert, **my_params)
    T_strong3_pert = cms.compute_strength_matrix_symetric(T_complete_pert, **my_params)
    
    f1_pert = identify_filtered_interfaces(T_strong1_pert, fp_pert['adjacencies'])
    f2_pert = identify_filtered_interfaces(T_strong2_pert, fp_pert['adjacencies'])
    f3_pert = identify_filtered_interfaces(T_strong3_pert, fp_pert['adjacencies'])
    
    
    

    
    
    import ipdb; ipdb.set_trace()
    
    
    
    
    