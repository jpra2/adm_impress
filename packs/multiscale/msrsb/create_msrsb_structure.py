import numpy as np
import time
import scipy.sparse as sp
import pandas as pd
import os
import networkx as nx
import metis
from rich import inspect, print as rprint
from typing import Tuple
import copy
from functools import reduce
import multiprocessing as mp
import time
from typing import Sequence
import inspect

from packs.manager import SuperArrayManager
from packs.utils import utils_old
from packs.utils.utils_old import time_func
from packs.fglobal import local_names, save_ijk_spe
from packs.multiscale.transmissibility_correction.enhanced import Enhanced
from packs.multiscale.transmissibility_correction.algorithimic_monotone import AlgorithimicMonotone
from packs.logs_func import create_operator_logger
from pyamg.strength import evolution_strength_of_connection

def is_sparse_symmetric(matrix: sp.csc_matrix, tol=1e-9):
    """
    Checks if a sparse SciPy matrix is symmetric.
    Returns:
        True if the matrix is symmetric, False otherwise.
    """
    if matrix.shape[0] != matrix.shape[1]:
        return False  # Not a square matrix, cannot be symmetric

    diff = matrix - matrix.transpose()
    return np.all(np.abs(diff.data) < tol)



def get_balanced_matrix(A: sp.csc_matrix) -> sp.csc_matrix:
    M = A + A.T
    v5 = M @ np.ones(M.shape[0])
    diagm = M.diagonal()
    v6 = diagm - v5
    M.setdiag(v6)
    M = 0.5*M
    return M


def get_dual_volumes(dual_id: np.ndarray, fine_faces_id: np.ndarray, fine_faces_of_faces: np.ndarray, ijks: np.ndarray, verify: np.ndarray):

    face_id = 0
    edge_id = 2
    vertice_id = 3

    dual_faces = fine_faces_id[dual_id == face_id]
    dual_others = fine_faces_id[dual_id != face_id]
    dual_volumes = []

    count_dual_volume = 0

    spe_data = load_spe_data()

    # save_ijk_spe(spe_data)

    # mesh_data = MeshData(dim=3, mesh_path='spe.h5m')
    # mesh_data.create_tag(tag_name='dual_id', data_type='int')
    dual_id_spe = np.repeat(-1, spe_data['ijk'].shape[0])
    mapvolumes = spe_data['map']
    test = mapvolumes > -1
    ids_spe = np.arange(spe_data['ijk'].shape[0])
    dual_id_spe[ids_spe[test]] = dual_id[mapvolumes[test]]
    # mesh_data.insert_tag_data('dual_id', dual_id, 'volumes')

    map_spe_to_fine_volumes(spe_data, ijks, verify)
    # plot_only_volumes_spe(fine_faces_id, spe_data, mesh_data, fine_faces_id)
    
    import pdb; pdb.set_trace()

    while dual_faces.shape[0] > 0:
        face0 = np.array([dual_faces[0]], dtype=np.int64)
        test = np.array([True], dtype=bool)
        
        while np.any(test):
            adjs = fine_faces_of_faces[np.isin(fine_faces_of_faces[:,0], face0)]
            faces_of_face0 = np.unique(np.concatenate(adjs))
            external_faces = np.setdiff1d(faces_of_face0, face0)
            test = dual_id[external_faces] == face_id
            face0 = np.setdiff1d(faces_of_face0, dual_others)

        adjs2 = fine_faces_of_faces[np.isin(fine_faces_of_faces[:,0], faces_of_face0)]
        faces_of_face0_v2 = np.unique(np.concatenate(adjs2))
        only_faces = np.setdiff1d(faces_of_face0_v2, faces_of_face0)

        adjs3 = np.unique(
            np.concatenate(
                fine_faces_of_faces[np.isin(fine_faces_of_faces[:,0], only_faces)]
            )
        )

        only_edges = adjs3[dual_id[adjs3] == edge_id]

        adjs4 = np.unique(
            np.concatenate(
                fine_faces_of_faces[np.isin(fine_faces_of_faces[:,0], only_edges)]
            )
        )

        only_vertices = adjs4[dual_id[adjs4] == vertice_id]
        dual_volume = np.concatenate([
            faces_of_face0,
            only_faces,
            only_edges,
            only_vertices
        ])

        # adjs3 = fine_faces_of_faces[np.isin(fine_faces_of_faces[:,0], faces_of_face0_v2)]


        # vertices_to_get = faces_of_face0_v2[
        #     dual_id[faces_of_face0_v2] == vertice_id
        # ]
        # faces_of_face0 = np.union1d(faces_of_face0, vertices_to_get)
        # # dual_volume = np.unique(np.concatenate(dual_volume))
        dual_faces = np.setdiff1d(dual_faces, faces_of_face0)
        # dual_volumes.append(faces_of_face0)
        dual_volumes.append(dual_volume)

        # mesh_data.export_only_the_elements('test_elements', 'volumes', dual_volume)

        import pdb; pdb.set_trace()
        
        
        unique_values, counts = np.unique(dual_id[dual_volume], return_counts=True)
        print(unique_values)
        print(counts)
        print(count_dual_volume)
        print()
        count_dual_volume += 1

    dual_volumes = np.array(dual_volumes, dtype='O')

    return dual_volumes

def create_structure(dual_id: np.ndarray, fine_id: np.ndarray, adjacencies: np.ndarray, ijks: np.ndarray, verify: np.ndarray):
    
    dual = dual_id.astype(np.int64)
    fineid = fine_id.astype(np.int64)
    adjacencies_local = adjacencies.astype(np.int64)

    dual_volumes = get_dual_volumes(dual, fineid, adjacencies_local, ijks, verify)
    np.save('dual_volumes.npy', dual_volumes)


def load_spe_data():
    my_dict_names = local_names()
    spe_data = SuperArrayManager(my_dict_names['spe_data_name'])
    spe_data.load_data()
    return spe_data

    

        

def map_spe_to_fine_volumes(spe_data: SuperArrayManager, ijks, verify):
    
    ijk_spe = spe_data['ijk']
    
    map_spe_to_fine_volumes = np.repeat(-1, ijk_spe.shape[0])
    range_ijks = np.arange(ijks.shape[0])

    # import pdb; pdb.set_trace()
    
    for i in range_ijks[verify]:
        
        # ijk = ijks[i]
        
        # t0 = ijk_spe[:, 0] == ijk[0]
        # t1 = ijk_spe[:, 1] == ijk[1]
        # t2 = ijk_spe[:, 2] == ijk[2]
        # test = t0 & t1 & t2
        # test2 = np.all(ijk_spe == ijks[i], axis=1)
        
        # import pdb; pdb.set_trace()
        map_spe_to_fine_volumes[np.all(ijk_spe == ijks[i], axis=1)] = i
        
        
        print(i)
    
    spe_data.insert_or_update_data({
        'map_spe_to_fine_volumes': map_spe_to_fine_volumes
    })
    spe_data.export_data()
    
    import pdb; pdb.set_trace()
        


        
def mount_graph(A: sp.csc_matrix) -> nx.Graph:
    A2 = mount_matrix_for_graph(A)
    G = nx.from_scipy_sparse_array(A2)
    return G

def mount_graph_with_weight(A: sp.csc_matrix) -> nx.Graph:
    A2 = mount_matrix_for_graph_with_weight(A)
    G = nx.from_scipy_sparse_array(A2)
    return G

def mount_matrix_for_graph(A: sp.csc_matrix) -> sp.csc_matrix:
    # A2: sp.csc_matrix = A.copy()
    A2: sp.csc_matrix = copy.deepcopy(A)
    A2.setdiag(np.zeros(A.shape[0]))
    A2.eliminate_zeros()
    A2.data[:] = 1.0
    return A2.tocsr()

def mount_matrix_for_graph_with_weight(A: sp.csc_matrix) -> sp.csc_matrix:
    A2: sp.csc_matrix = A.copy()
    A2.setdiag(np.zeros(A.shape[0]))
    A2.eliminate_zeros()
    return A2.tocsr()

def create_partition(A: sp.csc_matrix, nparts: int, disjointed=False, nvols_mean=0) -> np.ndarray:
    G = mount_graph(A)
    # G = nx.from_scipy_sparse_array(A, edge_attribute='weight')
    
    if disjointed==True:
        if nvols_mean == 0:
            raise ValueError
        
        partitions = np.repeat(-1, A.shape[0])
        components = nx.connected_components(G)
        subgraphs: Sequence[nx.Graph] = []
        comps = []
        for c in components:
            comps.append(np.array(list(c)))
            subgraphs.append(G.subgraph(c).copy())
        
        partition_id = 0
        
        for sg in subgraphs:
            nnodes = len(sg)
            nodes = np.array([j for j in sg.nodes()])
            
            if nnodes <= nvols_mean:
                part = np.repeat(0, nnodes)
            else:
                nparts2 = int(nnodes/nvols_mean) + 1
                part = np.array(metis.part_graph(sg, nparts=nparts2, contig=True)[1])
            
            part += partition_id
            partitions[nodes] = part
            partition_id = part.max() + 1
            
        return utils_old.test_primal_id(partitions)
    
    
    
    
    # G = mount_graph_with_weight(A)
    # G.graph['edge_weight_attr']='weight'
    
    [cost, part_vert] = metis.part_graph(G, nparts=nparts, recursive=True, contig=True)
    
    # Gm = metis.networkx_to_metis(G)
    # [cost, part_vert] = metis.part_graph(
    #     Gm, 
    #     nparts=nparts, 
    #     recursive=False, 
    #     contig=True,
    #     minconn=True,
    #     objtype='cut'
    # )
    
    return utils_old.test_primal_id(np.array(part_vert))

def connect_graph_for_metis(A_filtered, A_original, eps=1e-10):
    """
    Garante que o grafo seja conexo para o METIS sem alterar a física.
    A_filtered: Matriz após filtro 
    A_original: Matriz de transmissibilidade completa.
    eps: Peso minúsculo para as conexões reintegradas.
    """
    # 1. Identificar onde a matriz original tem conexões mas a filtrada não tem
    # Usamos a estrutura (indices) da original
    S_full = (A_original != 0).astype(bool)
    S_filt = (A_filtered != 0).astype(bool)
    
    # 2. Criar uma máscara das conexões que foram removidas pelo filtro
    S_gap = S_full > S_filt
    
    avg_trans = np.mean(np.abs(A_original.data))
    low_weight = avg_trans * eps
    
    A_background = S_gap.astype(float)
    A_background.data *= low_weight
    
    # 4. Somar as duas: O METIS agora vê um grafo único, mas as conexões 
    # de Watanabe são ordens de grandeza mais "atraentes" para não serem cortadas.
    A_for_metis = A_filtered + A_background
    
    return A_for_metis

    
    
    
    
    
    
    # if debug == True:
    #     import matplotlib.pyplot as plt
    #     pos = nx.spring_layout(G, seed=4)
    #     nx.draw(
    #         G, pos, 
    #         node_color=part_vert, 
    #         with_labels=True, 
    #         cmap=plt.cm.tab10,      # Use a qualitative colormap for distinct groups
    #         node_size=600, 
    #         edge_color="black", 
    #         alpha=0.8
    #     )

    #     plt.savefig('graph.png')
    #     import pdb; pdb.set_trace()
    
    return utils_old.test_primal_id(np.array(part_vert))

def create_partition_v0dep(A: sp.csc_matrix, nparts: int, disjointed=False, nvols_mean=0) -> np.ndarray:
    G = mount_graph(A)
    
    if disjointed==True:
        if nvols_mean == 0:
            raise ValueError
        
        partitions = np.repeat(-1, A.shape[0])
        components = nx.connected_components(G)
        subgraphs = []
        comps = []
        for c in components:
            comps.append(np.array(list(c)))
            subgraphs.append(G.subgraph(c).copy())
        n2 = len(subgraphs)
        n3 = int(nparts/n2)
        if n3 < 2:
            n3 = 2   
        # comps = [np.array(list(c)) for c in components]
        # subgraphs = [G.subgraph(c).copy() for c in components]
        partition_id = 0
        for sg in subgraphs:
            nnodes = len(sg)
            nodes = np.array([j for j in sg.nodes()])
            # test_cont = True
            # cont_local = 0
            # while cont_local < n3-1 and test_cont == True:
            #     prox = cont_local + 1
            #     if cont_local*nparts < nnodes <= prox*nparts:
            #         test_cont = False
            #         part = np.array(metis.part_graph(sg, nparts=prox, recursive=True, contig=True)[1])
                    
            if nnodes > 2*nvols_mean and nnodes <= 3*nvols_mean:
                part = np.array(metis.part_graph(sg, nparts=3, recursive=True, contig=True)[1])
                # part += partition_id
                # partitions[nodes] = part
                # partition_id = part.max() + 1
            elif nnodes > nvols_mean and nnodes <= 2*nvols_mean:
                part = np.array(metis.part_graph(sg, nparts=2, recursive=True, contig=True)[1])
                # part += partition_id
                # partitions[nodes] = part
                # partition_id = part.max() + 1
            elif nnodes <= nvols_mean:
                part = np.repeat(0, nnodes)
                # part += partition_id
                # partitions[nodes] = part
                # partition_id = part.max() + 1
                
                # part = nx.community.asyn_fluidc(sg, k=3, seed=1)
                # for pid, nodes in enumerate(part):
                #     partitions[np.array(list(nodes))] = partition_id
                #     partition_id += 1
            else:
                part = np.array(metis.part_graph(sg, nparts=n3, recursive=True, contig=True)[1])
                # part += partition_id
                # partitions[nodes] = part
                # partition_id = part.max() + 1
            
            part += partition_id
            partitions[nodes] = part
            partition_id = part.max() + 1
        
        # fineids = np.unique(np.concatenate(comps))
        # new_partition = fineids.copy()
        # new_partition[:] = -1
        
        # j = -1
        # for i, part in enumerate(partitions):
        #     unique_local_parts = np.unique(part)
        #     local_fine_ids = comps[i]
        #     for id2 in unique_local_parts:
        #         j += 1
        #         test = part == id2
        #         new_partition[local_fine_ids[test]] = j
        # return new_partition
        return utils_old.test_primal_id(partitions)
    
    
    
    
    # G = mount_graph_with_weight(A)
    # G.graph['edge_weight_attr']='weight'
    # para criar a particao os pesos devem ser inteiros
    [cost, part_vert] = metis.part_graph(G, nparts=nparts, recursive=True, contig=True)
    
    
    
    
    
    
    # if debug == True:
    #     import matplotlib.pyplot as plt
    #     pos = nx.spring_layout(G, seed=4)
    #     nx.draw(
    #         G, pos, 
    #         node_color=part_vert, 
    #         with_labels=True, 
    #         cmap=plt.cm.tab10,      # Use a qualitative colormap for distinct groups
    #         node_size=600, 
    #         edge_color="black", 
    #         alpha=0.8
    #     )

    #     plt.savefig('graph.png')
    #     import pdb; pdb.set_trace()
    
    return utils_old.test_primal_id(np.array(part_vert))

def create_vertices(primal_id: np.ndarray, A: sp.csc_matrix) -> np.ndarray:
    
    cids = np.unique(primal_id)
    fine_ids = np.arange(A.shape[0])
    diagonal_term = np.zeros(A.shape[0])
    all_vertices = np.repeat(-1, cids.shape[0])
    
    for cid in cids:
        fine_ids_in_cid = fine_ids[primal_id==cid]
        local_matrix = utils_old.get_local_matrix(fine_ids_in_cid, A, diagonal_term)
        G_local = mount_graph(local_matrix)
        # nao ponderado: nx.closeness_centrality(G) - usa BFS.
        # Ponderado: nx.closeness_centrality(G, distance="weight") - usa Dijkstra.
        guess = nx.closeness_centrality(G_local)
        # chaves = np.array(list(guess.keys()))
        values = np.array(list(guess.values()))
        closest = fine_ids_in_cid[values<=values.min()][0]
        all_vertices[cid] = closest
        # guess2 = nx.group_closeness_centrality(G, fine_ids_in_cid)
    
    return all_vertices

def create_local_support_region(data):
    fine_ids, primal_id, A2, n_levels_adj, cid = data
    fine_ids_in_cid = fine_ids[primal_id==cid]
    lines_adj = A2[fine_ids_in_cid].toarray()
    support_region = reduce(np.union1d, [fine_ids[i] for i in lines_adj])
    support_region = np.setdiff1d(support_region, fine_ids_in_cid)
    
    for i in range(n_levels_adj-1):
        lines_adj = A2[support_region].toarray()
        support_region = reduce(np.union1d,  [fine_ids[i] for i in lines_adj])
        
    support_region = np.union1d(support_region, fine_ids_in_cid)

    lines_adj = A2[support_region].toarray()
    # interaction = np.unique(np.concatenate([fine_ids[i] for i in lines_adj]))
    interaction = reduce(np.union1d, [fine_ids[i] for i in lines_adj])
    boundary = np.setdiff1d(interaction, support_region)
    
    return np.array([[cid], support_region, interaction, boundary], dtype='O')

def create_local_support_region_adj(data):
    
    fine_ids, primal_id, A2, n_levels_adj, cid = data
    fine_ids_in_cid = fine_ids[primal_id==cid]
    # support_region = reduce(np.union1d, A2[fine_ids_in_cid])
    support_region = np.unique(np.concatenate(A2[fine_ids_in_cid]))
    support_region = np.setdiff1d(support_region, fine_ids_in_cid)
    
    if support_region.shape[0] > 0:
        for i in range(n_levels_adj-1):
            # support_region = reduce(np.union1d,  A2[support_region])
            support_region = np.unique(np.concatenate(A2[support_region]))
        
    support_region = np.union1d(support_region, fine_ids_in_cid)

    # interaction = reduce(np.union1d, A2[support_region])
    interaction = np.unique(np.concatenate(A2[support_region]))
    boundary = np.setdiff1d(interaction, support_region)
    
    # print(f'PID: {os.getpid()}, CID: {cid}, dt: {dt} \n')
    
    # print(f' PID: {os.getpid()} \n')
    # print(f' CID: {cid} \n')
    
    return np.array([[cid], support_region, interaction, boundary], dtype='O')
     
def args_support_generator(fine_ids, primal_id, A2, n_levels_adj, cids):
    
    for cid in cids:
        yield (fine_ids, primal_id, A2, n_levels_adj, cid)

def args_volume_volume_adj_generator(fine_ids, A2):
    for i in fine_ids:
        yield (A2, i, fine_ids)

def create_volume_volume_adj(data):
    A2, i, fine_ids = data
    return fine_ids[A2[i].toarray().flatten()]
    
def get_boundary_of_suport_region(suporte1d: np.ndarray, ptr_suporte: np.ndarray, interacao1d: np.ndarray, ptr_interacao: np.ndarray):
    ncoarse_vols = len(ptr_suporte)-1
    tam_int = np.diff(ptr_interacao)
    tam_sup = np.diff(ptr_suporte)
    tamanhos_bordas = tam_int - tam_sup
    ptr_borda = np.zeros(len(ptr_interacao), dtype=int)
    ptr_borda[1:] = np.cumsum(tamanhos_bordas)
    dados_borda_1d = np.empty(ptr_borda[-1], dtype=int)
    
    for vol in range(ncoarse_vols):
        sup = suporte1d[ptr_suporte[vol]:ptr_suporte[vol+1]]
        inter = interacao1d[ptr_interacao[vol]:ptr_interacao[vol+1]]
        dados_borda_1d[ptr_borda[vol] : ptr_borda[vol+1]] = np.setdiff1d(inter, sup)
    
    boundary = np.array(np.split(dados_borda_1d, ptr_borda[1:-1]), dtype='O')
    
    return boundary

def get_boundary_of_suport_region_v2(suporte1d: np.ndarray, ptr_suporte: np.ndarray, interacao1d: np.ndarray, ptr_interacao: np.ndarray):
    ncoarse_vols = len(ptr_suporte)-1
    tam_int = np.diff(ptr_interacao)
    tam_sup = np.diff(ptr_suporte)
    tamanhos_bordas = tam_int - tam_sup
    ptr_borda = np.zeros(len(ptr_interacao), dtype=int)
    ptr_borda[1:] = np.cumsum(tamanhos_bordas)
    dados_borda_1d = np.empty(ptr_borda[-1], dtype=int)
    
    for vol in range(ncoarse_vols):
        sup = suporte1d[ptr_suporte[vol]:ptr_suporte[vol+1]]
        inter = interacao1d[ptr_interacao[vol]:ptr_interacao[vol+1]]
        dados_borda_1d[ptr_borda[vol] : ptr_borda[vol+1]] = np.setdiff1d(inter, sup)
    
    return dados_borda_1d, ptr_borda

@time_func
def create_support_region_and_boundary(primal_id: np.ndarray, volume_adjacencies: np.ndarray, n_levels_adj=3, ext='', **kwargs) -> dict:
    cids = np.unique(primal_id)
    fine_ids = np.arange(primal_id.shape[0])
    all_support_region = []
    all_boundary_region = []
    all_interaction_region = []
    
    
    
    ncpus = mp.cpu_count()
    ncpus = 5
    # rprint(f'[green] N cpus: {ncpus}')
    if ncpus > 1:
        ncpus -= 1
        
    # tasks_generator_adj = args_volume_volume_adj_generator(fine_ids, A2)
    
    # with mp.Pool(processes=ncpus) as pool:
    #     results = pool.imap(create_volume_volume_adj, tasks_generator_adj)
    #     for res in results:
    #         volume_volume_adjacency.append(res)
    
    # volume_volume_adjacency = np.array(volume_volume_adjacency, dtype='O')
    
    # ext_npy = '.npy'
    all_support_region_str = 'all_support_region' + ext
    all_interaction_region_str = 'all_interaction_region' + ext
    all_boundary_region_str = 'all_boundary_region' + ext
    
    # task_generator = args_support_generator(fine_ids, primal_id, volume_volume_adjacency, n_levels_adj, cids)
    # task_generator = args_support_generator(fine_ids, primal_id, A2, n_levels_adj, cids)
    task_generator = args_support_generator(fine_ids, primal_id, volume_adjacencies, n_levels_adj, cids)

    all_times = []
    count_loop = 0
    for task in task_generator:
        start_time = time.perf_counter()
        res = create_local_support_region_adj(task)
        end_time = time.perf_counter()
        all_support_region.append(res[1])
        all_interaction_region.append(res[2])
        all_boundary_region.append(res[3])
        elapsed_time = end_time - start_time
        all_times.append(elapsed_time)
        print(f'CID: {count_loop}, Time: {elapsed_time:.4f} seconds')
        count_loop += 1
    
        
    #     if count_loop >= 20:
    #         break
    
    # avg_time = sum(all_times) / len(all_times)
    
    # import pdb; pdb.set_trace()
    
    
    
    
   
    
    # all_cids = []
    # with mp.Pool(processes=ncpus) as pool:
    #     # results = pool.imap_unordered(create_local_support_region_adj, task_generator)
    #     # results = pool.starmap_async(create_local_support_region, task_generator).get()
    #     # results = pool.starmap(create_local_support_region_adj, task_generator)
    #     # results = pool.imap(create_local_support_region_adj, task_generator)
    #     results = pool.map(create_local_support_region_adj, task_generator)
        
    #     # all_results = list(results)
        
    #     # for res in results:
    #     #     all_support_region.append(res[1])
    #     #     all_interaction_region.append(res[2])
    #     #     all_boundary_region.append(res[3])
    #         # all_cids.append(res[0][0])
            
    
    # # import pdb; pdb.set_trace()
    
    # all_support_region = np.array([res[1] for res in results], dtype='O')
    # all_boundary_region = np.array([res[3] for res in results], dtype='O')
    # all_interaction_region = np.array([res[2] for res in results], dtype='O')

    # for cid in track(cids, description='[red] Computing support regions ...[/red]'):
    # # for cid in cids:
    #     fine_ids_in_cid = copy.deepcopy(fine_ids[primal_id==cid])
    #     lines_adj = A2[fine_ids_in_cid].toarray()
    #     support_region = reduce(np.union1d, [fine_ids[j] for j in lines_adj])
    #     support_region = np.setdiff1d(support_region, fine_ids_in_cid)
        
    #     for i in range(n_levels_adj-1):
    #         lines_adj = A2[support_region].toarray()
    #         support_region = reduce(np.union1d, [fine_ids[j] for j in lines_adj])
        
    #     support_region = np.union1d(support_region, fine_ids_in_cid)

    #     lines_adj = A2[support_region].toarray()
    #     interaction = np.unique(np.concatenate([fine_ids[i] for i in lines_adj]))
    #     boundary = np.setdiff1d(interaction, support_region)
        
    #     all_support_region.append(support_region)
    #     all_interaction_region.append(interaction)
    #     all_boundary_region.append(boundary)
    # #     # rprint('[red] Computing support regions ...[/red] \n')
    # #     # print(f'{(cid/cids.max()):.5f}%')
    
    all_support_region = np.array(all_support_region, dtype='O')
    all_interaction_region = np.array(all_interaction_region, dtype='O')
    all_boundary_region = np.array(all_boundary_region, dtype='O')
    
    # np.save(all_support_region_str + ext_npy, all_support_region)
    # np.save(all_interaction_region_str + ext_npy, all_interaction_region)
    # np.save(all_boundary_region_str + ext_npy, all_boundary_region)
    
    resp = {
        all_support_region_str: all_support_region,
        all_interaction_region_str: all_interaction_region,
        all_boundary_region_str: all_boundary_region
    }
    
    return resp

def create_support_region_and_boundary_v2(A: sp.csc_matrix, primal_id: np.ndarray, n_levels_adj=3, ext='', **kwargs) -> dict:
    
    all_support_region_str = 'all_support_region' + ext
    all_interaction_region_str = 'all_interaction_region' + ext
    all_boundary_region_str = 'all_boundary_region' + ext
    
    A2: sp.csr_matrix = copy.deepcopy(A).tocsr()
    A2.data[:] = 1.0

    expansao = get_OR_finite_volume(primal_id).tocsr()
    for _ in range(n_levels_adj):
        expansao = expansao @ A2
    
    suporteM: sp.csr_matrix = copy.deepcopy(expansao)
    all_support_region = np.array(np.split(expansao.indices, expansao.indptr[1:-1]), dtype='O')
    expansao = expansao @ A2
    interactionM = expansao
    all_interaction_region = np.array(np.split(expansao.indices, expansao.indptr[1:-1]), dtype='O')
    all_boundary_region = get_boundary_of_suport_region(suporteM.indices, suporteM.indptr, interactionM.indices, interactionM.indptr)
        
    resp = {
        all_support_region_str: all_support_region,
        all_interaction_region_str: all_interaction_region,
        all_boundary_region_str: all_boundary_region
    }
    
    return resp

def create_support_region_and_boundary_v3(A: sp.csc_matrix, primal_id: np.ndarray, n_levels_adj=3, ext='', **kwargs) -> dict:
    
    A2: sp.csr_matrix = copy.deepcopy(A).tocsr()
    A2.data[:] = 1.0

    expansao = get_OR_finite_volume(primal_id).tocsr()
    for _ in range(n_levels_adj):
        expansao = expansao @ A2
    
    suporteM: sp.csr_matrix = copy.deepcopy(expansao)
    expansao = expansao @ A2
    interactionM: sp.csr_matrix = expansao
    ind_boundary, ptr_boundary = get_boundary_of_suport_region_v2(suporteM.indices, suporteM.indptr, interactionM.indices, interactionM.indptr)
        
    resp = {
        'ind_support': copy.deepcopy(suporteM.indices),
        'ptr_support': copy.deepcopy(suporteM.indptr),
        'ind_interaction': copy.deepcopy(interactionM.indices),
        'ptr_interaction': copy.deepcopy(interactionM.indptr),
        'ind_boundary': ind_boundary,
        'ptr_boundary': ptr_boundary
    }
    
    return resp

@time_func
def create_support_region_and_boundary_sparse(primal_id: np.ndarray, A: sp.csc_matrix, n_levels_adj=3, ext='', debug=False, **kwargs) -> dict:
    cids = np.unique(primal_id)
    fine_ids = np.arange(primal_id.shape[0])
    all_support_region = []
    all_boundary_region = []
    all_interaction_region = []
    A2 = mount_matrix_for_graph(A).astype(np.bool_)
    
    
    
    ncpus = mp.cpu_count()
    ncpus = 5
    # rprint(f'[green] N cpus: {ncpus}')
    if ncpus > 1:
        ncpus -= 1
        
    # tasks_generator_adj = args_volume_volume_adj_generator(fine_ids, A2)
    
    # with mp.Pool(processes=ncpus) as pool:
    #     results = pool.imap(create_volume_volume_adj, tasks_generator_adj)
    #     for res in results:
    #         volume_volume_adjacency.append(res)
    
    # volume_volume_adjacency = np.array(volume_volume_adjacency, dtype='O')
    
    # ext_npy = '.npy'
    all_support_region_str = 'all_support_region' + ext
    all_interaction_region_str = 'all_interaction_region' + ext
    all_boundary_region_str = 'all_boundary_region' + ext
    
    # task_generator = args_support_generator(fine_ids, primal_id, volume_volume_adjacency, n_levels_adj, cids)
    task_generator = args_support_generator(fine_ids, primal_id, A2, n_levels_adj, cids)
    # task_generator = args_support_generator(fine_ids, primal_id, volume_adjacencies, n_levels_adj, cids)
    
    # all_times = []
    # count_loop = 0
    # for task in task_generator:
    #     start_time = time.perf_counter()
    #     res = create_local_support_region_adj(task)
    #     end_time = time.perf_counter()
    #     all_support_region.append(res[1])
    #     all_interaction_region.append(res[2])
    #     all_boundary_region.append(res[3])
    #     elapsed_time = end_time - start_time
    #     all_times.append(elapsed_time)
    #     count_loop += 1
        
    #     if count_loop >= 20:
    #         break
    
    # avg_time = sum(all_times) / len(all_times)
    
    # import pdb; pdb.set_trace()
    
    
    
    
   
    
    # all_cids = []
    with mp.Pool(processes=ncpus) as pool:
        # results = pool.imap_unordered(create_local_support_region_adj, task_generator)
        # results = pool.starmap_async(create_local_support_region, task_generator).get()
        # results = pool.starmap(create_local_support_region_adj, task_generator)
        # results = pool.imap(create_local_support_region_adj, task_generator)
        results = pool.map(create_local_support_region, task_generator)
        
        # all_results = list(results)
        
        # for res in results:
        #     all_support_region.append(res[1])
        #     all_interaction_region.append(res[2])
        #     all_boundary_region.append(res[3])
            # all_cids.append(res[0][0])
            
    
    # import pdb; pdb.set_trace()
    
    all_support_region = np.array([res[1] for res in results], dtype='O')
    all_boundary_region = np.array([res[3] for res in results], dtype='O')
    all_interaction_region = np.array([res[2] for res in results], dtype='O')

    # for cid in track(cids, description='[red] Computing support regions ...[/red]'):
    # # for cid in cids:
    #     fine_ids_in_cid = copy.deepcopy(fine_ids[primal_id==cid])
    #     lines_adj = A2[fine_ids_in_cid].toarray()
    #     support_region = reduce(np.union1d, [fine_ids[j] for j in lines_adj])
    #     support_region = np.setdiff1d(support_region, fine_ids_in_cid)
        
    #     for i in range(n_levels_adj-1):
    #         lines_adj = A2[support_region].toarray()
    #         support_region = reduce(np.union1d, [fine_ids[j] for j in lines_adj])
        
    #     support_region = np.union1d(support_region, fine_ids_in_cid)

    #     lines_adj = A2[support_region].toarray()
    #     interaction = np.unique(np.concatenate([fine_ids[i] for i in lines_adj]))
    #     boundary = np.setdiff1d(interaction, support_region)
        
    #     all_support_region.append(support_region)
    #     all_interaction_region.append(interaction)
    #     all_boundary_region.append(boundary)
    # #     # rprint('[red] Computing support regions ...[/red] \n')
    # #     # print(f'{(cid/cids.max()):.5f}%')
    
    # all_support_region = np.array(all_support_region, dtype='O')
    # all_interaction_region = np.array(all_interaction_region, dtype='O')
    # all_boundary_region = np.array(all_boundary_region, dtype='O')
    
    # np.save(all_support_region_str + ext_npy, all_support_region)
    # np.save(all_interaction_region_str + ext_npy, all_interaction_region)
    # np.save(all_boundary_region_str + ext_npy, all_boundary_region)
    
    resp = {
        all_support_region_str: all_support_region,
        all_interaction_region_str: all_interaction_region,
        all_boundary_region_str: all_boundary_region
    }
    
    return resp

def test_same_support_regions(support1: np.ndarray, support2: np.ndarray, debug=False):
    
    if debug == True:
        assert len(support1) == len(support2)
        rprint('[blue on white]Regioes de suporte tem o mesmo tamanho[/blue on white]')
        # rg: 1464542 SDS
        sup1 = []
        sup2 = []
        for s1, s2 in zip(support1, support2):
            sup1.append(np.unique(s1))
            sup2.append(np.unique(s2))
        
        s1 = np.array(sup1, dtype='O')
        s2 = np.array(sup2, dtype='O')
        
        v1 = np.concatenate(s1)
        v2 = np.concatenate(s2)
        assert np.all(v1 == v2)
        rprint('[blue on white]Sao as mesmas regioes de suporte[/blue on white]')
        
     

def create_support_region_and_boundary_barrier(primal_id: np.ndarray, fine_ids_in_barrier: np.ndarray, A: sp.csc_matrix, n_levels_adj=3, ext='') -> dict:
    cids = np.unique(primal_id)
    fine_ids = np.arange(A.shape[0])
    # diagonal_term = np.zeros(A.shape[0])
    all_support_region = []
    all_boundary_region = []
    all_interaction_region = []
    A2 = mount_matrix_for_graph(A).astype(np.bool_)

    n_levels_adj_iterate = n_levels_adj
    
    ext_npy = '.npy'
    all_support_region_str = 'all_support_region' + ext
    all_interaction_region_str = 'all_interaction_region' + ext
    all_boundary_region_str = 'all_boundary_region' + ext
    
    for cid in cids:
        support_region = fine_ids[primal_id==cid].copy()
        test_in = np.isin(support_region, fine_ids_in_barrier)
        if np.any(test_in):
            n_levels_adj_iterate = 5
        else:
            n_levels_adj_iterate = n_levels_adj
        for i in range(n_levels_adj_iterate):
            lines_adj = A2[support_region].toarray()
            support_region = np.unique(np.concatenate([fine_ids[j] for j in lines_adj]))

        lines_adj = A2[support_region].toarray()
        interaction = np.unique(np.concatenate([fine_ids[i] for i in lines_adj]))
        boundary = np.setdiff1d(interaction, support_region)
        
        all_support_region.append(support_region)
        all_interaction_region.append(interaction)
        all_boundary_region.append(boundary)
        print(cid)
    
    all_support_region = np.array(all_support_region, dtype='O')
    all_interaction_region = np.array(all_interaction_region, dtype='O')
    all_boundary_region = np.array(all_boundary_region, dtype='O')
    
    np.save(all_support_region_str + ext_npy, all_support_region)
    np.save(all_interaction_region_str + ext_npy, all_interaction_region)
    np.save(all_boundary_region_str + ext_npy, all_boundary_region)
    
    resp = {
        all_support_region_str: all_support_region,
        all_interaction_region_str: all_interaction_region,
        all_boundary_region_str: all_boundary_region
    }
    
    return resp

def get_enhanced_matrix(T: sp.csc_matrix):
    n = T.shape[0]
    all_data = sp.find(T)
    test_neg = all_data[2] < 0
    test_diag = all_data[0] == all_data[1]
    off_diag = ~test_diag
    
    test = off_diag & test_neg
    
    lines_neg = all_data[0][test]
    values_neg = all_data[2][test]
    cols_neg = all_data[1][test]
    
    index, idx = np.unique(lines_neg, return_inverse=True)
    soma = np.bincount(idx, weights=values_neg)
    
    lines = np.concatenate([lines_neg, index])
    cols = np.concatenate([cols_neg, index])
    data = np.concatenate([values_neg, -1*soma])
    
    # matrix = sp.csc_matrix((data, (lines, cols)), shape=(n, n))
    # return matrix
    
    return sp.csc_matrix((data, (lines, cols)), shape=(n, n))
           
def _get_D_matrix(T: sp.csc_matrix) -> sp.csc_matrix:
    n = T.shape[0]
    diag = T.diagonal()
    test = np.absolute(diag) <= 1e-14
    if np.any(test):
        print(diag[test])
        raise ValueError('Diagonal has zero or near zero entries, cannot compute D matrix.')
    D = sp.spdiags(1/diag, 0, n, n).tocsc()        
    return D

def _get_Ematrix(T: sp.csc_matrix, omega: float) -> sp.csr_matrix:
    m1: sp.csc_matrix = -omega*(_get_D_matrix(T)@T)
    m1.setdiag(m1.diagonal() + 1)
    return m1.tocsr(copy=True)

def get_OR_finite_volume(primal_id: np.ndarray) -> sp.csc_matrix:
    lines = primal_id
    cols = np.arange(primal_id.shape[0])
    data = np.repeat(1.0, cols.shape[0])
    n_coarse = np.unique(primal_id).shape[0]
    n_fine = primal_id.shape[0]
    OR = sp.csc_matrix((data, (lines, cols)), shape=(n_coarse, n_fine))
    return OR

def mount_indices_support_regions(primal_id, support_regions, debug=False):
    cids = np.unique(primal_id)
    coarse_ids_support_regions = []
    for cid in cids:
        support_region = support_regions[cid]
        coarse_ids_support_regions.append(np.repeat(cid, support_region.shape[0]))
    
    coarse_ids_support_regions = np.array(coarse_ids_support_regions, dtype='O')
    return coarse_ids_support_regions
    
    
def _mount_local_op_it(support_regions: np.ndarray, boundary_regions: np.ndarray, all_boundary_regions: np.ndarray, coarse_ids_support_regions: np.ndarray, coarse_ids_boundary_regions: np.ndarray, Ematrix: sp.csc_matrix, OP0: sp.csc_matrix) -> sp.csc_matrix:
    
    new_OP = Ematrix@OP0

    lines = np.concatenate(support_regions)
    cols = np.concatenate(coarse_ids_support_regions)

    new_OP2 = sp.lil_matrix((new_OP.shape))
    new_OP2[lines, cols] = new_OP[lines, cols]
    new_OP2 = new_OP2.tocsc()
    soma = np.array(new_OP2.sum(axis=1)).flatten()
    soma[:] = 1/soma
    
    indices_all_boundary_region = np.isin(new_OP2.indices, all_boundary_regions)
    # soma_all_boundary_regions = soma[new_OP.indices[indices_boundary_region]]

    new_OP2.data[indices_all_boundary_region] *= soma[new_OP2.indices[indices_all_boundary_region]]
    error = np.absolute((OP0 - new_OP2)[lines, cols].data)
    
    return new_OP2, error.max()

def _mount_local_op_it_strong(support_regions: np.ndarray, boundary_regions: np.ndarray, all_boundary_regions: np.ndarray, coarse_ids_support_regions: np.ndarray, coarse_ids_boundary_regions: np.ndarray, Ematrix: sp.csc_matrix, OP0: sp.csc_matrix, debug=False) -> sp.csc_matrix:
    
    new_OP = Ematrix@OP0

    lines = np.concatenate(support_regions)
    cols = np.concatenate(coarse_ids_support_regions)
    
    new_OP2 = sp.lil_matrix((new_OP.shape))
    new_OP2[lines, cols] = new_OP[lines, cols]
    new_OP2: sp.csc_matrix = new_OP2.tocsc()
    new_OP2.eliminate_zeros()
    soma = np.array(new_OP2.sum(axis=1)).flatten()
    soma[:] = 1/soma
    
    # indices_all_boundary_region = np.isin(new_OP2.indices, all_boundary_regions)
    # soma_all_boundary_regions = soma[new_OP.indices[indices_boundary_region]]

    new_OP2.data *= soma[new_OP2.indices]
    error = np.absolute((OP0 - new_OP2)[lines, cols].data)
    # error = np.absolute((OP0 - new_OP2).data)
    
    return new_OP2, error.max()

def _restrict_P_to_support_region(P: sp.csc_matrix, ind_support: np.ndarray, ptr_support: np.ndarray):
    """
        P: Matriz do operador de prolongamento iterada.
        ind_support, ptr_support : Estrutura 1D + ptr da região de suporte.
    """
    novos_data = []
    novos_indices = []
    novos_indptr = [0]
    ncoarse = len(ptr_support) - 1
    
    for i in range(ncoarse):
        # Suporte do volume grosso i
        suporte_i = ind_support[ptr_support[i] : ptr_support[i+1]]
        col_i = P.getcol(i).tocoo() # COO para pegar indices/data
        
        # Filtro
        mask = np.isin(col_i.row, suporte_i)
        
        novos_data.append(col_i.data[mask])
        novos_indices.append(col_i.row[mask])
        novos_indptr.append(len(col_i.data[mask]))
    
    novos_data = np.concatenate(novos_data)
    novos_indices = np.concatenate(novos_indices)
    novos_indptr = np.cumsum(novos_indptr)
    
    P_new = sp.csc_matrix((novos_data, novos_indices, novos_indptr), shape=P.shape)
    return P_new

def _mount_local_op_it_strong_v3(Msupport: sp.csc_matrix, Ematrix: sp.csc_matrix, OP0: sp.csc_matrix, debug=False) -> sp.csc_matrix:
    
    new_OP = Msupport.multiply(Ematrix@OP0)
    soma = np.array(new_OP.sum(axis=1)).flatten()
    soma[:] = 1/soma

    new_OP.data *= soma[new_OP.indices]
    error = np.absolute((OP0 - new_OP).data)
    
    return new_OP, error.max()

def _mount_local_op_it_strong_v2(ind_support: np.ndarray, ptr_support: np.ndarray, Ematrix: sp.csc_matrix, OP0: sp.csc_matrix, debug=False) -> sp.csc_matrix:
    
    new_OP = Ematrix@OP0
    new_OP = _restrict_P_to_support_region(new_OP.tocsc(), ind_support, ptr_support)
    soma = np.array(new_OP.sum(axis=1)).flatten()
    soma[:] = 1/soma

    new_OP.data *= soma[new_OP.indices]
    error = np.absolute((OP0 - new_OP).data)
    
    return new_OP, error.max()

def _mount_local_op_it_strong_cor(support_regions: np.ndarray, boundary_regions: np.ndarray, all_boundary_regions: np.ndarray, coarse_ids_support_regions: np.ndarray, coarse_ids_boundary_regions: np.ndarray, Ematrix: sp.csc_matrix, OP0: sp.csc_matrix, debug=False) -> sp.csc_matrix:
    
    new_OP = (Ematrix@OP0).tocsc()

    lines = np.concatenate(support_regions)
    cols = np.concatenate(coarse_ids_support_regions)
    
    soma = np.array(new_OP.sum(axis=1)).flatten()
    soma[:] = 1/soma
    new_OP.data *= soma[new_OP.indices]
    
    new_OP2 = sp.lil_matrix((new_OP.shape))
    new_OP2[lines, cols] = new_OP[lines, cols]
    
    error = np.absolute((OP0 - new_OP2).data)
    
    return new_OP2.tocsc(), error.max()

def get_modified_matrix(A: sp.csc_matrix) -> sp.csc_matrix:
    '''
    A Multiscale Restriction-Smoothed Basis
    Method for Compressible Black-Oil Models
    '''
    soma = np.array(A.sum(axis=1)).flatten()
    test = soma == 1
    test = ~test
    indices = np.arange(A.shape[0])[test]
    A2 = utils_old.get_matrix_slice(A, indices)
    A2_balanced = get_balanced_matrix(A2)
    new_index = np.arange(A2_balanced.shape[0])
    return A2_balanced, new_index, indices



def define_strong_coupled(A: sp.csc_matrix, eps: float=0.05, **kwargs) -> sp.csc_matrix:
    
    diagA = A.diagonal()

    a_data = sp.find(A)
    test_diag = a_data[0] == a_data[1]
    off_diag = ~test_diag

    lines_off = a_data[0][off_diag]
    values_off = a_data[2][off_diag]
    cols_off = a_data[1][off_diag]

    values_off_abs = np.abs(values_off)
    value_for_test = np.sqrt(np.absolute(eps*diagA[lines_off]*diagA[cols_off]))
    test_strong = values_off_abs > value_for_test
    
    # test_weak = ~test_strong
    vec1 = np.array([lines_off, cols_off]).T
    
    pairs1_strong = vec1[test_strong]
    
    # # view2 = pairs2_weak.view(dt2)
    view2 = np.zeros(pairs1_strong.shape[0], dtype=[('lines', np.int64), ('cols', np.int64)])
    view2['lines'][:] = pairs1_strong[:, 0]
    view2['cols'][:] = pairs1_strong[:, 1]
    
    view3 = np.zeros(lines_off.shape[0], dtype=[('lines', np.int64), ('cols', np.int64)])
    view3['lines'][:] = cols_off
    view3['cols'][:] = lines_off
    
    test_strong2 = np.isin(view3, view2)
    
    test_strong[:] = test_strong | test_strong2
    
    # test_strong = ~test_weak

    new_lines = lines_off[test_strong]
    new_cols = cols_off[test_strong]
    new_data = values_off[test_strong]

    A2 = sp.csc_matrix((new_data, (new_lines, new_cols)), shape=(A.shape[0], A.shape[1]))
    soma = np.array(A2.sum(axis=1)).flatten()
    test = np.absolute(soma) <= 1e-14
    soma[test] = -1.0
    A2.setdiag(-soma)   
    return A2

def define_strong_coupled_v2(A: sp.csc_matrix, eps: float=0.05, **kwargs) -> sp.csc_matrix:
    
    diagA = np.abs(A.diagonal())
    soma0 = np.array(A.sum(axis=1)).flatten()
    A2 = A.tocsr(copy=True)
    A2.setdiag(0)
    A2.eliminate_zeros()
    A2.data[:] = np.absolute(A2.data)
    
    row, col = A2.nonzero()
    is_strong = A2.data > (eps*np.sqrt(diagA[row]*diagA[col]))
    
    M = sp.csr_matrix((is_strong, (row, col)), shape=A.shape, dtype=bool)
    M = M + M.transpose()
    M.data[:] = 1.0
    # A_filtered = M.multiply(A2)
    A3: sp.csr_matrix = A.copy()
    A3.setdiag(0)
    A3.eliminate_zeros()
    A_filtered = M.multiply(A3)
    soma = -np.array(A_filtered.sum(axis=1)).flatten() + soma0
    soma[soma == 0] = -1.0
    A_filtered.setdiag(-soma)
    
    return A_filtered

def define_strong_coupled_v3(A: sp.csc_matrix, theta: float=0.25, **kwargs) -> sp.csc_matrix:
    
    """Novo filtro determinando a matriz com cada strong conection na face
        e pegando o maior valor na linha
    """
    diagA = np.abs(A.diagonal())
    soma0 = np.array(A.sum(axis=1)).flatten()
    
    A2 = A.tocsr(copy=True)
    A2.setdiag(0)
    A2.eliminate_zeros()
    A2.data[:] = np.absolute(A2.data)
    
    row, col = A2.nonzero()
    A2.data = A2.data/(np.sqrt(diagA[row]*diagA[col]))
    
    rowmax = A2.max(axis=1).toarray().flatten()
    is_strong = A2.data >= (theta * rowmax[row])
    
    M = sp.csr_matrix((is_strong, (row, col)), shape=A.shape, dtype=bool)
    M = M + M.transpose()
    M.data[:] = 1.0
    # A_filtered = M.multiply(A2)
    A3: sp.csr_matrix = A.copy()
    A3.setdiag(0)
    A3.eliminate_zeros()
    A_filtered = M.multiply(A3)
    soma = -np.array(A_filtered.sum(axis=1)).flatten() + soma0
    soma[soma == 0] = -1.0
    A_filtered.setdiag(-soma)
    
    return A_filtered



def compute_directional_strength_matrix(A, epsilon=1e-10, gamma=0.9):
    """
    Implementa o algoritmo de manter as conexoes pelo valor acumulado por linha
    
    Parâmetros:
    A : scipy.sparse matrix
        Matriz esparsa de transmissibilidade (coeficientes a_ij).
    epsilon : float
        Pequeno valor para evitar divisão por zero no passo 2.
    gamma : float
        Limiar para o filtro acumulativo (0 < gamma <= 1). Ex: 0.9 mantém 90% da "massa" dos pesos.
        
    Retorna:
    S_filtered : scipy.sparse.csr_matrix
        Matriz de strong connections simetrizada e filtrada.
    """
    # Garantir formato CSR para acesso eficiente por linhas
    A = sp.csr_matrix(A)
    n_rows, n_cols = A.shape
    
    # ---------------------------------------------------------
    # 1. Obter o stencil MPFA (elementos não nulos)
    # ---------------------------------------------------------
    rows, cols = A.nonzero()
    
    # Filtrar diagonal (j != i), pois o algoritmo define N_i com j != i
    off_diag_mask = rows != cols
    r = rows[off_diag_mask]
    c = cols[off_diag_mask]
    
    # Valores absolutos dos coeficientes
    abs_A_ij = np.abs(A.data[off_diag_mask])
    abs_diag = np.abs(A.diagonal())
    
    # ---------------------------------------------------------
    # 2. Calcular os pesos b_ij
    # ---------------------------------------------------------
    # b_ij = |a_ij| / (sqrt(|a_ii| * |a_jj|) + epsilon)
    denom = np.sqrt(abs_diag[r] * abs_diag[c]) + epsilon
    b_ij = abs_A_ij / denom
    
    # Criar matriz esparsa B apenas com os off-diagonais (pesos)
    B = sp.csr_matrix((b_ij, (r, c)), shape=(n_rows, n_cols))
    
    # ---------------------------------------------------------
    # 3. Determinar o maior peso m_i
    # ---------------------------------------------------------
    # Max por linha. B.max(axis=1) retorna matriz coluna, convertemos para array 1D
    m_i = np.array(B.max(axis=1)).flatten()
    
    # Evitar divisão por zero se uma linha tiver apenas zeros (célula isolada)
    m_i[m_i == 0] = 1.0 
    
    # ---------------------------------------------------------
    # 4. Calcular a força dirigida s_{i->j}
    # ---------------------------------------------------------
    # s_{i->j} = b_ij / m_i
    # Multiplicamos B pela inversa da diagonal de m_i
    D_inv = sp.diags(1.0 / m_i)
    S_directed = D_inv @ B 
    
    # ---------------------------------------------------------
    # 5. Fazer a simetrização S_ij = max[s_{i->j}, s_{j->i}]
    # ---------------------------------------------------------
    S_sym = S_directed.maximum(S_directed.T)
    
    # ---------------------------------------------------------
    # 6. Aplicar o filtro acumulativo
    # ---------------------------------------------------------
    # "Mantenha os maiores b_ij até atingir sum(b_mantidos)/sum(b_total) >= gamma"
    # Usamos B (pesos originais) para decidir o corte, e aplicamos a máscara em S_sym.
    
    S_sym_csr = S_sym.tocsr()
    
    S_filtered_data = []
    S_filtered_indices = []
    S_filtered_indptr = [0]
    
    for i in range(n_rows):
        # Dados da linha i na matriz de pesos B
        start_b = B.indptr[i]
        end_b = B.indptr[i+1]
        
        b_vals = B.data[start_b:end_b]
        b_cols = B.indices[start_b:end_b]
        
        if len(b_vals) == 0:
            S_filtered_indptr.append(len(S_filtered_data))
            continue
            
        # Ordenar b_vals em ordem decrescente para pegar os maiores primeiro
        sort_idx = np.argsort(b_vals)[::-1]
        sorted_vals = b_vals[sort_idx]
        sorted_cols = b_cols[sort_idx]
        
        total_sum = np.sum(sorted_vals)
        if total_sum == 0:
            S_filtered_indptr.append(len(S_filtered_data))
            continue
            
        cumsum = np.cumsum(sorted_vals)
        threshold = gamma * total_sum
        
        # Encontrar quantos elementos manter para atingir o threshold
        # searchsorted retorna o índice onde threshold seria inserido.
        # Queremos manter os elementos até que a soma acumulada >= threshold.
        keep_count = np.searchsorted(cumsum, threshold, side='right') + 1
        keep_count = min(keep_count, len(sorted_cols)) # Garantir que não estoure o array
        
        cols_to_keep = set(sorted_cols[:keep_count])
        
        # Filtrar a matriz simetrizada S_sym para a linha i
        start_s = S_sym_csr.indptr[i]
        end_s = S_sym_csr.indptr[i+1]
        
        s_vals = S_sym_csr.data[start_s:end_s]
        s_cols = S_sym_csr.indices[start_s:end_s]
        
        # Manter apenas as colunas que foram selecionadas pelo filtro de B
        # (Nota: geralmente não filtramos a diagonal se ela existir em S, mas S vem de off-diag)
        mask = np.array([c_idx in cols_to_keep for c_idx in s_cols])
        
        S_filtered_data.extend(s_vals[mask])
        S_filtered_indices.extend(s_cols[mask])
        S_filtered_indptr.append(len(S_filtered_data))
        
    S_filtered = sp.csr_matrix((S_filtered_data, S_filtered_indices, S_filtered_indptr), shape=(n_rows, n_cols))
    
    return S_filtered



    

@time_func
def get_msrsb_prolongation_operator(primal_id: np.ndarray, support_regions: np.ndarray, boundary_regions: np.ndarray, A: sp.csc_matrix, OP0: sp.csc_matrix, tol_op: float=0.05, maxit: int=1000, omega: float=2/3, debug: bool=False, **kwargs) -> Tuple[sp.csc_matrix, int]:
    
    count_it = 0
    Ematrix = _get_Ematrix(A, omega)
    emax = 1e4
    OP = OP0.copy()

    all_boundary_regions = np.unique(np.concatenate(boundary_regions))

    coarse_ids_boundary_regions = np.array([np.repeat(i, len(boundary_regions[i])) for i in range(boundary_regions.shape[0])], dtype='O')
    coarse_ids_support_regions = mount_indices_support_regions(
        primal_id,
        support_regions
    )

    while count_it < maxit and emax > tol_op:
        # OP, emax = _mount_local_op_it(
        #     support_regions,
        #     boundary_regions,
        #     all_boundary_regions,
        #     coarse_ids_support_regions,
        #     coarse_ids_boundary_regions,
        #     Ematrix,
        #     OP
        # )
        OP, emax = _mount_local_op_it_strong(
            support_regions,
            boundary_regions,
            all_boundary_regions,
            coarse_ids_support_regions,
            coarse_ids_boundary_regions,
            Ematrix,
            OP,
            debug=debug
        )
        # OP, emax = _mount_local_op_it_strong_cor(
        #     support_regions,
        #     boundary_regions,
        #     all_boundary_regions,
        #     coarse_ids_support_regions,
        #     coarse_ids_boundary_regions,
        #     Ematrix,
        #     OP,
        #     debug=debug
        # )
        count_it += 1
        
        print(count_it)
        print(emax)
        print()  
    
    soma = np.array(OP.sum(axis=1)).flatten()
    soma[:] = 1/soma
    OP.data *= soma[OP.indices]
    
    return OP, count_it

def get_msrsb_prolongation_operator_v3(ind_support: np.ndarray, ptr_support: np.ndarray, A: sp.csr_matrix, OP0: sp.csc_matrix, op_name: str='', op_logs_path: str='', tol_op: float=0.05, op_maxit: int=10, omega: float=2/3, debug: bool=False, **kwargs) -> Tuple[sp.csc_matrix, int]:
    
    maxit = op_maxit
    count_it = 0
    Ematrix = _get_Ematrix(A, omega)
    emax = 1e4
    OP = copy.deepcopy(OP0)
    Msupport = sp.csr_matrix((np.repeat(1.0, ind_support.shape[0]), ind_support, ptr_support), shape=OP0.transpose().shape).transpose()
    log_local = create_operator_logger(op_logs_path, op_name)
    log_local.info(f'Iniciando {inspect.currentframe().f_code.co_name}')
    log_local.info(f'OP_shape:  {OP.shape}')

    while count_it <= maxit and emax > tol_op:
        OP, emax = _mount_local_op_it_strong_v3(
            Msupport,
            Ematrix,
            OP,
            debug=debug
        )
        count_it += 1
        
        log_local.info(f'Iteracao: {count_it}, emax: {emax:.6f}')  
    
    return OP, count_it, emax

def get_msrsb_prolongation_operator_v2(ind_support: np.ndarray, ptr_support: np.ndarray, A: sp.csr_matrix, OP0: sp.csc_matrix, tol_op: float=0.05, maxit: int=1000, omega: float=2/3, debug: bool=False, **kwargs) -> Tuple[sp.csc_matrix, int]:
    
    count_it = 0
    Ematrix = _get_Ematrix(A, omega)
    emax = 1e4
    OP = copy.deepcopy(OP0)

    while count_it < maxit and emax > tol_op:
        OP, emax = _mount_local_op_it_strong_v2(
            ind_support, 
            ptr_support,
            Ematrix,
            OP,
            debug=debug
        )
        count_it += 1
        
        print(count_it)
        print(emax)
        print()  
    
    return OP, count_it      

@time_func
def enhanced_matrix(A: sp.csc_matrix, **kwargs) -> sp.csc_matrix:
    enhanced = Enhanced()
    return enhanced.get_enhanced_matrix(A)
    
@time_func
def algorithimic_monotone_matrix(A: sp.csc_matrix, epsilon_algorithimic: float=0.001, weight_algorithimic: float=1.0, **kwargs) -> sp.csc_matrix:
    algomonotone = AlgorithimicMonotone()
    M = algomonotone.get_monotone_matrix(A, epsilon=epsilon_algorithimic, w=weight_algorithimic)
    return M
    

def identify(xyz: np.ndarray):
    
    # testx = xyz[:, 0] < 0
    # testy = xyz[:, 1] < 0
    # testz = xyz[:, 2] < 0
    # verify  = testx | testy | testz
    
    # test = xyz < 0
    verify = np.any(xyz < 0, axis=1)
    verify = ~verify
    
    # spe_data = preprocess_spe_mesh()
    # save_ijk_spe(spe_data)
    # spe_data = load_spe_data()
    
    ijk = (xyz - 0.5).astype(np.int64)
    return ijk, verify

def filter_stuben_matrix(A: sp.csc_matrix, theta: float=0.2) -> sp.csr_matrix:
    """
        Mantém a conexão entre ij e ji se:
        |A_ij| >= theta * max_k(|A_ik|)  OU  |A_ji| >= theta * max_k(|A_jk|)
    """
    
    Aoff: sp.csr_matrix = A.tocsr(copy=True)
    Aoff.setdiag(0)
    Aoff.data[Aoff.data > 0] = 0
    Aoff.eliminate_zeros()
    Aoff.data[:] = np.absolute(Aoff.data)
    
    rowmax = Aoff.max(axis=1).toarray().flatten()
    rows, cols = Aoff.nonzero()
    
    is_strong = Aoff.data >= (theta * rowmax[rows])
    
    M = sp.csr_matrix((is_strong, (rows, cols)), shape=A.shape, dtype=bool)
    M = M + M.transpose()
    M.data[:] = 1.0
    
    A_filtered = M.multiply(A)
    soma = A_filtered.sum(axis=1).toarray().flatten()
    soma[soma == 0] = -1
    A_filtered.setdiag(-soma)
    
    return A_filtered


def evolution_filter(A: sp.csr_matrix, epsilon=2.0, k=2) -> sp.csr_matrix:
    """
    Retorna uma matriz que contém apenas as conexões 'fortes' 
    segundo o critério de evolucao.
    """
    # 1. Calcula a matriz de força de conexão (SoC)
    # O PyAMG retorna uma matriz onde os valores indicam a força relativa
    soc: sp.csr_matrix = evolution_strength_of_connection(A, epsilon=epsilon, k=k)
    soc.data[:] = 1
    
    A_filtrada: sp.csr_matrix = soc.multiply(A)
    A_filtrada.setdiag(0)
    soma = np.array(A_filtrada.sum(axis=1)).flatten()
    soma[soma == 0] = -1
    A_filtrada.setdiag(-soma)
    A_filtrada.eliminate_zeros()
    
    return A_filtrada
    
    
def get_S_matrix(A: sp.csc_matrix,**kwargs) -> sp.csc_matrix:
    S = copy.deepcopy(A)
    S.setdiag(0)
    soma = np.array(S.sum(axis=1)).flatten()
    soma[soma == 0] = -1
    S.setdiag(-soma)
    S.eliminate_zeros()
    return S
    
    
    
    
    
    
