from packs.multiscale.unstructured.operators.prolongation.dual_interaction_region import DualInteractionRegion
from packs.multiscale.unstructured.operators.prolongation.amsu import AmsU

from typing import Sequence
import scipy.sparse as sp
import numpy as np

def update_global_op_from_amsu(list_dual_interaction_region: Sequence[DualInteractionRegion], OP_AMSU: sp.csc_matrix):
    amsu = AmsU()
    OP_AMSU = OP_AMSU.tolil()
    all_global_edges = []
    all_local_edges = []

    for dual_i in list_dual_interaction_region:
        
        local_op_edges, local_edges = amsu.get_local_op(
            dual_i.local_map,
            dual_i.local_vertice,
            dual_i.local_internal_path,
            dual_i.local_boundary,
            dual_i.local_initial_cc,
            dual_i.local_transmissibility,
            dual_i.local_diagonal_term,
            dual_i.local_dual_id
        )
        global_edges = dual_i.region[local_edges]
        all_global_edges.append(global_edges)
        all_local_edges.append(local_edges)

        OP_AMSU[global_edges, dual_i.coarse_id] = local_op_edges
    
    OP_AMSU = update_edges_solution(OP_AMSU.tocsc())

    for i, dual_i in enumerate(list_dual_interaction_region):
        local_op_edges = OP_AMSU[all_global_edges[i], dual_i.coarse_id].toarray().flatten()
        local_op_faces = amsu.get_local_op_faces(
            dual_i.local_map,
            dual_i.local_vertice,
            dual_i.local_internal_path,
            dual_i.local_boundary,
            dual_i.local_initial_cc,
            dual_i.local_transmissibility,
            dual_i.local_diagonal_term,
            dual_i.local_dual_id,
            all_local_edges[i],
            local_op_edges
        )

        OP_AMSU[dual_i.region, dual_i.coarse_id] = local_op_faces

    return OP_AMSU.tocsc()

def update_edges_solution(OP_AMSU: sp.csc_matrix) -> sp.lil_matrix:
    OP_AMSU.eliminate_zeros()
    soma = np.array(OP_AMSU.sum(axis=1)).flatten()
    all_data = sp.find(OP_AMSU)
    lines = all_data[0]
    cols = all_data[1]
    data = all_data[2]
    soma2 = soma[lines]
    data = data/soma2
    OP_AMSU = OP_AMSU.tolil()
    OP_AMSU[lines, cols] = data
    return OP_AMSU
