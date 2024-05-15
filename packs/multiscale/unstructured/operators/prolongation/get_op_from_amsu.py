from packs.multiscale.unstructured.operators.prolongation.dual_interaction_region import DualInteractionRegion
from packs.multiscale.unstructured.operators.prolongation.amsu import AmsU

from typing import Sequence
import scipy.sparse as sp
import numpy as np

def update_global_op_from_amsu(list_dual_interaction_region: Sequence[DualInteractionRegion], OP_AMSU: sp.csc_matrix):
    amsu = AmsU()

    for dual_i in list_dual_interaction_region:
        local_op = amsu.get_local_op(
            dual_i.local_map,
            dual_i.local_vertice,
            dual_i.local_internal_path,
            dual_i.local_boundary,
            dual_i.local_initial_cc,
            dual_i.local_transmissibility,
            dual_i.local_diagonal_term
        )

        OP_AMSU[dual_i.region, dual_i.coarse_id] = local_op

    soma = np.array(OP_AMSU.sum(axis=1)).flatten()

    all_data = sp.find(OP_AMSU)
    lines = all_data[0]
    cols = all_data[1]
    data = all_data[2]
    soma2 = soma[lines]
    data = data/soma2
    OP_AMSU[lines, cols] = data
        
