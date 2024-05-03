from packs.multiscale.unstructured.operators.prolongation.dual_interaction_region import DualInteractionRegion
from packs.multiscale.unstructured.operators.prolongation.amsu import AmsU

from typing import Sequence
import scipy.sparse as sp

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
    import pdb; pdb.set_trace()
    soma = OP_AMSU.sum(axis=1).toarray().flatten()
    for i in range(OP_AMSU.shape[1]):
        OP_AMSU[:, i] = OP_AMSU[:, i].toarray().flatten()/soma
    
    OP_AMSU.eliminate_zeros()
    import pdb; pdb.set_trace()




        
