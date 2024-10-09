import numpy as np
import scipy.sparse as sp

def define_fine_levels_from_beta(
    fine_faces_adjacencies,
    OP: sp.csc_matrix,
    beta_lim: float=3.0
) -> np.ndarray:
    
    data_OP = sp.find(OP)
    beta = (1 - data_OP[2])/data_OP[2]
    beta[data_OP[2] < 0] = np.inf
    
    test = beta > beta_lim
    fine_level_ids = np.unique(data_OP[0][test])
    fine_level_ids = np.unique(fine_faces_adjacencies[fine_level_ids])
    return fine_level_ids
    
     
    
    