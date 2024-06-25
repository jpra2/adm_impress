import scipy.sparse as sp
import numpy as np

def define_fine_levels_from_alpha(
    OR: sp.csc_matrix,
    OP: sp.csc_matrix,
    T: sp.csc_matrix,
    alpha_lim: float=1.0,
    **kwargs
) -> np.ndarray:
    
    tau_coarse = get_tau_coarse(OR, OP, T)
    data_T_alpha = get_data_T_alpha(T, OP)
    
    test_neg = data_T_alpha[2] < 0
    test_diag = data_T_alpha[0] == data_T_alpha[1]
    test_offdiag = ~test_diag
    
    test = test_offdiag & test_neg
    values_neg = data_T_alpha[2][test]
    lines_neg = data_T_alpha[0][test]
    cols_neg = data_T_alpha[1][test]
    
    alpha_for_test = np.absolute(values_neg)/tau_coarse[cols_neg]
    alpha_test = alpha_for_test >= alpha_lim
    
    return np.unique(lines_neg[alpha_test])

def get_tau_coarse(
    OR: sp.csc_matrix,
    OP: sp.csc_matrix,
    T: sp.csc_matrix,
):
    
    T_coarse: sp.csc_matrix = OR@T@OP
    tau_coarse = T_coarse.diagonal()
    return tau_coarse

def get_data_T_alpha(T: sp.csc_matrix, OP: sp.csc_matrix):
     T_alpha = T@OP
     data_T_alpha = sp.find(T_alpha)
     return data_T_alpha
     
    
    
    
    
    
    
    
    