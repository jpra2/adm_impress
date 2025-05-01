import scipy.sparse as sp
import numpy as np

def define_fine_levels_from_alpha(
    OR: sp.csc_matrix,
    OP: sp.csc_matrix,
    T: sp.csc_matrix,
    primal_ids: np.ndarray,
    alpha_lim: float=1.0,
    **kwargs
) -> np.ndarray:
    
    n = T.shape[0]
    tau_coarse = get_tau_coarse(OR, OP, T)
    data_T_alpha = get_data_T_alpha(T, OP)
    alphai = np.zeros(n)
    
    
    for i in range(n):
        test = data_T_alpha[0] == i
        datai = data_T_alpha[2][test]
        colsi = data_T_alpha[1][test]

        cid = primal_ids[i]
        
        # ##########
        # test2 = datai > 0
        # if test2.sum() == 0:
        #     continue
        # datai = datai[test2]
        # colsi = colsi[test2]
        # ############

        test2 = colsi != cid
        if test2.sum() > 0:
            datai2 = datai[test2]
            colsi2 = colsi[test2]

            max_arg = np.argmax(datai2)
            alphai[i] = datai2[max_arg]/tau_coarse[colsi2[max_arg]]
    
    # test_neg = data_T_alpha[2] < 0
    # test_diag = data_T_alpha[0] == data_T_alpha[1]
    # test_offdiag = ~test_diag
    # test_neg = test_offdiag
    
    # test = test_offdiag & test_neg
    # values_neg = data_T_alpha[2][test]
    # lines_neg = data_T_alpha[0][test]
    # cols_neg = data_T_alpha[1][test]
    
    # alpha_for_test = np.absolute(values_neg)/tau_coarse[cols_neg]
    # alpha_test = alpha_for_test >= alpha_lim
    
    test = alphai > alpha_lim
    return np.arange(n)[test]

def get_tau_coarse(
    OR: sp.csc_matrix,
    OP: sp.csc_matrix,
    T: sp.csc_matrix,
):
    
    T_coarse: sp.csc_matrix = OR@T@OP
    tau_coarse = np.array(T_coarse.diagonal()).flatten()
    return tau_coarse

def get_data_T_alpha(T: sp.csc_matrix, OP: sp.csc_matrix):
     T_alpha = T@OP
     data_T_alpha = sp.find(T_alpha)
     return data_T_alpha
     
def get_alpha_lim_finescale(T: sp.csc_matrix):
    data = sp.find(T)
    diagonal = np.array(T.diagonal()).flatten()
    off_diagonal = data[0] != data[1]
    pos_test = data[2] > 0
    alphas = data[2]/diagonal[data[0]]
    alphas = np.absolute(alphas[(off_diagonal) & (pos_test)])
    return alphas.max()
    
    
      
    
    
    
    
    
    
    