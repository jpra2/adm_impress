import numpy as np
import scipy.sparse as sp

class Enhanced:
    
    def get_enhanced_matrix(self, T: sp.csc_matrix):
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        import pdb; pdb.set_trace()
        
