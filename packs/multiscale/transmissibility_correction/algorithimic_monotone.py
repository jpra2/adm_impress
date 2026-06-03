import scipy.sparse as sp
import numpy as np
import copy

class AlgorithimicMonotone:

    def get_B_matrix(self, T: sp.csc_matrix, epsilon=0.001, w=1.0, lines_to_modify=np.array([])) -> sp.csc_matrix:
        """
        Algorithimic monotone multiscale
        """

        diagonal = T.diagonal()
        n = diagonal.shape[0]

        if w == 0:
            B = sp.lil_matrix((n, n))
            return B.tocsc()


        all_data = sp.find(T)
        lines = all_data[0]
        cols = all_data[1]
        data = all_data[2]

        if lines_to_modify.shape[0] == 0:
            to_iterate = range(n)
        else:
            to_iterate = lines_to_modify

        test_pos = data > 0
        test_dif = lines == cols
        test_dif = ~test_dif
        test_pos_off_diagonal = (test_pos) & (test_dif)
        if test_pos_off_diagonal.sum() > 0:
            pass
        else:
            B = sp.lil_matrix((n, n))
            return B.tocsc()

        lines_B = []
        cols_B = []
        data_B = []

        for i in to_iterate:
            test_line = lines == i
            test = (test_line) & (test_pos) & (test_dif)
            local_data = data[test]
            columns = cols[test]
            xi = local_data/diagonal[i]
            test2 = xi > epsilon
            selected_columns = columns[test2]
            selected_data = local_data[test2]
            for col, ac_data in zip(selected_columns, selected_data):
                data_B.append([-w*ac_data, w*ac_data, -w*ac_data, w*ac_data])
                lines_B.append([i, i, col, col])
                cols_B.append([col, i, i, col])
        
        lines_B = np.concatenate(lines_B)
        cols_B = np.concatenate(cols_B)
        data_B = np.concatenate(data_B)

        B = sp.csc_matrix((data_B,(lines_B,cols_B)), shape=(n, n))
        return B

    def get_complete_B_matrix(self, T: sp.csc_matrix, epsilon=0.001, w=1.0, **kwargs):
        T2: sp.coo_matrix = copy.deepcopy(T).tocoo()
        diag = np.array(T2.diagonal()).flatten()
        n = T2.shape[0]
        
        test_off_diagonal = T2.row != T2.col
        test_pos = T2.data > 0
        test_off_diagonal_pos = test_off_diagonal & test_pos
        
        if test_off_diagonal_pos.sum() > 0:
            pass
        else:
            B = sp.lil_matrix((n, n))
            return B.tocsc()
        
        values = T2.data[test_off_diagonal_pos]
        lines = T2.row[test_off_diagonal_pos]
        cols = T2.col[test_off_diagonal_pos]
        values_test = np.absolute(values/diag[lines])
        test_eps = values_test > epsilon
        
        if test_eps.sum() > 0:
            pass
        else:
            B = sp.lil_matrix((n, n))
            return B.tocsc()
        
        values_eps = w*values[test_eps]
        lines_eps = lines[test_eps]
        cols_eps = cols[test_eps]
        
        lines_B = np.concatenate([lines_eps, lines_eps, cols_eps, cols_eps])
        cols_B = np.concatenate([cols_eps, lines_eps, lines_eps, cols_eps])
        data_B = np.concatenate([-values_eps, values_eps, -values_eps, values_eps])
        
        B = sp.csc_matrix((data_B,(lines_B,cols_B)), shape=(n, n))
        return B
    
    def get_monotone_matrix(self, T: sp.csc_matrix, epsilon=0.001, w=1.0, lines_to_modify=np.array([])) -> sp.csc_matrix:
        # B = self.get_B_matrix(T, epsilon=epsilon, w=w, lines_to_modify=lines_to_modify)
        B = self.get_complete_B_matrix(T, epsilon=epsilon, w=w, lines_to_modify=lines_to_modify)
        return T + B