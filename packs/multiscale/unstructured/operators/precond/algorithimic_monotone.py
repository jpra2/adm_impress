import scipy.sparse as sp
import numpy as np

class AlgorithimicMonotone:

    def get_B_matrix(self, T: sp.csc_matrix, epsilon=0.001, w=1, lines_to_modify=np.array([])):
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

    def get_monotone_matrix(self, T: sp.csc_matrix, epsilon=0.001, w=1, lines_to_modify=np.array([])):
        B = self.get_B_matrix(T, epsilon=epsilon, w=w, lines_to_modify=lines_to_modify)
        return T + B
    pass