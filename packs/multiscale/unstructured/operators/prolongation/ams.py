from packs.manager import SuperArrayManager
from packs import defnames
import numpy as np
import scipy.sparse as sp
from typing import Sequence
from packs.multiscale.operators.prolongation.AMS.ams_mpfa import AMSMpfa
from packs.utils import utils_old

class Unstructured2DAmsOperator(SuperArrayManager):

    dual_volumes_str = defnames.dual_volumes_str
    primal_id_str = defnames.fine_primal_id
    dual_id_str = defnames.fine_dual_id

    local_dual_volumes_str = 'local_dual_volumes'
    local_primal_id_str = 'local_primal_id'
    local_dual_id_str = 'local_dual_id'
    local_fine_map_str = 'local_fine_map'
    local_coarse_map_str = 'local_coarse_map'
    coarse_ids_dual_volumes_str = 'coarse_ids_dual_volumes'

    ##[local_matrices, local_ams] other_data

    def insert_ams_data(self, data: dict):
        f"""Insert the ams data for prolongation calculation.
        
        data is a dict with

            data = dict(
                {self.dual_volumes_str} = fine ids in with coarse dual volumes,
                {self.primal_id_str} = coarse primal id of fine faces,
                {self.dual_id_str} = fine dual ids with is inside packs.defnames.dual_ids
            )   
        """

        self.insert_data(data)
    
    def preprocess_ams_data(self):
        local_dual_volumes = []
        local_primal_ids = []
        local_coarse_maps = []
        local_dual_ids = []

        for dual_volume in self[self.dual_volumes_str]:
            local_dual_volume = np.arange(dual_volume.shape[0])
            local_dual_id = self[self.dual_id_str][dual_volume]
            primal_ids_local = self[self.primal_id_str][dual_volume]
            coarse_ids_dual_volume = np.unique(primal_ids_local)
            local_coarse_ids = np.arange(coarse_ids_dual_volume.shape[0])
            local_coarse_map = coarse_ids_dual_volume
            local_primal_id = np.array([local_coarse_ids[coarse_ids_dual_volume == i][0] for i in primal_ids_local])

            local_dual_volumes.append(local_dual_volume)
            local_primal_ids.append(local_primal_id)
            local_coarse_maps.append(local_coarse_map)
            local_dual_ids.append(local_dual_id)
        
        local_dual_volumes = np.array(local_dual_volumes, dtype='O')
        local_primal_ids = np.array(local_primal_ids, dtype='O')
        local_coarse_maps = np.array(local_coarse_maps, dtype='O')
        local_dual_ids = np.array(local_dual_ids, dtype='O')
        
        self.insert_data(
            {
                self.local_dual_volumes_str: local_dual_volumes,
                self.local_primal_id_str: local_primal_ids,
                self.local_coarse_map_str: local_coarse_maps,
                self.local_dual_id_str: local_dual_ids
            }
        )
        
    def get_local_transmissibility_matrix(self, list_of_volumes: Sequence[np.ndarray], T: sp.csc_matrix, diagonal_term: np.ndarray):
        local_matrices = []
        for local_volumes in list_of_volumes:
            local_matrices.append(self.get_local_matrix(local_volumes, T, diagonal_term))
        
        return local_matrices
    
    def get_local_matrix(self, local_volumes, T: sp.csc_matrix, diagonal_term):
        
        return utils_old.get_local_matrix(local_volumes, T, diagonal_term)

    def get_local_ams_op(self, list_of_volumes: Sequence[np.ndarray], T: sp.csc_matrix, diagonal_term: np.ndarray):

        local_prolongation = []
        for i, local_volumes in enumerate(list_of_volumes):
            local_matrix = self.get_local_matrix(local_volumes, T, diagonal_term)
            
            ams = AMSMpfa(
                self[self.local_dual_volumes_str][i].astype(np.int64),
                self[self.local_primal_id_str][i].astype(np.int64),
                self[self.local_dual_id_str][i].astype(np.int64)
            )
            local_coarse_map = self[self.local_coarse_map_str][i]
            
            op = ams.run(local_matrix)
            all_data = sp.find(op)
            lines = local_volumes[all_data[0]]
            cols = local_coarse_map[all_data[1]]
            data = all_data[2]
            local_prolongation.append([lines, cols, data])
        
        return local_prolongation

    def insert_data_in_global_OP(self, list_of_volumes: Sequence[np.ndarray], T: sp.csc_matrix, diagonal_term: np.ndarray, OP: sp.csc_matrix):
        local_op_data = self.get_local_ams_op(list_of_volumes, T, diagonal_term)
        for local_data in local_op_data:
            OP[local_data[0], local_data[1]] = local_data[2]

    def get_B_matrix(self, T: sp.csc_matrix, epsilon=0.001, w=1, lines_to_modify=np.array([])):
        """
        Algorithimic monotone multiscale
        """

        diagonal = T.diagonal()
        n = diagonal.shape[0]
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
    
    def get_monotone_matrix_v2(self, T: sp.csc_matrix, epsilon=0.001, w=1, lines_to_modify=np.array([])):
        
        diagonal = T.diagonal()
        n = diagonal.shape[0]
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

        for i in to_iterate:
            
            test_line = lines == i
            test = (test_line) & (test_pos) & (test_dif)

            columns = cols[test]

            data[test] = data[test] - w*data[test]
            data[(lines == i) & (cols == i)] = 0
            to_diagonal = -data[test_line].sum()
            data[(lines == i) & (cols == i)] = to_diagonal

            for column in columns:
                test_line_2 = lines == column
                t1 = (lines==column) & (cols==i)
                if t1.sum() > 0:
                    ac_data = w*data[t1]
                    data[t1] = data[t1] - ac_data
                    data[(lines==column) & (cols==column)] = 0
                    data[(lines==column) & (cols==column)] = data[test_line_2].sum()
            
        non_zero_data = np.nonzero(data)[0]

        import pdb; pdb.set_trace()

        lines_B = lines[non_zero_data]
        cols_B = cols[non_zero_data]
        data_B = data[non_zero_data]


        B = sp.csc_matrix((data_B,(lines_B,cols_B)), shape=(n, n))
        return B
        
    def get_finite_volume_restriction_operator(self, fine_ids, primal_ids):

        cids = np.unique(primal_ids)
        data_or = np.ones(fine_ids.shape[0], dtype=np.float64)
        OR = sp.csc_matrix((data_or, (primal_ids, fine_ids)), shape=(cids.shape[0], fine_ids.shape[0]))
        return OR

    def get_global_op(self, coarse_faces: np.ndarray, fine_faces: np.ndarray):
        OP = sp.lil_matrix((fine_faces.shape[0], coarse_faces.shape[0]))
        return OP




