import numpy as np
import scipy.sparse as sp
from packs.utils import utils_old
from scipy.sparse.linalg import spsolve

class AmsU:

    def get_local_op(self, region, vertice, internal_path, boundary, initial_cc, local_transmissibility, local_diagonal_term) -> np.ndarray:
        initial_boundary_conditions = self.get_edges_solution(internal_path, vertice, initial_cc, local_transmissibility, local_diagonal_term)
        all_solution = self.get_faces_solution(
            region,
            boundary,
            local_transmissibility,
            local_diagonal_term,
            internal_path,
            initial_boundary_conditions
        )

        return all_solution

    @staticmethod
    def get_edges_solution(internal_path, vertice, initial_cc, local_transmissibility, local_diagonal_term) -> np.ndarray:
        
        A = utils_old.get_local_matrix(internal_path, local_transmissibility, local_diagonal_term)
        local_map = np.arange(internal_path.shape[0])
        local_initial_cc = np.array([local_map[internal_path == i][0] for i in initial_cc])
        local_vertice = local_map[internal_path == vertice][0]

        A[local_initial_cc, :] = 0
        A[local_initial_cc, local_initial_cc] = 1
        A[local_vertice,:] = 0
        A[local_vertice, local_vertice] = 1
        A.eliminate_zeros()
        b1 = np.zeros(internal_path.shape[0])
        b1[local_initial_cc] = 0
        b1[local_vertice] = 1

        resp = spsolve(A, b1)
        resp[local_initial_cc] = 0
        resp[local_vertice] = 1
        return resp
    
    @staticmethod
    def get_faces_solution(region, boundary, local_transmissibility: sp.csc_matrix, local_diagonal_term, internal_path, initial_boundary_conditions) -> np.ndarray:

        A: sp.csc_matrix = local_transmissibility.copy()
        A.setdiag(A.diagonal() + local_diagonal_term)
        A[internal_path, :] = 0
        A[internal_path, internal_path] = 1
        A.eliminate_zeros()

        b1 = np.zeros(region.shape[0])
        b1[internal_path] = initial_boundary_conditions
        b1[boundary] = 0
        

        resp = spsolve(A, b1)
        resp[boundary] = 0
        resp[internal_path] = initial_boundary_conditions
        return resp
    



