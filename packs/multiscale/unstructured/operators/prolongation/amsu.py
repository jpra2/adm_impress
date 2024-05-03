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
        try:
            A = utils_old.get_local_matrix(internal_path, local_transmissibility, local_diagonal_term)
        except TypeError:
            import pdb; pdb.set_trace()
        A[initial_cc, :] = 0
        A[initial_cc, initial_cc] = 1
        A.eliminate_zeros()
        b1 = np.zeros(internal_path.shape[0])
        b1[initial_cc] = 0
        b1[vertice] = 1

        resp = spsolve(A, b1)
        return resp
    
    @staticmethod
    def get_faces_solution(region, boundary, local_transmissibility: sp.csc_matrix, local_diagonal_term, internal_path, initial_boundary_conditions) -> np.ndarray:

        A: sp.csc_matrix = local_transmissibility.copy()
        A.setdiag(A.diagonal() + local_diagonal_term)
        A[internal_path, :] = 0
        A[internal_path, internal_path] = 1
        A.eliminate_zeros()

        b1 = np.zeros(region.shape[0])
        b1[boundary] = 0
        b1[internal_path] = initial_boundary_conditions

        resp = spsolve(A, b1)
        return resp



