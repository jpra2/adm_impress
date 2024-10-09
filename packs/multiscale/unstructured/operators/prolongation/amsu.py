import numpy as np
import scipy.sparse as sp
from packs.utils import utils_old
from scipy.sparse.linalg import spsolve
from packs import defnames



class AmsU:

    def get_local_op(self, region, vertice, internal_path, boundary, initial_cc, local_transmissibility, local_diagonal_term, local_dual_id) -> np.ndarray:
        
        initial_boundary_conditions, local_edges = self.get_edges_solution(internal_path, vertice, initial_cc, local_transmissibility, local_diagonal_term, local_dual_id, boundary)
        # all_solution = self.get_faces_solution(
        #     region,
        #     boundary,
        #     local_transmissibility,
        #     local_diagonal_term,
        #     internal_path,
        #     initial_boundary_conditions
        # )

        return initial_boundary_conditions, local_edges
    

    def get_local_op_faces(self, region, vertice, internal_path, boundary, initial_cc, local_transmissibility, local_diagonal_term, local_dual_id, local_edges, initial_solution_edges) -> np.ndarray:
        
        all_solution = self.get_faces_solution_v1(
            region,
            boundary,
            local_transmissibility,
            local_diagonal_term,
            internal_path,
            initial_solution_edges,
            local_edges
        )

        return all_solution

    @staticmethod
    def get_edges_solution_v0_dep(internal_path, vertice, initial_cc, local_transmissibility, local_diagonal_term, local_dual_id, boundary) -> np.ndarray:
        
        A = utils_old.get_local_matrix(internal_path, local_transmissibility, local_diagonal_term).tolil()
        local_map = np.arange(internal_path.shape[0])
        local_initial_cc = np.array([local_map[internal_path == i][0] for i in initial_cc])
        local_vertice = local_map[internal_path == vertice][0]

        A[local_initial_cc, :] = 0
        A[local_initial_cc, local_initial_cc] = 1
        A[local_vertice,:] = 0
        A[local_vertice, local_vertice] = 1
        A = A.tocsc()
        A.eliminate_zeros()
        b1 = np.zeros(internal_path.shape[0])
        b1[local_initial_cc] = 0
        b1[local_vertice] = 1

        resp = spsolve(A, b1)
        resp[local_initial_cc] = 0
        resp[local_vertice] = 1
        return resp
    
    @staticmethod
    def get_edges_solution(internal_path, vertice, initial_cc, local_transmissibility, local_diagonal_term, local_dual_id, boundary) -> np.ndarray:
        
        ids_toget_solution = np.unique(np.concatenate([internal_path, boundary]))
        A = utils_old.get_local_matrix(ids_toget_solution, local_transmissibility, local_diagonal_term).tolil()
        local_map = np.arange(ids_toget_solution.shape[0])
        new_initial_cc = np.array([i for i in initial_cc if local_dual_id[i] == defnames.dual_ids('vertice_id')])
        local_initial_cc = np.array([local_map[ids_toget_solution == i][0] for i in new_initial_cc])
        local_vertice = local_map[ids_toget_solution == vertice][0]
        local_boundary = np.array([local_map[ids_toget_solution == i][0] for i in boundary])

        A[local_initial_cc, :] = 0
        A[local_initial_cc, local_initial_cc] = 1
        A[local_vertice,:] = 0
        A[local_vertice, local_vertice] = 1
        A = A.tocsc()
        A.eliminate_zeros()
        b1 = np.zeros(ids_toget_solution.shape[0])
        # b1[local_initial_cc] = 0
        b1[local_vertice] = 1

        resp = spsolve(A, b1)
        resp[local_boundary] = 0
        resp[local_vertice] = 1
        return resp, ids_toget_solution
    
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
    
    @staticmethod
    def get_faces_solution_v1(region, boundary, local_transmissibility: sp.csc_matrix, local_diagonal_term, internal_path, initial_boundary_conditions, local_edges) -> np.ndarray:

        A: sp.csc_matrix = local_transmissibility.copy()
        A.setdiag(A.diagonal() + local_diagonal_term)
        A = A.tolil()
        A[local_edges, :] = 0
        A[local_edges, local_edges] = 1
        A = A.tocsc()
        A.eliminate_zeros()

        b1 = np.zeros(region.shape[0])
        b1[local_edges] = initial_boundary_conditions
        
        resp = spsolve(A, b1)
        
        return resp

    @staticmethod
    def get_Ts_2d(
        fine_transmissibility: sp.csc_matrix,
        permutation: sp.csc_matrix,
        dual_id: np.ndarray,
        *args,
        **kwargs
    ):
        nf = (dual_id == defnames.dual_ids('face_id')).sum()
        ne = (dual_id == defnames.dual_ids('edge_id')).sum()
        nv = (dual_id == defnames.dual_ids('vertice_id')).sum()
        
        permT = permutation.transpose().copy()
        
        Twire = permutation*fine_transmissibility*permT
        
        Tff = Twire[0:nf, 0:nf]
        Tfe = Twire[0:nf, nf:nf+ne]
        Tfv = Twire[0:nf, nf+ne:nf+ne+nv]
        
        Tef = Twire[nf:nf+ne, 0:nf]
        Tee = Twire[nf:nf+ne, nf:nf+ne]
        Tev = Twire[nf:nf+ne, nf+ne:nf+ne+nv]
        
        # Tvf = Twire[nf+ne:nf+ne+nv, 0:nf]
        # Tve = Twire[nf+ne:nf+ne+nv, nf:nf+ne]
        # Tvv = Twire[nf+ne:nf+ne+nv, nf+ne:nf+ne+nv]
        
        resp = {
            'ff': Tff,
            'fe': Tfe,
            'fv': Tfv,
            'ef': Tef,
            'ee': Tee,
            'ev': Tev
        }
        
        return resp, permT, nf, ne, nv
    
    
    @staticmethod 
    def get_As_2d(
        fine_transmissibility: sp.csc_matrix,
        permutation: sp.csc_matrix,
        dual_id: np.ndarray,
        *args,
        **kwargs
    ):
        
        Ts, permT, nf, ne, nv = AmsU.get_Ts_2d(fine_transmissibility, permutation, dual_id)
        Aee: sp.csc_matrix = Ts['ee'].copy()
        
        soma = Ts['ef'].sum(axis=1)
        d1 = np.matrix(Aee.diagonal()).reshape([ne, 1])
        d1 += soma
        Aee.setdiag(d1)
        
        Ts.update({'ee': Aee, 'ef': None})
        
        return Ts, permT, nf, ne, nv
        
    @staticmethod
    def get_Mee_corr_inv(OP_ev, Aev, ne):
        AevT = Aev.transpose().copy()
        Mev: sp.csc_matrix = (Aev*AevT)        
        Mev_inv = spsolve(Mev, sp.identity(ne))
        Mee_cor_inv = (-OP_ev*Aev)*Mev_inv
        return Mee_cor_inv
        
        
        
        
        
        
        
        
        
        
        
    
    @staticmethod
    def get_correction_function(OP_AMSU: sp.csc_matrix, permutation: sp.csc_matrix, fine_source: np.ndarray, fine_transmissibility: sp.csc_matrix, dual_id: np.ndarray):
        As, permT, nf, ne, nv = AmsU.get_As_2d(fine_transmissibility, permutation, dual_id)
        fine_source_perm = permutation*fine_source
        OP_perm = permutation*OP_AMSU
        Mee_corr_inv = AmsU.get_Mee_corr_inv(OP_perm[nf:nf+ne], As['ev'], ne)
        
        pcorr_ee = Mee_corr_inv*fine_source_perm[nf:nf+ne]
        pcorr_fe = -spsolve(As['ff'], As['fe']*pcorr_ee)
        
        pcorr_ff = spsolve(As['ff'], fine_source[0:nf])
        
        pcorr = np.zeros(fine_source.shape[0])
        pcorr[nf:nf+ne] = pcorr_ee
        pcorr[0:nf] = pcorr_ff + pcorr_fe
        
        return permT*pcorr
        
        
        
        
        
        
        
        
        
        