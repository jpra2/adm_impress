import scipy.sparse as sp
import numpy as np
from typing import Sequence, Tuple
from packs.utils.utils_old import get_local_matrix

class MsRSB:
    
    def get_OP(
        self,
        faces: np.ndarray,
        T: sp.csc_matrix,
        diagonal_term: np.ndarray,
        interation_regions: Sequence[np.ndarray],
        interation_boundaries: Sequence[np.ndarray],
        vertices: np.ndarray,
        dual_edges: np.ndarray,
        dual_faces: np.ndarray,
        coarse_ids: np.ndarray,
        OR_fv: sp.csc_matrix,
        omega: float=2/3,
        etol: float=1e-3,
        maxit: int=200,
        **kwargs
    ) -> sp.csc_matrix:
        """Get the prolongation operator from MsRSB method

        Args:
            faces (np.ndarray): all fine faces
            T (sp.csc_matrix): fine transmissibility from flux only
            diagonal_term (np.ndarray): aditional fine diagonal term from T 
            interation_regions (Sequence[np.ndarray]):  dual iteration regions sequence 
            interation_boundaries (Sequence[np.ndarray]): boundaries of interation regions
            vertices (np.ndarray): vertices of interation with solution = 1
            dual_edges (np.ndarray): global dual edges
            dual_faces (np.ndarray): global dual faces
            coarse_ids (np.ndarray): coarse ids array
            OR_fv (sp.csc_matrix): finite volume restriction operator 
            omega (float): relaxation
            etol (float): tolerance for iteration
            maxit: (int): max iteration

        Returns:
            sp.csc_matrix: prolongation operator
        """
        
        OP = OR_fv.transpose().tolil()
        OP_toget = OP.copy()
        Ematrix = -omega*(self.get_D_matrix(T)@T)
        emax = 1e4
        toget_Dij = self.mount_toget_Dij(interation_regions, interation_boundaries, coarse_ids, faces)

        cont = True
        count = 1
        while emax >= etol and count < maxit and cont:
            for i in range(10):
                newemax = self._mount_local_op_it(
                    coarse_ids, 
                    OP, 
                    dual_edges, 
                    dual_faces, 
                    vertices, 
                    Ematrix,
                    toget_Dij,
                )
                count += 1
                if newemax < emax:
                    emax = newemax
                    OP_toget = OP.copy()
                else:
                    
                    print(f'MsRSB -- Loop: {count} and emax: {newemax} \n')
                print(f'MsRSB -- Loop: {count} and emax: {newemax} \n')
        
        print(f'MsRSB Operator With {count} iterations and tol={emax} \n')
        
        OP = OP_toget.copy()
        del OP_toget
        soma_op = np.array(OP.sum(axis=1)).flatten()
        data_OP = sp.find(OP)
        data = data_OP[2]/soma_op[data_OP[0]]
        OP2 = sp.csc_matrix((data, (data_OP[0], data_OP[1])), shape=OP.shape)
        
        return OP2
    
    def mount_toget_Dij(self, interation_regions, interation_boundaries, coarse_ids, faces) -> sp.csc_matrix:
        lines = []
        cols = []

        n = faces.shape[0]
        m = coarse_ids.shape[0]
        
        for i, coarse_id in enumerate(coarse_ids):
            region = interation_regions[i]
            boundary = interation_boundaries[i]
            internals = np.setdiff1d(region, boundary, assume_unique=True)
            
            lines.append(internals)
            cols.append(np.repeat(coarse_id, internals.shape[0]))
        
        lines = np.concatenate(lines)
        cols = np.concatenate(cols)
        data = np.ones_like(lines)
        
        Matrix = sp.csc_matrix((data, (lines, cols)), shape=(n,m))
        return Matrix
            
    def _update_Dij_dual_edges(self, Dij: sp.lil_matrix, dual_edges: np.ndarray, OP: sp.csc_matrix) -> None:
              
        soma_dual_edges = np.array(Dij.sum(axis=1)[dual_edges]) 
        # soma_dual_edges = soma_dual_edges.reshape(soma_dual_edges.shape[0], 1)
        Pij_dual_edges = OP[dual_edges].toarray()
        dij_dual_edges = Dij[dual_edges].toarray()
        dij_dual_edges = (dij_dual_edges - Pij_dual_edges*soma_dual_edges)/(1 + soma_dual_edges)
        Dij[dual_edges] = dij_dual_edges
        
    def _mount_local_op_it(self, coarse_ids, OP: sp.lil_matrix, dual_edges, dual_faces, vertices, Ematrix: sp.csc_matrix, toget_Dij):
        
        
        Dij: sp.lil_matrix = (Ematrix@OP).tolil()
        Dij[vertices, coarse_ids] = 0

        self._update_Dij_dual_edges(Dij, dual_edges, OP.tocsc())
        Dij = Dij.tocsc()
        Dij.eliminate_zeros()
        
        dij_final = Dij.multiply(toget_Dij)
        OP[:] = OP[:] + dij_final
        emax = abs(Dij[dual_faces]).max()
        
        return emax
        
    def get_D_matrix(self, T: sp.csc_matrix) -> sp.csc_matrix:
        n = T.shape[0]
        D = sp.spdiags(1/T.diagonal(), 0, n, n).tocsc()        
        return D
        
    
        
            