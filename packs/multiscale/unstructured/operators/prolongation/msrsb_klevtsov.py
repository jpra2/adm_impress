import scipy.sparse as sp
import numpy as np
from typing import Sequence, Tuple
from packs.utils.utils_old import get_local_matrix

class MsRSB:
    '''
    Klevtsov, Sergey
    MULTISCALE SOLVER FOR SUBSURFACE POROMECHANICAL PROBLEMS
    '''
    
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

        ###### so para testar
        support_regions, replicate_coarse_ids = self.get_support_regions(interation_regions, interation_boundaries, coarse_ids)

        OP_0 = OR_fv.transpose().tocsc()
        Ematrix = self.get_Ematrix(omega, T)
        emax = 1e4

        cont = True
        count = 1
        while emax >= etol and count < maxit and cont:
            for i in range(10):
                OP_0, emax = self._mount_local_op_it(
                    coarse_ids,
                    OP_0, 
                    dual_edges, 
                    dual_faces, 
                    vertices, 
                    Ematrix,
                    support_regions
                )
                count += 1
                print(f'MsRSB Klevtsov -- Loop: {count} and emax: {emax} \n')


        soma_op = 1/(np.array(OP_0.sum(axis=1)).flatten())
        OP_0.data *= soma_op[OP_0.indices]

        return OP_0
        
    def _mount_local_op_it(self, coarse_ids, OP_0: sp.csc_matrix, dual_edges, dual_faces, vertices, Ematrix: sp.csc_matrix, support_regions):
        
        
        new_OP = Ematrix@OP_0
        data_new_OP = sp.find(new_OP)
        data_new = []
        lines_new = []
        cols_new = []
        for i in coarse_ids:
            test1 = np.isin(data_new_OP[0], support_regions[i])
            test2 = data_new_OP[1] == i
            test3 = test1 & test2
            data_new.append(data_new_OP[2][test3])
            lines_new.append(data_new_OP[0][test3])
            cols_new.append(data_new_OP[1][test3])
        
        data_new = np.concatenate(data_new)
        lines_new = np.concatenate(lines_new)
        cols_new = np.concatenate(cols_new)

        new_OP = sp.csc_matrix((data_new, (lines_new, cols_new)), shape=OP_0.shape)        
        soma = np.array(new_OP.sum(axis=1)).flatten()
        soma[:] = 1/soma
        soma[dual_faces] = 1

        new_OP.data *= soma[new_OP.indices]

        op_diff = new_OP - OP_0
        difference = np.absolute(op_diff.data[np.isin(op_diff.indices, dual_faces)]).max()

        return new_OP, difference

        
    def get_D_matrix(self, T: sp.csc_matrix) -> sp.csc_matrix:
        n = T.shape[0]
        D = sp.spdiags(1/T.diagonal(), 0, n, n).tocsc()        
        return D
        
    def get_Ematrix(self, omega, T) -> sp.csc_matrix:
        m1: sp.csc_matrix = -omega*(self.get_D_matrix(T)@T)
        m1.setdiag(m1.diagonal() + 1)
        return m1

    def get_support_regions(self, interation_regions, interation_boundaries, coarse_ids):
        support_regions = []
        replicate_coarse_ids = []
        for i in range(len(interation_regions)):
            int_region = interation_regions[i]
            bound_region = interation_boundaries[i]
            support_region = np.setdiff1d(int_region, bound_region)
            support_regions.append(support_region)
            replicate_coarse_ids.append(np.repeat(coarse_ids[i], support_region.shape[0]))
        
        support_regions = np.array(support_regions, dtype='O')
        replicate_coarse_ids = np.array(replicate_coarse_ids, dtype='O')
        
        return support_regions, replicate_coarse_ids
