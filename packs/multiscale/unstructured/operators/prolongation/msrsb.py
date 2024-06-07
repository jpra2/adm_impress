import scipy.sparse as sp
import numpy as np
from typing import Sequence, Tuple
from packs.utils.utils_old import get_local_matrix

class MsRSB:
    
    def get_OP(
        self,
        T: sp.csc_matrix,
        diagonal_term: np.ndarray,
        interation_regions: Sequence[np.ndarray],
        interation_boundaries: Sequence[np.ndarray],
        vertices: np.ndarray,
        coarse_ids: np.ndarray,
        OR_fv: sp.csc_matrix,
        **kwargs
    ) -> sp.csc_matrix:
        """Get the prolongation operator from MsRSB method

        Args:
            T (sp.csc_matrix): fine transmissibility from flux only
            diagonal_term (np.ndarray): aditional fine diagonal term from T 
            interation_regions (Sequence[np.ndarray]):  dual iteration regions sequence 
            interation_boundaries (Sequence[np.ndarray]): boundaries of interation regions
            vertices (np.ndarray): vertices of interation with solution = 1
            coarse_ids (np.ndarray): coarse ids array
            OR_fv (sp.csc_matrix): finite volume restriction operator 

        Returns:
            sp.csc_matrix: prolongation operator
        """
        
        OP = OR_fv.transpose().copy()
        Dij = sp.lil_matrix(OP.shape)
        local_interation_region, local_boundaries, local_vertices = self.preprocess(interation_regions, interation_boundaries, vertices)
        local_matrices = self._get_matrices(T, interation_regions, diagonal_term, local_vertices)
        
        
    def _mount_local_op_it(self, local_matrix, local_region, local_boundary, local_vertice_id, local_D, omega):
        pass
        
        
    def get_D_matrices(self, local_matrices: Sequence[sp.csc_matrix]):
        pass
        
      
    def preprocess(self, interation_regions, interation_boundaries, vertices) -> Tuple[Sequence[np.ndarray], Sequence[np.ndarray], np.ndarray]:
        local_maps = []
        local_boundaries = []
        local_vertices = np.arange(vertices.shape[0])
        for i, region in enumerate(interation_regions):
            local_map = np.arange(region.shape[0])
            local_maps.append(local_map)
            local_boundaries.append(np.array([local_map[region==j][0] for j in interation_boundaries[i]]))
            local_vertices[i] = local_map[region==vertices[i]]
        
        local_maps = np.array(local_maps)
        local_boundaries = np.array(local_boundaries)
        
        return local_maps, local_boundaries, local_vertices
    
    
    def _get_matrices(self, T: sp.csc_matrix, dual_volumes: Sequence[np.ndarray], diagonal_term: np.ndarray, local_vertices: np.ndarray) -> Sequence[sp.csc_matrix]:
        matrices = []
        for i, dual_volume in enumerate(dual_volumes):
            T_local = get_local_matrix(dual_volume, T, diagonal_term)
            local_vertice = local_vertices[i]
            T_local[local_vertice,:] = 0
            T_local[local_vertice, local_vertice] = 1.0
            T_local.eliminate_zeros()
            matrices.append(T_local)
        
        return matrices
        
            