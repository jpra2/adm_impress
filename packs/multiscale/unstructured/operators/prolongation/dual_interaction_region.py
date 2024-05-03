import numpy as np
import scipy.sparse as sp
from typing import Sequence

from packs.utils import utils_old

class DualInteractionRegion:

    def __init__(
            self, 
            region: np.ndarray, 
            vertice: int, 
            coarse_id: int, 
            internal_edge_path_to_vertice: np.ndarray, 
            boundary: np.ndarray, 
            initial_cc: np.ndarray,
            global_transmissibility: sp.csc_matrix,
            global_diagonal_term: np.ndarray
        ):
        
        self.region = region
        self.vertice = vertice
        self.coarse_id = coarse_id
        self.internal_edge_path_to_vertice = internal_edge_path_to_vertice
        self.boundary = boundary
        self.initial_cc = initial_cc
        self.local_diagonal_term = global_diagonal_term[region]
        self.local_transmissibility = self.get_local_transmissibility(region, global_transmissibility, np.zeros(global_diagonal_term.shape[0]))
        self.preprocess_data()

    def get_local_transmissibility(self, local_ids, global_transmissibility, diagonal_term):
        self.local_transmissibility = utils_old.get_local_matrix(local_ids, global_transmissibility, diagonal_term)
    
    def preprocess_data(self):
        self.local_map = np.arange(self.region.shape[0])
        self.local_internal_path = np.array([self.local_map[self.region == i] for i in self.internal_edge_path_to_vertice])
        self.local_boundary = np.array([self.local_map[self.region == i] for i in self.boundary])
        self.local_initial_cc = np.array([self.local_map[self.region == i] for i in self.initial_cc])
        self.local_vertice = self.local_map[self.region == self.vertice][0]


def create_dual_interaction_regions(
        list_regions, 
        list_vertices, 
        list_coarse_id,
        list_internal_path,
        list_boundarys,
        list_initial_cc,
        global_transmissibility,
        global_diagonal_term

    ) -> Sequence[DualInteractionRegion]:
    
    n = len(list_regions)
    all_dual_interaction = []
    for i in range(n):
        all_dual_interaction.append(
            DualInteractionRegion(
                list_regions[i],
                list_vertices[i],
                list_coarse_id[i],
                list_internal_path[i],
                list_boundarys[i],
                list_initial_cc[i],
                global_transmissibility,
                global_diagonal_term
            )
        )
    
    return all_dual_interaction
