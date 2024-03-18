from packs.manager import SuperArrayManager
from packs import defnames
import numpy as np
import scipy.sparse as sp
from typing import Sequence

class Unstructured2DAmsProlongation(SuperArrayManager):

    dual_volumes_str = defnames.dual_volumes_str
    primal_id_str = defnames.fine_primal_id
    dual_id_str = defnames.fine_dual_id

    local_dual_volumes_str = 'local_dual_volumes'
    local_primal_id_str = 'local_primal_id'
    local_dual_id_str = 'local_dual_id'
    local_fine_map_str = 'local_fine_map'
    local_coarse_map_str = 'local_coarse_map'
    coarse_ids_dual_volumes_str = 'coarse_ids_dual_volumes'

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
            T2 = T[local_volumes][:,local_volumes].copy()
            data = np.array(T2.sum(axis=1).transpose())[0]
            data2 = T2.diagonal()
            data2 -= data
            T2.setdiag(data2)
            T2[local_volumes, local_volumes] += diagonal_term[local_volumes]
            local_matrices.append(T2)
        
        return local_matrices
    














