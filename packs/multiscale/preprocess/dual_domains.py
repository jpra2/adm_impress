import scipy.sparse as sp
import numpy as np
from packs.multiscale.ms_utils.multiscale_functions import update_local_transmissibility, map_global_id_to_local_id
from packs.multiscale.operators.prolongation.AMS.ams_tpfa import AMSTpfa
from collections.abc import Sequence
from scipy.sparse.linalg import splu

class DualSubdomain:
    
    def __init__(self, volumes=None, dual_id=None, coarse_id=None, create=True):
        
        # if not create:
        #     return None
    
        n = len(volumes)
        self.gids = volumes
        self.dual_id = dual_id
        self.coarse_id = coarse_id        
        Tlocal = sp.lil_matrix((n, n))
        self.Tlocal = Tlocal.tocsc()
        self.local_update = np.full(n, False, dtype=bool)
        self.local_ids = np.arange(n)
        self.local_coarse_id, self.rmap_lcid_cid = map_global_id_to_local_id(coarse_id)
        rmap_cid_to_lid = np.zeros(self.coarse_id.max() + 1, dtype=int)
        for i, j in zip(self.local_coarse_id, self.rmap_lcid_cid):
             rmap_cid_to_lid[j] = i
        
        self.ams_solver = AMSTpfa(
            self.local_ids[self.dual_id == 0],
            self.local_ids[self.dual_id == 1],
            self.local_ids[self.dual_id == 2],
            self.local_ids[self.dual_id == 3],
            self.local_ids,
            rmap_cid_to_lid[self.coarse_id]
        )
        
        self.As = {
            'Aee': sp.lil_matrix((0, 0)).tocsc(),
            'Aff': sp.lil_matrix((0, 0)).tocsc(),
            'Aii': sp.lil_matrix((0, 0)).tocsc()
        }
        self.lu_matrices = dict()
        self.local_source_term = np.zeros(len(self.gids))
        
    def update_t_local(self, Tglobal, diagonal_term):
        
        self.Tlocal[:] = update_local_transmissibility(Tglobal, self.gids, diagonal_term)
    
    def update_as(self, Tlocal):
        As = self.ams_solver.get_as(self.ams_solver.get_twire(Tlocal))
        self.As.update(As)
        self.update_lu_matrices()

    def update_lu_matrices(self):
        
        self.lu_matrices.update({
            'Aee': splu(self.As['Aee']),
            'Aff': splu(self.As['Aff']),
            'Aii': splu(self.As['Aii'])
        })
    
    def update_local_source_term(self, global_source_term):
        self.local_source_term[:] = global_source_term[self.gids]

    def set_update(self, global_update):
        self.local_update[:] = global_update[self.gids]
    
    def reinitialize_local_update(self):
        self.local_update[:] = False

    def test_update(self):
        if np.any(self.local_update):
            return True
        else:
            return False


class DualSubdomainMethods:
    
    @staticmethod
    def get_bool_update_dual_subdomains(dual_subdomains):
        
        dual_subdomains: Sequence[DualSubdomain]
        updated_dual_subdomains = np.array([dual.test_update() for dual in dual_subdomains], dtype=bool)
        return updated_dual_subdomains

    @staticmethod
    def update_local_source_terms_dual_subdomains(dual_subdomains, global_source_term):
        
        dual_subdomains: Sequence[DualSubdomain]
        for dual in dual_subdomains:
            dual.update_local_source_term(global_source_term)
            
    @staticmethod
    def update_matrices_dual_subdomains(dual_subdomains, Tglobal, global_diagonal_term, test=True):
        
        dual_subdomains: Sequence[DualSubdomain]
        if test:    
            for dual in dual_subdomains:
                if dual.test_update():    
                    dual.update_t_local(Tglobal, global_diagonal_term[dual.gids])
                    dual.update_as(dual.Tlocal)
        else:
            for dual in dual_subdomains:
                dual.update_t_local(Tglobal, global_diagonal_term[dual.gids])
                dual.update_as(dual.Tlocal)
    
    @staticmethod
    def get_subdomains_to_update(dual_subdomains):
        
        dual_subdomains: Sequence[DualSubdomain]
        array_bool = DualSubdomainMethods.get_bool_update_dual_subdomains(dual_subdomains)
        return dual_subdomains[array_bool]

    @staticmethod
    def set_local_update(dual_subdomains, global_update):
        
        for dual in dual_subdomains:
            dual: DualSubdomain
            dual.set_update(global_update)
            
    @staticmethod
    def reinitialize_update(dual_subdomains):
        
        for dual in dual_subdomains:
            dual: DualSubdomain
            dual.reinitialize_local_update()

    @staticmethod
    def reordenate_dual_subdomains(dual_subdomains):
        array_bool = DualSubdomainMethods.get_bool_update_dual_subdomains(dual_subdomains)
        ids_array = np.arange(len(array_bool))
        ids_updated = ids_array[array_bool]
        others = np.setdiff1d(ids_array, ids_updated)
        n_subdomains_not_update = len(others)
        ids_reordenated = np.concatenate([others, ids_updated])
        dual_subdomains_reordenated = dual_subdomains[ids_reordenated]
        return dual_subdomains_reordenated, n_subdomains_not_update
        
        
def create_dual_subdomains(dual_volumes, global_flag_dual_id, global_coarse_id):
        
        dual_domains = []
        
        for volumes in dual_volumes:
            dual_domains.append(DualSubdomain(volumes, global_flag_dual_id[volumes], global_coarse_id[volumes]))
        
        dual_domains = np.array(dual_domains)
        return dual_domains