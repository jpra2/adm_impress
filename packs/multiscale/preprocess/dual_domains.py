import scipy.sparse as sp
import numpy as np
from packs.multiscale.ms_utils.multiscale_functions import update_local_transmissibility, map_global_id_to_local_id
from packs.multiscale.operators.prolongation.AMS.ams_tpfa import AMSTpfa
from collections.abc import Sequence

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
        self.local_coarse_id = map_global_id_to_local_id(coarse_id)
        
        self.ams_solver = AMSTpfa(
            self.local_ids[self.dual_id == 0],
            self.local_ids[self.dual_id == 1],
            self.local_ids[self.dual_id == 2],
            self.local_ids[self.dual_id == 3],
            self.local_ids,
            self.local_coarse_id
        )
        
        self.As = dict()
        self.local_source_term = np.zeros(len(self.gids))
        
    def update_t_local(self, Tglobal, diagonal_term):
        
        self.Tlocal[:] = update_local_transmissibility(Tglobal, self.gids, diagonal_term)
    
    def update_as(self, Tlocal):
        As = self.ams_solver.get_as(self.ams_solver.get_twire(Tlocal))
        self.As = As

    def update_local_source_term(self, global_source_term):
        self.local_source_term = global_source_term[self.gids]

    def set_update(global_update):
        self.local_update[:] = global_update[self.gids]
    
    def reinitialize_local_update(self):
        self.local_update[:] = False

    

def get_bool_updated_dual_subdomains(dual_subdomains):
    
    updated_dual_subdomains = np.full(len(dual_subdomains), False, dtype=bool)
    
    for i, dual in enumerate(dual_domains):
        if np.any(dual.local_update):
            updated_dual_subdomains[i] = True
    
    return np.array(updated_dual_subdomains)
        

def create_dual_subdomains(dual_volumes, global_flag_dual_id, global_coarse_id):
    
    dual_domains = []
    
    for volumes in dual_volumes:
        dual_domains.append(DualSubdomain(volumes, global_flag_dual_id[volumes], global_coarse_id[volumes]))
    
    dual_domains = np.array(dual_domains)
    
    return dual_domains


def update_local_source_terms_dual_subdomains(dual_subdomains: Sequence[DualSubdomain], global_source_term):
    for dual in dual_subdomains:
        dual.update_local_source_term(global_source_term)