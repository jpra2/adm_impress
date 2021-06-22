from packs.utils import constants as ctes
import numpy as np


class Euler:
    @staticmethod
    def update_composition(Nk, q, Fk_vols_total, delta_t):
        Nk = Nk + delta_t * (q + Fk_vols_total)
        z = Nk[0:ctes.Nc,:] / np.sum(Nk[0:ctes.Nc,:], axis = 0)
        return Nk, z

class RK3:
    @staticmethod
    def update_composition_RK3_1(Nk, q, Fk_vols_total, delta_t):
        Nk = Nk + delta_t * (q + Fk_vols_total)
        return Nk

    @staticmethod
    def update_composition_RK3_2(Nk_old, q, Nk, Fk_vols_total, delta_t):
        Nk = 1*Nk/4 + 3*Nk_old/4 + 1/4*delta_t * (q + Fk_vols_total)
        return Nk
        
    @staticmethod
    def update_composition_RK3_3(Nk_old, q, Nk, Fk_vols_total, delta_t):
        Nk = 2*Nk/3 + 1*Nk_old/3 + 2/3*delta_t * (q + Fk_vols_total)
        return Nk
