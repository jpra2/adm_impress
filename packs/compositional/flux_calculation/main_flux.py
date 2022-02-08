import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes


if data_loaded['compositional_data']['component_data']['constant_K']:
    from packs.compositional.Kflash import StabilityCheck
else:
    from packs.compositional.stability_check import StabilityCheck


def compute_flux(M, fprop, wells, ft_internal, P_old, Nk_old, Nk_SP_old, Pot_hid, \
    delta_t, t, G):
    if ctes.MUSCL:
        from .MUSCL import MUSCL
        wave_velocity, Fk_vols_total = MUSCL().run(M, fprop, wells, P_old, \
            ft_internal, Pot_hid) #trocar ordem da saida
    elif ctes.FR:
        from .FR_CPR import FR
        wave_velocity, fprop.Nk, fprop.z, fprop.Nk_SP, Fk_vols_total = FR().run(M, fprop, wells,
            ft_internal, Nk_SP_old, P_old, delta_t, t)
    else:
        from .first_order import FirstOrder
        Fk_vols_total, wave_velocity = FirstOrder().run(M, fprop, \
            ft_internal, P_old, Nk_old, G)

    return Fk_vols_total, wave_velocity
