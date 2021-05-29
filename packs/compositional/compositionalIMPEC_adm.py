from .pressure_solver import TPFASolver
from .flux_calculation import FOUM, MUSCL
from .update_time import delta_time
import numpy as np
from ..utils import constants as ctes
from ..directories import data_loaded

from packs.compositional.adm_tpfa_compositional_solver import AdmTpfaCompositionalSolver
from packs.utils.test_functions import test_kwargs_keys
from packs.compositional.compositionalIMPEC import CompositionalFVM

class CompositionalFvmADM(CompositionalFVM):
    
    _kwargs_keys = {
        '__call__': [
            'multilevel_data',
            'multilevel_operators'
        ]
    }

    def __call__(self, M, wells, fprop, delta_t, **kwargs):
        # test_kwargs_keys(CompositionalFVM._kwargs_keys['__call__'], kwargs.keys())
        # import pdb; pdb.set_trace()
        params=kwargs.get('params')
        self.update_gravity_term(fprop)
        if ctes.MUSCL: self.get_faces_properties_average(fprop)
        else: self.get_faces_properties_upwind(fprop)
        self.get_phase_densities_internal_faces(fprop)
        r = 0.8 # enter the while loop
        # psolve = TPFASolver(fprop)
        psolve = AdmTpfaCompositionalSolver(fprop)
        # params['dVtdP'] = AdmTpfaCompositionalSolver.dVtP
        # params['dVtdNk'] = AdmTpfaCompositionalSolver.dVtk
        P_old = np.copy(fprop.P)
        Nk_old = np.copy(fprop.Nk)
        while (r!=1.):
            fprop.Nk = np.copy(Nk_old)
            fprop.P, total_flux_internal_faces, self.q = psolve.get_pressure(M, wells, fprop, delta_t,
                                                                             multilevel_data=kwargs.get('multilevel_data'),
                                                                             multilevel_operators=kwargs.get('multilevel_operators'),
                                                                             params=params)
            #self.update_composition(fprop, delta_t)
            #wave_velocity = MUSCL().run(M, fprop, wells, P_old, total_flux_internal_faces)
            #self.update_composition_RK3_1(fprop, fprop.Nk, delta_t)
            if ctes.MUSCL:
                order = data_loaded['compositional_data']['MUSCL']['order']
                wave_velocity = MUSCL().run(M, fprop, wells, P_old, total_flux_internal_faces, order)
            else:
                FOUM().update_flux(fprop, total_flux_internal_faces,
                                     fprop.rho_j_internal_faces,
                                     fprop.mobilities_internal_faces)
                wave_velocity = []

            ''' For the composition calculation the time step might be different\
             because it treats composition explicitly and this explicit models \
             are conditionally stable - which can be based on the CFL parameter '''

            delta_t_new = delta_time.update_CFL(delta_t, wells, fprop, wave_velocity)
            r = delta_t_new/delta_t
            delta_t = delta_t_new

        self.update_composition(fprop, delta_t)
        #self.update_composition_RK3_2(fprop, fprop.Nk, delta_t)
        return delta_t
