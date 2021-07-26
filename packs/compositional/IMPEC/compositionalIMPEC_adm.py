from packs.compositional.IMPEC.flux_calculation import Flux, MUSCL
from packs.compositional.update_time import delta_time
import numpy as np
from packs.utils import constants as ctes
from packs.directories import data_loaded

from packs.compositional.IMPEC.adm_tpfa_compositional_solver import AdmTpfaCompositionalSolver
from packs.compositional.IMPEC.compositionalIMPEC import CompositionalFVM
from .composition_solver import Euler

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

        G = self.update_gravity_term(fprop)
        if ctes.MUSCL: self.get_faces_properties_average(fprop)
        else: self.get_faces_properties_upwind(fprop, G)
        self.get_phase_densities_internal_faces(fprop)
        r = 0.8 # enter the while loop
        # psolve = TPFASolver(fprop)
        dVjdNk, dVjdP = self.dVt_derivatives(fprop)

        # psolve = AdmTpfaCompositionalSolver(fprop)
        psolve = AdmTpfaCompositionalSolver(dVjdNk, dVjdP)
        params.update({
            'dVtdP': psolve.dVtP,
            'dVtdNk': psolve.dVtk
        })
        P_old = np.copy(fprop.P)
        Nk_old = np.copy(fprop.Nk)
        pressures = []
        while (r!=1.):
            print(f'\nr: {r}\n')
            fprop.Nk = np.copy(Nk_old)
            fprop.P, total_flux_internal_faces, q = psolve.get_pressure(M, wells, fprop, delta_t, P_old, **kwargs)
            pressures.append(fprop.P)
            #self.update_composition(fprop, delta_t)
            #wave_velocity = MUSCL().run(M, fprop, wells, P_old, total_flux_internal_faces)
            #self.update_composition_RK3_1(fprop, fprop.Nk, delta_t)
            if ctes.MUSCL:
                order = data_loaded['compositional_data']['MUSCL']['order']
                wave_velocity = MUSCL().run(M, fprop, wells, P_old, total_flux_internal_faces, order)
            else:
                UPW = Flux()
                fprop.Fk_vols_total = UPW.update_flux(M, fprop, total_flux_internal_faces,
                                                      fprop.rho_j_internal_faces, fprop.mobilities_internal_faces)
                wave_velocity = UPW.wave_velocity_upw(M, fprop, fprop.mobilities, fprop.rho_j, fprop.xkj,
                                                      fprop.Csi_j, total_flux_internal_faces)

            ''' For the composition calculation the time step might be different\
             because it treats composition explicitly and this explicit models \
             are conditionally stable - which can be based on the CFL parameter '''

            delta_t_new = delta_time.update_CFL(delta_t, wells, fprop, wave_velocity)
            r = delta_t_new/delta_t
            delta_t = delta_t_new

        ##### remove
        # self.update_composition(fprop, delta_t)
        # # self.update_composition_RK3_2(fprop, fprop.Nk, delta_t)
        # return delta_t
        #####

        ########## new
        if not ctes.FR:

            fprop.Nk, fprop.z = Euler().update_composition(fprop.Nk, q,
                                                           fprop.Fk_vols_total, delta_t)
            # wave_velocity = UPW.wave_velocity_upw(M, fprop, fprop.mobilities, fprop.rho_j, fprop.xkj,
            #    fprop.Csi_j, total_flux_internal_faces)

        else:
            fprop.Nk = Nk;
            fprop.z = z;
            fprop.Nk_SP = Nk_SP

        fprop.wave_velocity = wave_velocity
        fprop.total_flux_internal_faces = total_flux_internal_faces
        if any(fprop.xkj.sum(axis=0).flatten() > 1 + 1e-10): import pdb; pdb.set_trace()
        if len(fprop.Nk[fprop.Nk < 0]) > 0: import pdb; pdb.set_trace()
        # if fprop.P[0]<fprop.P[1] or fprop.P[1]<fprop.P[2]: import pdb; pdb.set_trace()
        # import pdb; pdb.set_trace()
        return delta_t
        import pdb; pdb.set_trace()
        ####### end new
