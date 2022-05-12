from packs.compositional.IMPEC.flux_calculation import Flux, MUSCL
from packs.compositional.update_time import delta_time
import numpy as np
from packs.utils import constants as ctes
from packs.directories import data_loaded
from packs.compositional.flux_calculation.main_flux import compute_flux

from packs.compositional.IMPEC.adm_tpfa_compositional_solver import AdmTpfaCompositionalSolver
from packs.compositional.IMPEC.compositionalIMPEC import CompositionalFVM
from .composition_solver import Euler, RK3
from .flux_calculation import FirstOrder, MUSCL, FR
import math

class CompositionalFvmADM(CompositionalFVM):
    
    _kwargs_keys = {
        '__call__': [
            'multilevel_data',
            'multilevel_operators'
        ]
    }

    def __call__(self, M, wells, fprop, delta_t, t, params, **kwargs):
        # test_kwargs_keys(CompositionalFVM._kwargs_keys['__call__'], kwargs.keys())
        # import pdb; pdb.set_trace()

        G = self.update_gravity_term(fprop)
        Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]
        '''if ctes.MUSCL or ctes.FR:
            self.get_faces_properties_weighted_average(fprop, G)
        else: self.get_faces_properties_upwind(fprop, G)'''
        self.get_faces_properties_upwind(fprop, G)
        self.get_phase_densities_internal_faces(fprop)

        r = 0.8 # enter the while loop

        dVjdNk, dVjdP = self.dVt_derivatives(fprop)
        
        psolve = AdmTpfaCompositionalSolver(dVjdNk, dVjdP)
        params.update({
            'dVtdP': psolve.dVtP,
            'dVtdk': psolve.dVtk
        })
        P_old = np.copy(fprop.P)
        Nk_old = np.copy(fprop.Nk)
        Nk_SP_old = np.copy(fprop.Nk_SP)
        
        while (r!=1.):
            
            fprop.Nk = np.copy(Nk_old)

            # fprop.P, total_flux_internal_faces, q = psolve.get_pressure(M, wells, fprop, P_old, delta_t, **kwargs)
            
            fprop.P, Ft_internal, fprop.qk_molar = psolve.get_pressure(M, wells,
                fprop, P_old, delta_t, params=params, **kwargs)
            if any(np.isnan(fprop.P)): import pdb; pdb.set_trace()

            '''total_flux_internal_faces = np.ones((1,ctes.n_internal_faces)) * 1/(24*60*60)
            q = np.zeros_like(fprop.Nk)
            frj = fprop.mobilities[:,...] / \
                np.sum(fprop.mobilities[:,...], axis = 1)
            frj[:,1,0] = 1
            frj[:,0,0] = 0
            q[:,wells['all_wells']] = np.sum(frj[:,:,wells['all_wells']] * fprop.Csi_j[:,:,wells['all_wells']]*\
                np.array([[1, 0, 0, 0, 0],[0,0.25, 0.25, 0.25, 0.25]]).T[:,np.newaxis,:] * \
                total_flux_internal_faces[:,0], axis=1)
            q[:,-1] = -1*q[:,-1]
            fprop.q_phase = total_flux_internal_faces[:,0][:,np.newaxis] * np.ones((1,2))
            '''
            # #import pdb; pdb.set_trace()
            # if ctes.MUSCL:
            #     wave_velocity, Fk_vols_total = MUSCL().run(M, fprop, wells, P_old, \
            #         total_flux_internal_faces, Pot_hid)

            # elif ctes.FR:
            #     wave_velocity, Nk, z, Nk_SP, Fk_vols_total = FR().run(M, fprop, wells,
            #         total_flux_internal_faces, Nk_SP_old, P_old, q, delta_t, t)
            # else:
            #     if ctes.RS['LLF']:
            #         Fk_vols_total, wave_velocity = FirstOrder().LLF(M, fprop, total_flux_internal_faces, P_old)
            #     elif ctes.RS['MDW']:
            #         Fk_vols_total, wave_velocity = FirstOrder().MDW(M, fprop, total_flux_internal_faces, P_old)
            #     elif ctes.RS['ROE']:
            #         Fk_vols_total, wave_velocity = FirstOrder().ROE(M, fprop, total_flux_internal_faces, P_old)
            #     else:
            #         self.get_faces_properties_upwind(fprop, G)
            #         Fk_vols_total, wave_velocity = FirstOrder().FOU(M, fprop, total_flux_internal_faces)
            
            Fk_vols_total, wave_velocity = compute_flux(M, fprop, wells, Ft_internal, \
                P_old, Nk_old, Nk_SP_old, Pot_hid, delta_t, t, G)

            ''' For the composition calculation the time step might be different\
             because it treats composition explicitly and this explicit models \
             are conditionally stable - which can be based on the CFL parameter '''

            # delta_t_new = delta_time.update_CFL(delta_t, Fk_vols_total, fprop.Nk, wave_velocity)
            delta_t_new = delta_time.update_CFL(delta_t, fprop, wells, Fk_vols_total, fprop.Nk, wave_velocity)
            r = delta_t_new/delta_t
            # delta_t = delta_t_new
            r = 1
        
            # ########
            # r = 1
            # ########

        # if not ctes.FR:
        #     fprop.Nk, fprop.z = Euler().update_composition(fprop.Nk, q,
        #         Fk_vols_total, delta_t)
        # else:
        #     fprop.Nk = Nk; fprop.z = z; fprop.Nk_SP = Nk_SP
        
        if not ctes.FR:
            #Fk_vols_total *= 1/ctes.ds_faces
            fprop.Nk, fprop.z = Euler().update_composition(fprop.Nk, fprop.qk_molar,
                Fk_vols_total, delta_t)

        fprop.wave_velocity = wave_velocity
        # fprop.total_flux_internal_faces = total_flux_internal_faces
        fprop.Ft_internal = Ft_internal
        #import pdb; pdb.set_trace()
        # if fprop.P[2]<fprop.P[3]: import pdb; pdb.set_trace()
        if any(fprop.Nk.flatten()<0): import pdb; pdb.set_trace()
        if any(np.isnan(fprop.Nk).flatten()): import pdb; pdb.set_trace()
        # if any(total_flux_internal_faces.flatten()<-1e-6): import pdb; pdb.set_trace()
        #if (Nk_old[0,0]>fprop.Nk[0,0]): import pdb; pdb.set_trace()
        return delta_t
