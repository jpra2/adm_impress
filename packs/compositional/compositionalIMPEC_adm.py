from .pressure_solver import TPFASolver
from .flux_calculation import FOUM, MUSCL
from .update_time import delta_time
import numpy as np
from ..utils import constants as ctes
from ..directories import data_loaded

from packs.compositional.adm_tpfa_compositional_solver import AdmTpfaCompositionalSolver
from packs.utils.test_functions import test_kwargs_keys

class CompositionalFVM:
    
    _kwargs_keys = {
        '__call__': set(['multilevel_data',
                     'multilevel_operators'])
    }

    def __call__(self, M, wells, fprop, delta_t, **kwargs):
        test_kwargs_keys(CompositionalFVM._kwargs_keys['__call__'], kwargs.keys())
        self.update_gravity_term(fprop)
        if ctes.MUSCL: self.get_faces_properties_average(fprop)
        else: self.get_faces_properties_upwind(fprop)
        self.get_phase_densities_internal_faces(fprop)
        r = 0.8 # enter the while loop
        # psolve = TPFASolver(fprop)
        psolve = AdmTpfaCompositionalSolver(fprop)
        P_old = np.copy(fprop.P)
        Nk_old = np.copy(fprop.Nk)
        while (r!=1.):
            fprop.Nk = np.copy(Nk_old)
            fprop.P, total_flux_internal_faces, self.q = psolve.get_pressure(M, wells, fprop, delta_t,
                                                                             multilevel_data=kwargs.get('multilevel_data'),
                                                                             multilevel_operators=kwargs.get('multilevel_operators'))
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

    def update_gravity_term(self, fprop):
        self.G = ctes.g * fprop.rho_j * ctes.z

    def get_faces_properties_upwind(self, fprop):
        ''' Using one-point upwind approximation '''
        Pot_hid = fprop.P + fprop.Pcap - self.G[0,:,:]
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]

        fprop.mobilities_internal_faces = np.zeros([1, ctes.n_phases, ctes.n_internal_faces])
        mobilities_vols = fprop.mobilities[:,:,ctes.v0[:,0]]
        mobilities_vols_up = fprop.mobilities[:,:,ctes.v0[:,1]]
        fprop.mobilities_internal_faces[0,Pot_hidj_up <= Pot_hidj] = mobilities_vols[0,Pot_hidj_up <= Pot_hidj]
        fprop.mobilities_internal_faces[0,Pot_hidj_up > Pot_hidj] = mobilities_vols_up[0,Pot_hidj_up > Pot_hidj]

        fprop.Csi_j_internal_faces = np.zeros([1, ctes.n_phases, ctes.n_internal_faces])
        Csi_j_vols = fprop.Csi_j[:,:,ctes.v0[:,0]]
        Csi_j_vols_up = fprop.Csi_j[:,:,ctes.v0[:,1]]
        fprop.Csi_j_internal_faces[0,Pot_hidj_up <= Pot_hidj] = Csi_j_vols[0,Pot_hidj_up <= Pot_hidj]
        fprop.Csi_j_internal_faces[0,Pot_hidj_up > Pot_hidj] = Csi_j_vols_up[0,Pot_hidj_up > Pot_hidj]

        fprop.xkj_internal_faces = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_internal_faces])
        xkj_vols = fprop.xkj[:,:,ctes.v0[:,0]]
        xkj_vols_up = fprop.xkj[:,:,ctes.v0[:,1]]
        fprop.xkj_internal_faces[:,Pot_hidj_up <= Pot_hidj] = xkj_vols[:,Pot_hidj_up <= Pot_hidj]
        fprop.xkj_internal_faces[:,Pot_hidj_up > Pot_hidj] = xkj_vols_up[:,Pot_hidj_up > Pot_hidj]


    def get_faces_properties_average(self, fprop):
        fprop.mobilities_internal_faces = (fprop.Vp[ctes.v0[:,0]] * fprop.mobilities[:,:,ctes.v0[:,0]] +
                                                fprop.Vp[ctes.v0[:,1]] * fprop.mobilities[:,:,ctes.v0[:,1]]) /  \
                                                (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])
        fprop.Csi_j_internal_faces = (fprop.Vp[ctes.v0[:,0]] * fprop.Csi_j[:,:,ctes.v0[:,0]] +
                                                fprop.Vp[ctes.v0[:,1]] * fprop.Csi_j[:,:,ctes.v0[:,1]]) /  \
                                                (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])
        fprop.xkj_internal_faces = (fprop.Vp[ctes.v0[:,0]] * fprop.xkj[:,:,ctes.v0[:,0]] +
                                                fprop.Vp[ctes.v0[:,1]] * fprop.xkj[:,:,ctes.v0[:,1]]) /  \
                                                (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])

    def get_phase_densities_internal_faces(self, fprop):
        fprop.rho_j_internal_faces = (fprop.Vp[ctes.v0[:,0]] * fprop.rho_j[:,:,ctes.v0[:,0]] +
                                    fprop.Vp[ctes.v0[:,1]] * fprop.rho_j[:,:,ctes.v0[:,1]]) /  \
                                    (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])

    def update_composition(self, fprop, delta_t):
        fprop.Nk = fprop.Nk + delta_t * (self.q + fprop.Fk_vols_total)
        fprop.z = fprop.Nk[0:ctes.Nc,:] / np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0)

    def update_composition_RK2(self, fprop, Nk_old, delta_t):
        #fprop.Nk = Nk_old + delta_t * (self.q + fprop.Fk_vols_total)
        fprop.Nk = fprop.Nk/2 + Nk_old/2 + 1/2*delta_t * (self.q + fprop.Fk_vols_total)
        fprop.z = fprop.Nk[0:ctes.Nc,:] / np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0)

    def update_composition_RK3_1(self, fprop, Nk_old, delta_t):
        fprop.Nk = 1*fprop.Nk/4 + 3*Nk_old/4 + 1/4*delta_t * (self.q + fprop.Fk_vols_total)
        fprop.z = fprop.Nk[0:ctes.Nc,:] / np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0)

    def update_composition_RK3_2(self, fprop, Nk_old, delta_t):
        fprop.Nk = 2*fprop.Nk/3 + 1*Nk_old/3 + 2/3*delta_t * (self.q + fprop.Fk_vols_total)
        fprop.z = fprop.Nk[0:ctes.Nc,:] / np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0)

        #material balance error calculation:
        #Mb = (np.sum(fprop.Nk - Nk_n,axis=1) - np.sum(self.q,axis=1))/np.sum(self.q,axis=1)
