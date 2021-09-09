from .pressure_solver import TPFASolver
from .flux_calculation import Flux, MUSCL, FR
from ..update_time import delta_time
import numpy as np
from packs.utils import constants as ctes
from packs.directories import data_loaded
from .composition_solver import Euler, RK3
import math

class CompositionalFVM:

    def __call__(self, M, wells, fprop, delta_t, t):
        G = self.update_gravity_term(fprop)

        if ctes.MUSCL or ctes.FR: self.get_faces_properties_average(fprop)
        else: self.get_faces_properties_upwind(fprop, G)
        self.get_phase_densities_internal_faces(fprop)

        r = 0.8 # enter the while loop

        dVjdNk, dVjdP = self.dVt_derivatives(fprop)
        psolve = TPFASolver(dVjdNk, dVjdP)

        P_old = np.copy(fprop.P)
        Nk_old = np.copy(fprop.Nk)
        if ctes.FR: Nk_SP_old = np.copy(fprop.Nk_SP)
        import pdb; pdb.set_trace()
        while (r!=1.):

            fprop.Nk = np.copy(Nk_old)
            fprop.P, total_flux_internal_faces, q = psolve.get_pressure(M, wells, fprop, P_old, delta_t)

            if ctes.MUSCL:
                #order = data_loaded['compositional_data']['MUSCL']['order']
                wave_velocity = MUSCL().run(M, fprop, wells, P_old, total_flux_internal_faces)

            elif ctes.FR:
                wave_velocity, Nk, z, Nk_SP = FR().run(M, fprop, wells,
                    total_flux_internal_faces, Nk_SP_old, P_old, q, delta_t, t)
            else:
                UPW = Flux()
                fprop.Fk_vols_total, wave_velocity = UPW.update_flux(M, fprop, P_old, total_flux_internal_faces,
                                     fprop.rho_j_internal_faces,
                                     fprop.mobilities_internal_faces)



            ''' For the composition calculation the time step might be different\
             because it treats composition explicitly and this explicit models \
             are conditionally stable - which can be based on the CFL parameter '''
            #wave_velocity[:,wells['all_wells']] = 0
            delta_t_new = delta_time.update_CFL(delta_t, fprop.Fk_vols_total, fprop.Nk, wave_velocity)
            r = delta_t_new/delta_t
            delta_t = delta_t_new


        dd = q
        if any(fprop.P<0): import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()
        if not ctes.FR:
            fprop.Nk, fprop.z = Euler().update_composition(fprop.Nk, q,
                fprop.Fk_vols_total, delta_t)

        else:
            fprop.Nk = Nk; fprop.z = z; fprop.Nk_SP = Nk_SP

        fprop.wave_velocity = wave_velocity
        fprop.total_flux_internal_faces = total_flux_internal_faces
        if any(fprop.xkj.sum(axis=0).flatten()>1+1e-10): import pdb; pdb.set_trace()
        if len(fprop.Nk[fprop.Nk<0]>1): import pdb; pdb.set_trace()
        #if fprop.P[0]<fprop.P[1] or fprop.P[1]<fprop.P[2]: import pdb; pdb.set_trace()
        return delta_t

    def update_gravity_term(self, fprop):
        if any((ctes.z - ctes.z[0]) != 0):
            G = ctes.g * fprop.rho_j * ctes.z
        else:
            G = np.zeros_like(fprop.rho_j)
        return G

    def dVt_derivatives(self, fprop):
        dVjdNk = np.zeros((ctes.n_components, ctes.n_phases, ctes.n_volumes))
        dVjdP = np.empty((1, ctes.n_phases, ctes.n_volumes))
        if ctes.load_k:
            self.EOS = ctes.EOS_class(fprop.T)
            if not ctes.compressible_k:
                dVjdNk[0:ctes.Nc,0,:] = 1 / fprop.Csi_j[0,0,:]
                dVjdNk[0:ctes.Nc,1,:] =  np.zeros_like(dVjdNk[0:ctes.Nc,0])#1 / fprop.Csi_j[0,1,:]
                dVjdP[0,0,:] = np.zeros(ctes.n_volumes)
                dVjdP[0,1,:] = np.zeros_like(dVjdP[0,0,:])
            else:
                dVjdP[0,0,:], dVjdP[0,1,:], dVjdNk[0:ctes.Nc,0,:], dVjdNk[0:ctes.Nc,1,:] = self.EOS.get_all_derivatives(fprop)
                #dVtP = dVjdP.sum(axis=0)
        #else:
            #dVjdP[0,:] = np.zeros(ctes.n_volumes); dVjdP[1,:] = np.zeros(ctes.n_volumes)
            #dVjdNk[:,0,:] = np.zeros((ctes.n_components, ctes.n_volumes));
            #dVjdNk[:,1,:] =  np.zeros_like(dVldNk)

        if ctes.load_w:
            dVjdNk[ctes.n_components-1,2,:] = 1 / fprop.Csi_j[0,ctes.n_phases-1,:]
            dVjdP[0,2,:] = - fprop.Nk[ctes.n_components-1,:] * fprop.Csi_W0 * ctes.Cw / (fprop.Csi_W)**2

        #else: dVjdP[2,:] = np.zeros(ctes.n_volumes)
        return dVjdNk, dVjdP

    def get_faces_properties_upwind(self, fprop, G):
        ''' Using one-point upwind approximation '''
        Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]
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

        'TESTAR'
        '''a = (fprop.Nk - abs(fprop.Nk))/2
        a = a[:,ctes.v0]
        fprop.Csi_j_internal_faces[a>=0] = Csi_j_vols[a>=0]
        fprop.Csi_j_internal_faces[a<0] = Csi_j_vols_up[a<0]'''

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


        #material balance error calculation:
        #Mb = (np.sum(fprop.Nk - Nk_n,axis=1) - np.sum(self.q,axis=1))/np.sum(self.q,axis=1)
