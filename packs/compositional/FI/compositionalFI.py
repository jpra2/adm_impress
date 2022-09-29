from .pressure_solver import TPFASolver
from .newton_solver import NewtonSolver
from .flux_calculation import FirstOrder, MUSCL, FR
from ..update_time import delta_time
import numpy as np
from packs.utils import constants as ctes
from packs.directories import data_loaded
from .composition_solver import Euler, RK3
import math
import copy

class CompositionalFVM:

    def __call__(self, M, wells, fprop, delta_t, t, flash, StabilityCheck, p1):
        G = self.update_gravity_term(fprop)
        Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]
        '''if ctes.MUSCL or ctes.FR:
            self.get_faces_properties_weighted_average(fprop, G)
        else: self.get_faces_properties_upwind(fprop, G)'''
        self.get_faces_properties_upwind(fprop, G)
        #self.get_faces_properties_weighted_average(fprop, G)
        self.get_phase_densities_internal_faces(fprop)
        face_properties = self.get_faces_properties_upwind
        #face_properties = self.get_faces_properties_weighted_average
        phase_densities = self.get_phase_densities_internal_faces

        r = 0.8 # enter the while loop

        dVjdNk, dVjdP = self.dVt_derivatives(fprop) # muda -----------------------------------------
        psolve = TPFASolver(dVjdNk, dVjdP) # muda --------------------------------------------------

        P_old = np.copy(fprop.P)
        Nk_old = np.copy(fprop.Nk)
        z_old = np.copy(fprop.z)
        if ctes.FR: Nk_SP_old = np.copy(fprop.Nk_SP)
        while (r!=1.):
            #fprop.Nk = np.copy(Nk_old)
            #import pdb; pdb.set_trace()
            #fprop_aux = copy.deepcopy(fprop)

            #solve = NewtonSolver(fprop_aux)
            #Fk_vols_total, wave_velocity, total_flux_internal_faces = \
                #solve.solver(wells, fprop_aux, delta_t, Nk_old, G, flash, StabilityCheck, p1, M)

            solve = NewtonSolver(fprop)
            Fk_vols_total, wave_velocity, total_flux_internal_faces = \
                solve.solver(wells, fprop, delta_t, Nk_old, G, flash, StabilityCheck, p1, M, face_properties, phase_densities, psolve, dVjdNk, dVjdP, P_old)


            #delta_t_new = delta_time.update_CFL(delta_t, fprop, wells, Fk_vols_total, fprop.Nk, wave_velocity)
            #r = delta_t_new/delta_t
            #delta_t = delta_t_new
            r=1
        #import pdb; pdb.set_trace()
        #dd = q

        """if not ctes.FR:
            import pdb; pdb.set_trace()
            fprop.Nk, fprop.z = Euler().update_composition(fprop.Nk, q,
                Fk_vols_total, delta_t)
        else:
            fprop.Nk = Nk; fprop.z = z; fprop.Nk_SP = Nk_SP"""

        fprop.wave_velocity = wave_velocity
        fprop.total_flux_internal_faces = total_flux_internal_faces

        #if fprop.P[2]<fprop.P[3]: import pdb; pdb.set_trace()
        if fprop.P[0]<fprop.P[1]: import pdb; pdb.set_trace()
        if any(fprop.Nk.flatten()<0): import pdb; pdb.set_trace()
        if any(np.isnan(fprop.Nk).flatten()): import pdb; pdb.set_trace()
        if any(total_flux_internal_faces.flatten()<-1e-6): import pdb; pdb.set_trace()
        #if (Nk_old[0,0]>fprop.Nk[0,0]): import pdb; pdb.set_trace()


        #import pdb; pdb.set_trace()

        return delta_t#, fprop

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
                dVjdP[0,0,:], dVjdP[0,1,:], dVjdNk[0:ctes.Nc,0,:], dVjdNk[0:ctes.Nc,1,:] = \
                self.EOS.get_all_derivatives(fprop)

        if ctes.load_w:
            dVjdNk[ctes.n_components-1,2,:] = 1 / fprop.Csi_j[0,ctes.n_phases-1,:]
            dVjdP[0,2,:] = - fprop.Nk[ctes.n_components-1,:] * fprop.Csi_W0 * ctes.Cw / (fprop.Csi_W)**2

        return dVjdNk, dVjdP

    def harmonic_average(self, prop):
        prop_face = 2 * (prop[:,:,ctes.v0[:,0]] * prop[:,:,ctes.v0[:,1]]) /  \
            (prop[:,:,ctes.v0[:,0]] + prop[:,:,ctes.v0[:,1]])
        prop_face[(prop[:,:,ctes.v0[:,0]]==0) * (prop[:,:,ctes.v0[:,1]]==0)] = 0
        return prop_face

    def upwind(self, P, Pcap, G, prop):
        Pot_hid = P + Pcap - G[0,:,:]
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]

        prop_face = np.zeros([prop.shape[0], prop.shape[1], ctes.n_internal_faces])
        prop_vols = prop[:,:,ctes.v0[:,0]]
        prop_vols_up = prop[:,:,ctes.v0[:,1]]
        prop_face[:,Pot_hidj_up <= Pot_hidj] = prop_vols[:,Pot_hidj_up <= Pot_hidj]
        prop_face[:,Pot_hidj_up > Pot_hidj] = prop_vols_up[:,Pot_hidj_up > Pot_hidj]
        return prop_face

    def weighted_by_volume_average(self, fprop, prop):
        prop_face = (fprop.Vp[ctes.v0[:,0]] * prop[:,:,ctes.v0[:,0]] +
                    fprop.Vp[ctes.v0[:,1]] * prop[:,:,ctes.v0[:,1]]) /  \
                    (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])
        return prop_face

    def get_faces_properties_upwind(self, fprop, G):
        ''' Using one-point upwind approximation '''
        Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]

        fprop.mobilities_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.mobilities)
        fprop.Csi_j_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.Csi_j)
        fprop.xkj_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.xkj)
        'TESTAR'
        '''a = (fprop.Nk - abs(fprop.Nk))/2
        a = a[:,ctes.v0]
        fprop.Csi_j_internal_faces[a>=0] = Csi_j_vols[a>=0]
        fprop.Csi_j_internal_faces[a<0] = Csi_j_vols_up[a<0]'''

    def get_faces_properties_harmonic_average(self, fprop, G):
        #fprop.mobilities_internal_faces = self.harmonic_average(fprop.mobilities)
        fprop.mobilities_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.mobilities)
        fprop.Csi_j_internal_faces = self.harmonic_average(fprop.Csi_j)
        fprop.xkj_internal_faces = self.harmonic_average(fprop.xkj)

    def get_faces_properties_weighted_average(self, fprop, G):
        #fprop.mobilities_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.mobilities)
        fprop.mobilities_internal_faces = self.weighted_by_volume_average(fprop, fprop.mobilities)
        fprop.Csi_j_internal_faces = self.weighted_by_volume_average(fprop, fprop.Csi_j)
        fprop.xkj_internal_faces = self.weighted_by_volume_average(fprop, fprop.xkj)
        #import pdb; pdb.set_trace()
        #upwind no contorno
        fprop.mobilities_internal_faces[...,0] = fprop.mobilities[...,0]
        fprop.Csi_j_internal_faces[...,0] = fprop.Csi_j[...,0]
        fprop.xkj_internal_faces[...,0] = fprop.xkj[...,0]

        fprop.mobilities_internal_faces[...,1] = fprop.mobilities[...,-1]
        fprop.Csi_j_internal_faces[...,1] = fprop.Csi_j[...,-1]
        fprop.xkj_internal_faces[...,1] = fprop.xkj[...,-1]

    def get_phase_densities_internal_faces(self, fprop):
        fprop.rho_j_internal_faces = self.weighted_by_volume_average(fprop, fprop.rho_j)


        #material balance error calculation:
        #Mb = (np.sum(fprop.Nk - Nk_n,axis=1) - np.sum(self.q,axis=1))/np.sum(self.q,axis=1)
