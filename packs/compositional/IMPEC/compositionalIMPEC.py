from .pressure_solver import TPFASolver
from ..flux_calculation.main_flux import compute_flux
from ..update_time import delta_time
import numpy as np
from packs.utils import constants as ctes
from packs.directories import data_loaded
from .composition_solver import Euler, RK3
import math
import time

class CompositionalFVM:

    def __call__(self, M, wells, fprop, delta_t, t):
        G = self.update_gravity_term(fprop)

        Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]

        #self.get_faces_properties_upwind(fprop, G)
        '''if ctes.n_points>1:
            self.get_faces_properties_weighted_average(fprop, G)
        else: self.get_faces_properties_upwind(fprop, G)'''
        self.get_faces_properties_weighted_average(fprop, G)
        #self.get_faces_properties_harmonic_average(fprop, G)
        #self.get_faces_properties_upwind(fprop, G)
        self.get_phase_densities_internal_faces(fprop, G)


        r = 0.8 # enter the while loop

        dVjdNk, dVjdP = self.dVt_derivatives(fprop)
        psolve = TPFASolver(dVjdNk, dVjdP)

        P_old = np.copy(fprop.P)
        if ctes.FR: Nk_old = np.copy(fprop.Nk_SP)
        else: Nk_old = np.copy(fprop.Nk)

        while (r!=1.):
            #fprop.Nk = np.copy(Nk_old)
            fprop.P, Ft_internal, fprop.qk_molar = psolve.get_pressure(M, wells,
                fprop, P_old, delta_t)

            if any(fprop.P<0): import pdb; pdb.set_trace()

            fprop.qk_prod = fprop.qk_molar[:,wells['ws_prod']]
            #if fprop.qk_prod.sum()==0: import pdb; pdb.set_trace()

            Fk_vols_total, wave_velocity = compute_flux(M, fprop, wells, Ft_internal, \
                P_old, Nk_old, Pot_hid, delta_t, t, G)

            ''' For the composition calculation the time step might be different\
             because it treats composition explicitly and this explicit models \
             are conditionally stable - which can be based on the CFL parameter '''

            #delta_t_new = delta_time.update_CFL(delta_t, fprop, wells, Fk_vols_total, fprop.Nk, wave_velocity)
            #if any(Ft_internal[0]<0): delta_t_new=delta_t/2
            #r = delta_t_new/delta_t
            #delta_t = delta_t_new
            r=1

        if not ctes.FR:
            time_integration = getattr(self, ctes.time_integration)
            time_integration(M, fprop, wells, P_old, Pot_hid, G, Nk_old, \
                Ft_internal, Fk_vols_total, delta_t, t)

        #if fprop.z[-1,0] > 0: import pdb; pdb.set_trace()
        if any(fprop.Sg<0): import pdb; pdb.set_trace()
        fprop.wave_velocity = wave_velocity
        fprop.Ft_internal = Ft_internal
        #fprop.Nk[(fprop.Nk<0)*(abs(fprop.Nk)<1e-300)] = 0
        #if fprop.P[0]<fprop.P[1]: import pdb; pdb.set_trace()
        #fprop.Nk[abs(fprop.Nk)<1e-300] = abs(fprop.Nk[abs(fprop.Nk)<1e-300])
        if any(fprop.Nk.flatten()<0): import pdb; pdb.set_trace()
        #if any(fprop.Sw>1): import pdb; pdb.set_trace()
        #fprop.Nk[fprop.Nk<0] = 1e-30
        #if any(Ft_internal[0]<0): import pdb; pdb.set_trace()
        if any(np.isnan(fprop.Nk).flatten()): import pdb; pdb.set_trace()
        #if any(Ft_internal.flatten()<-1e-6): import pdb; pdb.set_trace()
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
            if ctes.miscible_w:
                dVjdP[0,0,:], dVjdP[0,1,:], dVjdP[0,2,:], \
                dVjdNk[:,0,:], dVjdNk[:,1,:], dVjdNk[:,2,:] = \
                self.EOS.get_all_derivatives(fprop)
            else:
                dVjdP[0,0,:], dVjdP[0,1,:], dVjdNk[0:ctes.Nc,0,:], dVjdNk[0:ctes.Nc,1,:] = \
                self.EOS.get_all_derivatives(fprop)

        if ctes.load_w and not ctes.miscible_w:
            dVjdNk[ctes.n_components-1,-1,:] = 1 / fprop.Csi_j[0,-1,:]
            dVjdP[0,-1,:] = - fprop.Nk[ctes.n_components-1,:] * fprop.Csi_W0 * ctes.Cw / (fprop.Csi_W)**2

        if data_loaded['compositional_data']['component_data']['constant_K']:
            dVjdNk[:,:,:] = 1/fprop.Csi_j
            dVjdP[0,:,:] = 0

        return dVjdNk, dVjdP

    def harmonic_average(self, Vl, prop):
        prop_face = 2 * (prop[:,:,ctes.v0[:,0]] * prop[:,:,ctes.v0[:,1]]) /  \
            (prop[:,:,ctes.v0[:,0]] + prop[:,:,ctes.v0[:,1]])
        prop_face[(prop[:,:,ctes.v0[:,0]]==0) * (prop[:,:,ctes.v0[:,1]]==0)] = 0
        prop_upw = prop[...,ctes.v0[:,0]]
        prop_face[:,Vl[0,:,ctes.v0[:,1]].T==0] = prop_upw[:,Vl[0,:,ctes.v0[:,1]].T==0]
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
        prop_upw = prop[...,ctes.v0[:,0]]#.flatten()
        prop_face[:,fprop.Vj[0,:,ctes.v0[:,1]].T==0] = prop_upw[:,fprop.Vj[0,:,ctes.v0[:,1]].T==0]
        return prop_face

    def get_faces_properties_upwind(self, fprop, G):
        ''' Using one-point upwind approximation '''

        fprop.mobilities_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.mobilities)
        fprop.Csi_j_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.Csi_j)
        fprop.xkj_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.xkj)

        '''fprop.mobilities_internal_faces = fprop.mobilities[...,ctes.v0[:,0]]
        fprop.Csi_j_internal_faces = fprop.Csi_j[...,ctes.v0[:,0]]
        fprop.xkj_internal_faces = fprop.xkj[...,ctes.v0[:,0]]'''


    def get_faces_properties_harmonic_average(self, fprop, G):
        fprop.mobilities_internal_faces = self.harmonic_average(fprop.Vj, fprop.mobilities)
        #fprop.mobilities_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.mobilities)
        fprop.Csi_j_internal_faces = self.harmonic_average(fprop.Vj, fprop.Csi_j)
        fprop.xkj_internal_faces = self.harmonic_average(fprop.Vj, fprop.xkj)

        '''fprop.mobilities_internal_faces[...,0] = fprop.mobilities[...,0]
        fprop.Csi_j_internal_faces[...,0] = fprop.Csi_j[...,0]
        fprop.xkj_internal_faces[...,0] = fprop.xkj[...,0]'''

        '''fprop.mobilities_internal_faces[...,1] = fprop.mobilities[...,-2]
        fprop.Csi_j_internal_faces[...,1] = fprop.Csi_j[...,-2]
        fprop.xkj_internal_faces[...,1] = fprop.xkj[...,-2]'''

    def get_faces_properties_weighted_average(self, fprop, G):
        #fprop.mobilities_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.mobilities)
        fprop.mobilities_internal_faces = self.weighted_by_volume_average(fprop, fprop.mobilities)
        fprop.Csi_j_internal_faces = self.weighted_by_volume_average(fprop, fprop.Csi_j)
        fprop.xkj_internal_faces = self.weighted_by_volume_average(fprop, fprop.xkj)
        #import pdb; pdb.set_trace()
        #upwind no contorno - ajeitar para 2D ... no idea ...
        fprop.mobilities_internal_faces[...,0] = fprop.mobilities[...,0]
        fprop.Csi_j_internal_faces[...,0] = fprop.Csi_j[...,0]
        fprop.xkj_internal_faces[...,0] = fprop.xkj[...,0]

        fprop.mobilities_internal_faces[...,1] = fprop.mobilities[...,-1]
        fprop.Csi_j_internal_faces[...,1] = fprop.Csi_j[...,-1]
        fprop.xkj_internal_faces[...,1] = fprop.xkj[...,-1]

    def get_phase_densities_internal_faces(self, fprop, G):
        #fprop.rho_j_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.rho_j)
        fprop.rho_j_internal_faces = self.weighted_by_volume_average(fprop, fprop.rho_j)

        #material balance error calculation:
        #Mb = (np.sum(fprop.Nk - Nk_n,axis=1) - np.sum(self.q,axis=1))/np.sum(self.q,axis=1)

    def Euler(self, M, fprop, wells, P_old, Pot_hid, G, Nk_old, Ft_internal, Fk_vols_total, delta_t, t):
        fprop.Nk, fprop.z = Euler().update_composition(np.copy(Nk_old), fprop.qk_molar,
            Fk_vols_total, delta_t)

    def RK3(self, M, fprop, wells, P_old, Pot_hid, G, Nk_old, Ft_internal, Fk_vols_total, delta_t, t):
        Nk1 = RK3.update_composition_RK3_1(np.copy(Nk_old), fprop.qk_molar,
            Fk_vols_total, delta_t)
        fprop.Nk = Nk1
        fprop.z = Nk1[:ctes.Nc]/np.sum(Nk1[:ctes.Nc],axis=0)
        #Perform flash calc
        #update properties
        Fk_vols_total, wave_velocity = compute_flux(M, fprop, wells, Ft_internal, \
            P_old, Nk1, Pot_hid, delta_t, t, G)
        Nk2 = RK3.update_composition_RK3_2(np.copy(Nk_old), fprop.qk_molar,
            np.copy(Nk1), Fk_vols_total, delta_t)
        fprop.Nk = Nk2
        fprop.z = Nk2[:ctes.Nc]/np.sum(Nk2[:ctes.Nc],axis=0)
        Fk_vols_total, wave_velocity = compute_flux(M, fprop, wells, Ft_internal, \
            P_old, Nk2, Pot_hid, delta_t, t, G)
        fprop.Nk = RK3.update_composition_RK3_3(np.copy(Nk_old), fprop.qk_molar,
            np.copy(Nk2),Fk_vols_total, delta_t)
        fprop.z = fprop.Nk[:ctes.Nc]/np.sum(fprop.Nk[:ctes.Nc],axis=0)
