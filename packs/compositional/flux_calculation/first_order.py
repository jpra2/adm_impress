from .riemann_solvers import RiemannSolvers
from packs.utils import constants as ctes
from packs.directories import data_loaded
import scipy.sparse as sp
import numpy as np
from .flux_volumes import Flux
import math


class FirstOrder:
    def __init__(self):
        pass

    def run(self,  M, fprop, Ft_internal, P_old, Nk_old, G):
        self.Nk = np.copy(Nk_old)
        if ctes.RS != 'UPW':
            Nk_face = Nk_old[:,ctes.v0]
            P_face = P_old[ctes.v0].sum(axis=-1)/2
            P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)
            Vp_face = fprop.Vp[ctes.v0].sum(axis=-1)/2

            RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)

            Fk_face = RS.get_Fk_face(fprop, M, Nk_face, P_face, Vp_face, Ft_internal)
            solver = getattr(RS, ctes.RS)
            Fk_internal_faces, wave_velocity = solver(M, fprop, Nk_face, P_face,
                Ft_internal, Fk_face, np.ones_like(Ft_internal[0],dtype=bool))
            Fk_internal_faces[:,[0,1]] = Fk_face[:,[0,1],0]

            Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)

        else:
            self.get_faces_properties_upwind(fprop, G)
            Fk_vols_total, wave_velocity = self.FOU(M, fprop, \
                Ft_internal, fprop.mobilities_internal_faces)
        return Fk_vols_total, wave_velocity

    def FOU(self, M, fprop, Ft_internal, mobilities_internal_faces):
        UPW = Flux()

        #BURGERS
        '''Nk_faces = self.Nk[:,ctes.v0]
        Nk_faces_contour = np.copy(Nk_faces[:,np.newaxis,0])
        Nk_faces_contour[:,0,0] = self.Nk[:,-1]
        Nk_faces_contour[:,0,1] = self.Nk[:,0]
        Nk_faces_all = np.concatenate((Nk_faces,Nk_faces_contour),axis=1)

        Fk_faces_all = (Nk_faces_all**2/2)/(1/ctes.n_volumes)
        dF_faces_all = Nk_faces_all/(1/ctes.n_volumes)
        A = (abs(Nk_faces_all[...,0]-Nk_faces_all[...,1])>10**(-10))
        B = ((dF_faces_all[...,0]*dF_faces_all[...,1])<0)
        a = dF_faces_all[...,0]
        a[A] = (Fk_faces_all[A,0] - Fk_faces_all[A,1])/(Nk_faces_all[A,0] - Nk_faces_all[A,1])
        beta = np.max(abs(Fk_faces_all),axis=-1)
        a[B] = beta[B]'''
        #a[a>1] = 1

        #comment for burger

        Fk_vols_total = UPW.update_flux(M, fprop, Ft_internal, fprop.rho_j_internal_faces, \
            mobilities_internal_faces)

        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)
        Nk_face = self.Nk[:,ctes.v0]#.sum(axis=-1)/2
        P_face = fprop.P[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)

        ponteiro = np.ones(ctes.n_internal_faces,dtype=bool)
        #dNkmax_small = np.max(abs(Nk_face[:,:,0]-Nk_face[:,:,1]),axis=0)<1e-20
        #ponteiro[dNkmax_small] = False
        wave_velocity = np.zeros((ctes.n_components, ctes.n_internal_faces))
        if any(ponteiro):
            wave_velocity[:,ponteiro],m = RS.medium_wave_velocity(M, fprop, Nk_face, P_face, \
            Ft_internal, ponteiro)

            #wave_velocity[:,ponteiro] = 1e-10

        #Fk_face = RS.get_Fk_face(fprop, M, Nk_face, P_face, fprop.Vp[ctes.v0].flatten(), Ft_internal)
        #wave_velocity_RH = (Fk_face[...,1] - Fk_face[...,0])/(Nk_face[...,1] - Nk_face[...,0])
        #e = 1e-5
        #wave_velocity[abs(Nk_face[...,1] - Nk_face[...,0])>e] = wave_velocity_RH[abs(Nk_face[...,1] - Nk_face[...,0])>e]
        #import pdb; pdb.set_trace()
        #wave_velocity = Fk_vols_total/self.Nk #np.max(abs(wave_velocity),axis=0)


        #burger
        '''a[:,0:-1] = wave_velocity
        Fk_face_all = 0.5*(Fk_faces_all.sum(axis=-1)-abs(a)*(Nk_faces_all[...,1] - Nk_faces_all[...,0]))
        Fk_face_all = Fk_faces_all[...,0]
        Fk_face_all[a<=0] = Fk_faces_all[a<=0,1]
        Fk_internal_faces = Fk_face_all[:,:-1]
        self.vols_vec = -np.ones((ctes.n_volumes,2),dtype=int)
        lines = np.arange(ctes.n_internal_faces)
        self.vols_vec[ctes.v0[:,0],1] = lines
        self.vols_vec[ctes.v0[:,1],0] = lines

        Fk_vols_RS_neig = Fk_internal_faces[:,self.vols_vec]
        Fk_vols_RS_neig[:,self.vols_vec<0] = Fk_face_all[:,-1] #FOR THE BURGERS
        Fk_vols_total = Fk_vols_RS_neig[...,0] - Fk_vols_RS_neig[...,1]'''

        return Fk_vols_total, wave_velocity


    def get_faces_properties_upwind(self, fprop, G):
        ''' Using one-point upwind approximation '''
        fprop.mobilities_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.mobilities)
        fprop.Csi_j_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.Csi_j)
        fprop.xkj_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.xkj)

    def upwind(self, P, Pcap, G, prop):
        Pot_hid = P + Pcap - G[0,:,:]
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]

        #Pot_hidj_up = -Pot_hidj
        prop_face = np.zeros([prop.shape[0], prop.shape[1], ctes.n_internal_faces])
        prop_vols = prop[:,:,ctes.v0[:,0]]
        prop_vols_up = prop[:,:,ctes.v0[:,1]]
        prop_face[:,Pot_hidj_up <= Pot_hidj] = prop_vols[:,Pot_hidj_up <= Pot_hidj]
        prop_face[:,Pot_hidj_up > Pot_hidj] = prop_vols_up[:,Pot_hidj_up > Pot_hidj]
        return prop_face
