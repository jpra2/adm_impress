from .riemann_solvers import RiemannSolvers
from packs.utils import constants as ctes
from packs.directories import data_loaded
import scipy.sparse as sp
import numpy as np
from .flux_volumes import Flux


class FirstOrder:
    def __init__(self):
        pass

    def run(self,  M, fprop, ft_internal, P_old, Nk_old, G):
        if ctes.RS['LLF']:
            Fk_vols_total, wave_velocity = self.LLF(M, fprop, \
            ft_internal, P_old, Nk_old)
        elif ctes.RS['DW']:
            Fk_vols_total, wave_velocity = self.DW(M, fprop, ft_internal, P_old)
        elif ctes.RS['MDW']:
            Fk_vols_total, wave_velocity = self.MDW(M, fprop, ft_internal, P_old)
        elif ctes.RS['ROE']:
            Fk_vols_total, wave_velocity = self.ROE(M, fprop, ft_internal, P_old)
        else:
            self.get_faces_properties_upwind(fprop, G)
            Fk_vols_total, wave_velocity = self.FOU(M, fprop, \
                ft_internal, fprop.mobilities_internal_faces)
        return Fk_vols_total, wave_velocity

    def FOU(self, M, fprop, Ft_internal, mobilities_internal_faces):
        UPW = Flux()
        Fk_vols_total = UPW.update_flux(M, fprop, Ft_internal,
                             fprop.rho_j_internal_faces, mobilities_internal_faces)
        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)
        Nk_face = fprop.Nk[:,ctes.v0]#.sum(axis=-1)/2
        P_face = fprop.P[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)

        ponteiro = np.ones(ctes.n_internal_faces,dtype=bool)
        dNkmax_small = np.max(abs(Nk_face[:,:,0]-Nk_face[:,:,1]),axis=0)<1e-20
        ponteiro[dNkmax_small] = False
        wave_velocity = np.zeros((ctes.n_components, ctes.n_internal_faces))
        if any(ponteiro):
            wave_velocity[:,ponteiro],m = RS.medium_wave_velocity(M, fprop, Nk_face, P_face, \
            Ft_internal, ponteiro)
        #wave_velocity[:,ponteiro] = 1e-10
        #Fk_face = RS.get_Fk_face(fprop, M, Nk_face, fprop.P[ctes.v0], Ft_internal)
        #wave_velocity_RH = (Fk_face[...,1] - Fk_face[...,0])/(Nk_face[...,1] - Nk_face[...,0])
        #e = 1e-5
        #wave_velocity[abs(Nk_face[...,1] - Nk_face[...,0])>e] = wave_velocity_RH[abs(Nk_face[...,1] - Nk_face[...,0])>e]
        #import pdb; pdb.set_trace()
        #wave_velocity = Fk_vols_total/fprop.Nk #np.max(abs(wave_velocity),axis=0)
        return Fk_vols_total, wave_velocity

    def LLF(self, M, fprop, Ft_internal, P_old, Nk_old):

        Nk_face = Nk_old[:,ctes.v0]
        P_face = P_old[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)
        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)

        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, P_face, Ft_internal)
        Fk_internal_faces, wave_velocity = RS.LLF(M, fprop, Nk_face, P_face,
            Ft_internal, Fk_face, np.ones_like(Ft_internal[0],dtype=bool))
        '''wave_velocity_LR = wave_velocity_5[...,:2]
        upw = (wave_velocity_LR[...,0] * wave_velocity_LR[...,1])<=0
        upw[(wave_velocity_LR[...,0]==0) * (wave_velocity_LR[...,1]==0)] = False
        upw_cond = np.sum(upw, axis=0, dtype=bool)'''
        # contour faces upwind clássico (organize this later)
        Fk_internal_faces[:,[0,1]] = Fk_face[:,[0,1],0]

        Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)
        return Fk_vols_total, wave_velocity

    def MDW(self, M, fprop, Ft_internal, P_old):

        Nk_face = fprop.Nk[:,ctes.v0]
        P_face = P_old[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)

        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)

        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, P_face, Ft_internal)
        ponteiro = np.zeros(ctes.n_internal_faces,dtype=bool)
        Fk_internal_faces, wave_velocity = RS.MDW(M, fprop, Nk_face, P_face,
            Ft_internal, Fk_face, ~ponteiro)

        # contour faces upwind clássico (organize this later)
        Fk_internal_faces[...,0] = Fk_face[:,0,0]
        #Fk_internal_faces[...,2] = Fk_face[:,2,0]
        #Fk_internal_faces[...,1] = Fk_face[:,1,0]

        Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)
        return Fk_vols_total, wave_velocity

    def DW(self, M, fprop, Ft_internal, P_old):

        Nk_face = fprop.Nk[:,ctes.v0]
        P_face = P_old[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)

        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)

        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, P_face, Ft_internal)
        ponteiro = np.zeros(ctes.n_internal_faces,dtype=bool)
        Fk_internal_faces, wave_velocity = RS.DW(M, fprop, Nk_face, P_face,
            Ft_internal, Fk_face)

        # contour faces upwind clássico (organize this later)
        Fk_internal_faces[...,0] = Fk_face[:,0,0]
        #Fk_internal_faces[...,1] = Fk_face[:,1,0]

        Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)
        return Fk_vols_total, wave_velocity

    def ROE(self, M, fprop, Ft_internal, P_old):
        Nk_face = fprop.Nk[:,ctes.v0]
        P_face = P_old[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)
        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)

        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, P_face, Ft_internal)
        Fk_internal_faces, wave_velocity = RS.ROE_MDW(M, fprop, Nk_face, P_face,
            Ft_internal, Fk_face)

        # contour faces upwind clássico (organize this later)
        Fk_internal_faces[...,0] = Fk_face[:,0,0]
        #Fk_internal_faces[...,1] = Fk_face[:,1,0]

        Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)
        return Fk_vols_total, wave_velocity

    def get_faces_properties_upwind(self, fprop, G):
        ''' Using one-point upwind approximation '''
        Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]

        fprop.mobilities_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.mobilities)
        fprop.Csi_j_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.Csi_j)
        fprop.xkj_internal_faces = self.upwind(fprop.P, fprop.Pcap, G, fprop.xkj)

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
