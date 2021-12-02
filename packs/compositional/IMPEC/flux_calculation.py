import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from .properties_calculation import PropertiesCalc
from .composition_solver import RK3, Euler
import scipy.sparse as sp
from packs.compositional import prep_FR as ctes_FR
from packs.compositional import prep_MUSCL as ctes_MUSCL
import math
import time

if data_loaded['compositional_data']['component_data']['constant_K']:
    from packs.compositional.Kflash import StabilityCheck
else:
    from packs.compositional.stability_check import StabilityCheck

'Todo esse código vai ainda ser ajeitado! As funções de Flux devem ser funçoes para \
calculo independente do fluxo. Funcoes de rotina. E algumas funçes do MUSCL devem também \
seguir essa logica (sao funçoes que nao pertencem ao metodo em si mas sao auxiliares a ele \
e acabam sendo necessarias a outros metodos tambem)'

class Flux:
    """ Class created for computing flux accondingly to the First Order Upwind \
    Method - actually, it's only Flux because of the choice of mobilities and \
    other properties to be taken as function of the potencial gradient. But here \
    flux is computed at the interfaces, and then, volumes, only that."""

    def update_flux(self, M, fprop, Ft_internal_faces, rho_j_internal_faces,
        mobilities_internal_faces):
        ''' Main function that calls others '''
        self.Nk = fprop.Nk
        Fj_internal_faces = self.update_Fj_internal_faces(Ft_internal_faces,
            rho_j_internal_faces, mobilities_internal_faces, fprop.Pcap[:,ctes.v0],
            ctes.z[ctes.v0], ctes.pretransmissibility_internal_faces)
        Fk_internal_faces = self.update_Fk_internal_faces(
            fprop.xkj_internal_faces, fprop.Csi_j_internal_faces, Fj_internal_faces)
        Fk_vols_total = self.update_flux_volumes(Fk_internal_faces)
        return Fk_vols_total

    def update_Fj_internal_faces(self, Ft_internal_faces, rho_j_internal_faces,
        mobilities_internal_faces, Pcap_face, z_face,
        pretransmissibility_internal_faces):
        ''' Function to calculate phase flux '''

        frj = mobilities_internal_faces[0,...] / \
            np.sum(mobilities_internal_faces[0,...], axis = 0)

        Fj_internal_faces = frj[np.newaxis,...] * (Ft_internal_faces +
            pretransmissibility_internal_faces * (np.sum(mobilities_internal_faces *
            (Pcap_face[:,:,1] - Pcap_face[:,:,0] - ctes.g * rho_j_internal_faces *
            (z_face[:,1] - z_face[:,0])), axis=1) - np.sum(mobilities_internal_faces,
            axis=1) * (Pcap_face[:,:,1] - Pcap_face[:,:,0] - ctes.g *
            rho_j_internal_faces * (z_face[:,1] - z_face[:,0]))))
        return Fj_internal_faces

    def update_Fk_internal_faces(self, xkj_internal_faces, Csi_j_internal_faces, Fj_internal_faces):
        ''' Function to compute component molar flux '''
        Fk_internal_faces = np.sum(xkj_internal_faces * Csi_j_internal_faces *
            Fj_internal_faces, axis = 1)
        return Fk_internal_faces

    def update_flux_volumes(self, Fk_internal_faces):
        ''' Function to compute component molar flux balance through the control \
        volume interfaces '''
        cx = np.arange(ctes.n_components)
        lines = np.array([np.repeat(cx,len(ctes.v0[:,0])), np.repeat(cx,len(ctes.v0[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(ctes.v0[:,0],ctes.n_components), np.tile(ctes.v0[:,1], ctes.n_components)]).flatten()
        data = np.array([-Fk_internal_faces, Fk_internal_faces]).flatten()
        Fk_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_components, ctes.n_volumes)).toarray()

        '''
        inds = np.array([0,-1])
        a = (self.Nk[:,inds] - abs(self.Nk[:,inds]))/2
        a = a[:,self.v0]
        if a>=0:
            Fk_contour_face = self.Nk[0,-1] **2/2*len(self.Nk[-1])
        else:
            Fk_contour_face = self.Nk[0,0]**2/2*len(self.Nk[0])
        Fk_contour_face = self.Nk[0,-1]**2/2*len(self.Nk[0])
        Fk_vols_total[:,0] = -(Fk_internal_faces[:,0] - Fk_contour_face)
        Fk_vols_total[:,-1] = -(Fk_contour_face - Fk_internal_faces[:,1])'''

        return Fk_vols_total

class RiemannSolvers:
    def __init__(self, v0, pretransmissibility):
        self.v0 = v0
        self.pretransmissibility = pretransmissibility

    def ROE_entropy_correction(self, wave_velocity):
        wave_velocity_LR = wave_velocity[:,:,:2]
        C1 = ((wave_velocity_LR[:,:,0]<=0) * (wave_velocity_LR[:,:,1]>=0))
        C11 = ~C1
        Cc1 = ~np.sum(C11,axis=0,dtype=bool)
        wave_velocity_M = wave_velocity[:,:,-1]
        Cc2 = self.umbilic_points(wave_velocity_M)
        return Cc2 + Cc1

    def umbilic_points(self, wave_velocity_M):
        difs = np.empty((ctes.n_components, ctes.n_components,ctes.n_internal_faces))
        ind = np.arange(ctes.n_components).astype(int)
        for k in range(ctes.n_components):
            difs[k] = abs(wave_velocity_M[k,:] - wave_velocity_M[ind,:])
            difs[k,k] = 1e5
        C2 = np.min(abs(difs),axis = 1)
        E = 1e-8 #1e-9
        C22 = np.min(C2,axis=0)
        Cc2 = C22 < E*np.max(abs(wave_velocity_M),axis=0)
        #Cc2 = ~np.sum(~Cc2, axis=0,dtype=bool)
        return Cc2

    def LLF(self, M, fprop, Nk_face, P_face, ftotal, Fk_face, ponteiro_LLF):
        ponteiro = np.ones_like(ftotal[0],dtype=bool)
        #dNk_small = abs(Nk_face[:,:,0]-Nk_face[:,:,1])<1e-25
        #dNkmax_small = np.max(abs(Nk_face[:,:,0]-Nk_face[:,:,1]),axis=0)<1e-20
        #ponteiro[dNkmax_small] = False
        alpha_5 = np.empty((ctes.n_components, ctes.n_internal_faces, 5))
        if any(ponteiro):
            alpha_5[:,ponteiro,:], eigvec_m = self.wave_velocity_LLF(M, fprop, Nk_face,
                                P_face, ftotal, ponteiro)
        #alpha_5[:,~ponteiro] = 0
        #alpha_5[dNk_small] = 0
        alpha_LLF = np.max(abs(alpha_5),axis=0) #* ctes.ds_faces[:,np.newaxis]

        Fk_internal_faces = np.zeros_like(Fk_face[...,0])

        Fk_internal_faces[...,ponteiro_LLF] = self.update_flux_LLF(Fk_face[:,ponteiro_LLF],
            Nk_face[:,ponteiro_LLF], alpha_LLF[ponteiro_LLF,:])
        #alpha_corr = self.Harten_entropy_corr(alpha_5)
        #Fk_internal_faces[...,~ponteiro_LLF] = self.update_flux_ROE(Fk_face[:,~ponteiro_LLF],
        #    Nk_face[:,~ponteiro_LLF], alpha_5[:,~ponteiro_LLF,-1], eigvec_m[...,~ponteiro_LLF])
        return Fk_internal_faces, alpha_LLF

    def MDW(self, M, fprop, Nk_face, P_face, ftotal, Fk_face, ponteiro_MDW, alpha_m = []):
        ponteiro = np.ones_like(ponteiro_MDW[ponteiro_MDW], dtype=bool)
        #ponteiro[~np.sum((Nk_face[:,ponteiro_MDW,0]!=Nk_face[:,ponteiro_MDW,1]),axis=0,dtype=bool)] = True
        #dNkmax_small = np.max(abs(Nk_face[:,:,0]-Nk_face[:,:,1]),axis=0)<1e-30
        #ponteiro[dNkmax_small] = False

        alpha_MDW_ = np.zeros((len(ponteiro_MDW[ponteiro_MDW]),4))
        if any(ponteiro):
            alpha_MDW_[ponteiro] = self.wave_velocity_MDW( M, fprop, Nk_face,
                    P_face, ftotal, Fk_face, np.copy(ponteiro))
            alpha_MDW = np.max(abs(alpha_MDW_),axis=-1)

        alphaLR = self.LR_wave_velocity(M, fprop, Nk_face, P_face, \
            ftotal, np.copy(ponteiro))
        alphaLRmax = np.max(abs(alphaLR),axis=0)
        alpha_max = np.max(alphaLRmax,axis=-1)
        alpha_min = np.min(alphaLRmax,axis=-1)
        alpha_MDWl = alpha_MDW[ponteiro]

        alpha_MDWl[abs(alpha_MDWl)<alpha_min] = alpha_min[abs(alpha_MDWl)<alpha_min]
        alpha_MDWl[abs(alpha_MDWl)>alpha_max] = alpha_max[abs(alpha_MDWl)>alpha_max]
        alpha_MDW[ponteiro] = alpha_MDWl
        alpha_MDW[ponteiro] = self.Harten_entropy_corr(alpha_MDWl, alphaLR)
        Fk_internal_faces = self.update_flux_MDW(Fk_face[:,ponteiro_MDW], \
            Nk_face[:,ponteiro_MDW], alpha_MDW)
        return Fk_internal_faces, alpha_MDW

    def DW(self, M, fprop, Nk_face, P_face, ftotal, Fk_face, ponteiro_DW):
        ponteiro = np.ones_like(ponteiro_DW[ponteiro_DW], dtype=bool)
        dNkmax_small = np.max(abs(Nk_face[:,:,0]-Nk_face[:,:,1]),axis=0)<1e-25
        ponteiro[dNkmax_small] = False

        alpha_DW = np.zeros((len(ponteiro_DW[ponteiro_DW])))
        if any(ponteiro):
            alpha_DW[ponteiro] = self.wave_velocity_DW(M, fprop, Nk_face,
                    P_face, ftotal, Fk_face, np.copy(ponteiro))
            #alpha_DW = np.max(abs(alpha_DW_),axis=-1)
        alphaLR = self.LR_wave_velocity(M, fprop, Nk_face, P_face, \
            ftotal, np.copy(ponteiro))
        alphaLRmax = np.max(abs(alphaLR),axis=0)
        alpha_max = np.max(alphaLRmax,axis=-1)
        alpha_min = np.min(alphaLRmax,axis=-1)
        alpha_DWl = alpha_DW[ponteiro]

        alpha_DWl[abs(alpha_DWl)<alpha_min] = alpha_min[abs(alpha_DWl)<alpha_min]
        alpha_DWl[abs(alpha_DWl)>alpha_max] = alpha_max[abs(alpha_DWl)>alpha_max]
        alpha_DW[ponteiro] = alpha_DWl
        alpha_DW[ponteiro] = self.Harten_entropy_corr(alpha_DWl, alphaLR)
        Fk_internal_faces = self.update_flux_MDW(Fk_face[:,ponteiro_DW], \
            Nk_face[:,ponteiro_DW], alpha_DW)
        return Fk_internal_faces, alpha_DW

    def ROE_MDW(self, M, fprop, Nk_face, P_face, ftotal, Fk_face):
        ponteiro = np.ones_like(ftotal[0],dtype=bool)
        Vpm = fprop.Vp[ctes.v0].sum(axis=-1)/2
        alpha_LR, alpha_m_or, eigvec_m = self.LRM_wave_velocity(M, fprop, Nk_face, P_face, \
            ftotal, ponteiro)
        alphas = alpha_m_or[...,np.newaxis]

        alpha_m = self.Harten_entropy_corr(alpha_m_or, alpha_LR)

        ponteiro_MDW = self.umbilic_points(alpha_m)

        Fk_internal_faces = np.empty_like(Fk_face[...,0])

        Fk_internal_faces[...,ponteiro_MDW], alpha_MDW = self.MDW(M, fprop, Nk_face, \
            P_face, ftotal, Fk_face, ponteiro_MDW, alphas[:,ponteiro_MDW])

        eigvec_m_ROE = eigvec_m[...,~ponteiro_MDW]
        real_eigvecs = np.isreal(eigvec_m_ROE)

        if any(~real_eigvecs.flatten()): import pdb; pdb.set_trace()

        eigvec_m_ROE = np.real(eigvec_m_ROE)

        Fk_internal_faces[...,~ponteiro_MDW] = self.update_flux_ROE(Fk_face[:,~ponteiro_MDW],
            Nk_face[:,~ponteiro_MDW], alpha_m[:,~ponteiro_MDW], eigvec_m_ROE)

        alpha = np.copy(alpha_m)
        alpha[:,ponteiro_MDW] = alpha_MDW
        return Fk_internal_faces, alpha

    def ROE_LLF(self, M, fprop, Nk_face, P_face, ftotal, Fk_face):
        ponteiro = np.ones_like(ftotal[0],dtype=bool)
        Vpm = fprop.Vp[ctes.v0].sum(axis=-1)/2
        alpha_5, eigvec_m = self.wave_velocity_LLF(M, fprop, Nk_face, P_face, \
            ftotal, ponteiro)
        alpha_m = alpha_5[...,-1]
        alpha_LLF = np.max(abs(alpha_5),axis=0)
        alpha_m = self.Harten_entropy_corr(alpha_m, alpha_5[...,[0,1]])
        ponteiro_LLF = self.umbilic_points(alpha_m)
        Fk_internal_faces = np.empty_like(Fk_face[...,0])

        Fk_internal_faces[...,ponteiro_LLF] = self.update_flux_LLF(Fk_face[:,ponteiro_LLF],
            Nk_face[:,ponteiro_LLF], alpha_LLF[ponteiro_LLF,:])

        eigvec_m_ROE = eigvec_m[...,~ponteiro_LLF]
        real_eigvecs = np.isreal(eigvec_m_ROE)

        if any(~real_eigvecs.flatten()): import pdb; pdb.set_trace()

        eigvec_m_ROE = np.real(eigvec_m_ROE)

        Fk_internal_faces[...,~ponteiro_LLF] = self.update_flux_ROE(Fk_face[:,~ponteiro_LLF],
            Nk_face[:,~ponteiro_LLF], alpha_m[:,~ponteiro_LLF], eigvec_m_ROE)

        alpha = np.copy(np.max(abs(alpha_m),axis=0))
        alpha[ponteiro_LLF] = np.max(abs(alpha_LLF[ponteiro_LLF,:]),axis=-1)
        return Fk_internal_faces, alpha

    def reshape_constant_property(self, y, ponteiro, v):

        y_aux0 = y[self.v0][:,ponteiro,0]
        y_aux1 = y[self.v0][:,ponteiro,1]

        y_reshaped1 = np.tile(y_aux0,(int(v/2) * \
            (np.sign(v - ctes.n_components)**2) + int(v)*(1 - np.sign(v -
            ctes.n_components)**2)))
        y_reshaped2 = np.tile(y_aux1,(int(v/2) * \
            (np.sign(v - ctes.n_components)**2) + int(v/(ctes.n_components+1))*(1 -
            np.sign(v-ctes.n_components))))

        y_reshaped = np.concatenate((y_reshaped1, y_reshaped2), axis=-1)

        if v==3:
            y_reshaped3 = (y_aux0 + y_aux1)*0.5
            y_reshaped = np.concatenate((y_reshaped, y_reshaped3))

        if ctes.FR:
            if self.v0.shape[1] >= 3:
                y_reshaped = y_aux0
                for i in range(1,self.v0.shape[1]):
                    try:
                        y_reshaped_i = y[self.v0][ponteiro,i]
                    except:
                        y_reshaped_i = y[:,self.v0][:,ponteiro,i]
                    y_reshaped = np.concatenate((y_reshaped, y_reshaped_i), axis=-1)
        return y_reshaped

    def get_extrapolated_properties(self, fprop, M, Nk_face, z_face, P_face, Vp_face, v, ponteiro):
        xkj_face = np.empty((ctes.n_components, ctes.n_phases, len(P_face)))
        Csi_j_face = np.empty((1, ctes.n_phases, len(P_face)))
        rho_j_face = np.empty_like(Csi_j_face)

        ''' Flash calculations and properties calculations at each side of the \
        interface '''
        if ctes.load_k:
            if ctes.compressible_k:
                L_face, V_face, xkj_face[0:ctes.Nc,0,...], xkj_face[0:ctes.Nc,1,...], \
                Csi_j_face[0,0,...], Csi_j_face[0,1,...], rho_j_face[0,0,...], \
                rho_j_face[0,1,...] = StabilityCheck(P_face, fprop.T).run_init(P_face, z_face)

            else:
                L_face = np.ones(len(P_face)); V_face = np.zeros(len(P_face))
                xkj_face[0:ctes.Nc,0:2,:] = 1

                rho_j_face[0,0,:] = np.tile(fprop.rho_j[0,0,self.v0[ponteiro,0]], v)
                rho_j_face[0,1,:] = np.tile(fprop.rho_j[0,1,self.v0[ponteiro,0]], v) #constante, independe da pressao
                #self.reshape_constant_property(fprop.rho_j[0,0:2,:], ponteiro, v)
                Csi_j_face[0,0,:] = np.tile(fprop.Csi_j[0,0,self.v0[ponteiro,0]], v)
                Csi_j_face[0,1,:] = np.tile(fprop.Csi_j[0,1,self.v0[ponteiro,0]], v)
                #self.reshape_constant_property(fprop.Csi_j[0,0:2,:], ponteiro, v)
        else: L_face = []; V_face = []

        if ctes.load_w:
            xkj_face[-1,-1,...] = 1
            xkj_face[-1,0:-1,...] = 0
            xkj_face[0:ctes.Nc,-1,...] = 0

            if data_loaded['compositional_data']['water_data']['mobility']:
                # Csi_W0 é constante independente de qualquer coisa(por teoria)
                Csi_W0_face = np.tile(fprop.Csi_W0[self.v0[ponteiro,0]], v)
                #self.reshape_constant_property(fprop.Csi_W0, ponteiro, v)

                Sw_face, Csi_j_face[0,-1,...], rho_j_face[0,-1,...] = \
                PropertiesCalc().update_water_saturation(fprop, Nk_face[-1,...],
                P_face, Vp_face, Csi_W0_face)

            else:
                # se não movel, prop constante (no varia com P)
                Sw_face = np.tile(fprop.Sw[self.v0[ponteiro,0]], v)
                #self.reshape_constant_property(fprop.Sw, ponteiro, v)
                rho_j_face[0,-1] = np.tile(fprop.rho_j[0,-1,self.v0[ponteiro,0]], v)
                #self.reshape_constant_property(fprop.rho_j[0,-1], ponteiro, v)
                Csi_j_face[0,-1] = np.tile(fprop.Csi_j[0,-1,self.v0[ponteiro,0]], v)
                #self.reshape_constant_property(fprop.Csi_j[0,-1], ponteiro, v)
        else:
            Sw_face = np.tile(fprop.Sw[self.v0[ponteiro,0]],v)

        So_face, Sg_face =  PropertiesCalc().update_saturations(Sw_face,
            Csi_j_face, L_face, V_face)

        mobilities_face = PropertiesCalc().update_mobilities(fprop, So_face,
            Sg_face, Sw_face, Csi_j_face, xkj_face)
        return mobilities_face, rho_j_face, Csi_j_face, xkj_face

    def Fk_from_Nk(self, fprop, M, Nk, P_face, Vp_face, ftotal, ponteiro):
        ''' Function to compute component flux based on a given composition (Nk) '''
        v = int(len(Nk[0,:])/len(ponteiro[ponteiro]))
        z = Nk[0:ctes.Nc] / np.sum(Nk[0:ctes.Nc], axis = 0)
        #P_face = np.tile(P_face,v)
        #P_face_resh = np.concatenate(np.hsplit(P_face,v),axis=0)[:,0]

        mobilities, rho_j, Csi_j, xkj = self.get_extrapolated_properties(fprop, M, Nk, z,
            P_face, Vp_face, v, ponteiro)

        f = Flux()

        Pcap_reshaped = np.concatenate(np.dsplit(np.tile(fprop.Pcap[:,self.v0][:,ponteiro],v),v),axis=1)
        z_reshaped = np.concatenate(np.hsplit(np.tile(ctes.z[self.v0][ponteiro],v),v),axis=0)
        Fj = f.update_Fj_internal_faces(ftotal, rho_j, mobilities, Pcap_reshaped,
            z_reshaped, np.tile(self.pretransmissibility[ponteiro],v))
        Fk = f.update_Fk_internal_faces(xkj, Csi_j, Fj)
        return Fk

    def Fk_and_rho_from_Nk(self, fprop, M, Nk, P_face, Vp_face, ftotal, ponteiro):
        ''' Function to compute component flux based on a given composition (Nk) '''
        v = int(len(Nk[0,:])/len(ponteiro[ponteiro]))
        z = Nk[0:ctes.Nc] / np.sum(Nk[0:ctes.Nc], axis = 0)
        #P_face = np.tile(P_face,v)
        #P_face_resh = np.concatenate(np.hsplit(P_face,v),axis=0)[:,0]

        mobilities, rho_j, Csi_j, xkj = self.get_extrapolated_properties(fprop, M, Nk, z,
            P_face, Vp_face, v, ponteiro)

        f = Flux()

        Pcap_reshaped = np.concatenate(np.dsplit(np.tile(fprop.Pcap[:,self.v0][:,ponteiro],v),v),axis=1)
        z_reshaped = np.concatenate(np.hsplit(np.tile(ctes.z[self.v0][ponteiro],v),v),axis=0)
        Fj = f.update_Fj_internal_faces(ftotal, rho_j, mobilities, Pcap_reshaped,
            z_reshaped, np.tile(self.pretransmissibility[ponteiro],v))
        Fk = f.update_Fk_internal_faces(xkj, Csi_j, Fj)
        return Fk, rho_j

    def get_Fk_face(self, fprop, M, Nk_face, P_face, ftotal):
        ft_Nks = np.tile(ftotal,2)
        Nks = np.concatenate((Nk_face[:,:,0], Nk_face[:,:,1]), axis=1)
        Vp_face = np.concatenate((fprop.Vp[ctes.v0[:,0]], fprop.Vp[ctes.v0[:,1]]), axis=0)
        P_faces = np.concatenate((P_face[:,0],P_face[:,1]))
        Fk_faces = self.Fk_from_Nk(fprop, M, Nks, P_faces, Vp_face, ft_Nks, np.ones_like(ftotal[0], dtype=bool))
        Fk_faceL, Fk_faceR = np.hsplit(Fk_faces, 2)
        Fk_face = np.concatenate((Fk_faceL[:,:,np.newaxis], Fk_faceR[:,:,np.newaxis]), axis=-1)
        return Fk_face

    def medium_wave_velocity(self, M, fprop, Nk_face, P_face, ftotal, ponteiro):
        Nkm = (Nk_face[:,ponteiro,1] + Nk_face[:,ponteiro,0])/2
        delta = 1e-3 * abs(Nkm)
        delta[Nkm==0] = 1e-13

        #ponteiro = np.ones_like(ftotal[0, ponteiro], dtype=bool)
        Nk_aux_matrix = np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        matrix_deltas = np.identity(ctes.n_components)[:,:,np.newaxis, np.newaxis] * Nk_aux_matrix
        delta_05 = delta[...,np.newaxis] * 0.5 * matrix_deltas

        Nkm_aux = Nkm[np.newaxis,:,:] * Nk_aux_matrix[...,0]
        Nkm_plus = np.copy(Nkm_aux)
        Nkm_minus = np.copy(Nkm_aux)
        Nkm_plus  += delta_05[...,0]
        Nkm_minus -= delta_05[...,0]
        Nkm_plus[Nkm_minus<=0] = Nkm_plus[Nkm_minus <= 0] + delta_05[Nkm_minus<=0,0]
        Nkm_minus[Nkm_minus<=0] = Nkm_aux[Nkm_minus <= 0]
        Nkms = np.concatenate((Nkm_plus[...,np.newaxis], Nkm_minus[...,np.newaxis]),axis=-1)
        Nkms = np.concatenate(np.split(Nkms, ctes.n_components),axis=2)[0,...]
        Nkms = np.concatenate(np.dsplit(Nkms, 2),axis=1)[...,0]
        Nkm_plus, Nkm_minus = np.hsplit(Nkms,2)

        ft_Nks = np.tile(ftotal[:,ponteiro],ctes.n_components*2)
        Vp = fprop.Vp[ctes.v0]
        Vpm = np.sum(Vp[ponteiro],axis=-1)*0.5
        Vps = np.tile(Vpm, ctes.n_components*2)
        P_faces = np.tile(P_face[ponteiro, 0], ctes.n_components*2)

        Fks = self.Fk_from_Nk(fprop, M, Nkms, P_faces, Vps, ft_Nks, ponteiro)

        Fk_plus, Fk_minus = np.hsplit(Fks,2)
        dFkdNk = (Fk_plus - Fk_minus)/ (Nkm_plus - Nkm_minus).sum(axis=0)
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,np.newaxis],ctes.n_components),axis = 2)
        dFkdNk[np.isnan(dFkdNk)] = 0
        dFkdNk = dFkdNk.transpose(1,0,2)

        eigval1, eigvec = np.linalg.eig(dFkdNk)
        dFkdNk_eigvalue = eigval1.T
        dFdNk_eigvec = eigvec.transpose(1,2,0)
        return dFkdNk_eigvalue, dFdNk_eigvec

    def LR_wave_velocity(self, M, fprop, Nk_face, P_face, ftotal, ponteiro):

        delta = 1e-3 * abs(Nk_face[:,ponteiro])
        delta[Nk_face[:,ponteiro]==0] = 1e-13

        Nk_aux_matrix = np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        matrix_deltas = np.identity(ctes.n_components)[:,:,np.newaxis, np.newaxis] * Nk_aux_matrix
        delta_05 = delta * 0.5 * matrix_deltas

        Nk_face_aux = Nk_face[np.newaxis,:,ponteiro,:] * Nk_aux_matrix
        Nk_face_plus = np.copy(Nk_face_aux)
        Nk_face_minus = np.copy(Nk_face_aux)

        Nk_face_plus += delta_05
        Nk_face_minus -= delta_05
        Nk_face_plus[Nk_face_minus <= 0] = Nk_face_plus[Nk_face_minus <= 0] + delta_05[Nk_face_minus<=0]
        Nk_face_minus[Nk_face_minus <= 0] = Nk_face_aux[Nk_face_minus <= 0]
        Nks = np.concatenate((Nk_face_plus,Nk_face_minus),axis=-1)
        Nks = np.concatenate(np.split(Nks, ctes.n_components),axis=2)[0,...]
        Nks = np.concatenate(np.dsplit(Nks, 4),axis=1)[...,0]
        Nk_face_plus, Nk_face_minus = np.hsplit(Nks,2)


        ft_Nks = np.tile(ftotal[:,ponteiro],ctes.n_components*4)
        Vp_faceL = np.tile(fprop.Vp[self.v0[ponteiro,0]], ctes.n_components)
        Vp_faceR = np.tile(fprop.Vp[self.v0[ponteiro,1]], ctes.n_components)
        Vp_plus = np.concatenate((Vp_faceL, Vp_faceR), axis=0)
        Vps = np.tile(Vp_plus, 2)
        P_face = np.tile(P_face[ponteiro,0],ctes.n_components*4)

        Fks = self.Fk_from_Nk(fprop, M, Nks, P_face, Vps, ft_Nks, ponteiro)
        Fk_faces_plus, Fk_faces_minus = np.hsplit(Fks,2)

        dFkdNk = (Fk_faces_plus -Fk_faces_minus)/(Nk_face_plus - Nk_face_minus).sum(axis=0)
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,np.newaxis],2),axis=2)
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,:,np.newaxis],ctes.n_components),axis=3)
        dFkdNk[np.isnan(dFkdNk)] = 0
        dFkdNk = dFkdNk.transpose(2,1,0,3)

        eigval1, v = np.linalg.eig(dFkdNk)
        dFkdNk_eigvalue = eigval1.T
        return dFkdNk_eigvalue

    def LRM_wave_velocity(self, M, fprop, Nk_face, P_face, ftotal, ponteiro):

        delta = 1e-3 * abs(Nk_face)
        delta[Nk_face==0] = 1e-13

        Nkm = (Nk_face[:,ponteiro,1] + Nk_face[:,ponteiro,0])/2
        deltam = 1e-3 * abs(Nkm)
        deltam[Nkm==0] = 1e-13

        Nk_aux_matrix = np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        matrix_deltas = np.identity(ctes.n_components)[:,:,np.newaxis, np.newaxis] * Nk_aux_matrix
        delta_05 = delta * 0.5 * matrix_deltas
        delta_05m = deltam[...,np.newaxis] * 0.5 * matrix_deltas

        Nk_face_aux = Nk_face[np.newaxis,:,ponteiro,:] * Nk_aux_matrix
        Nk_face_plus = np.copy(Nk_face_aux)
        Nk_face_minus = np.copy(Nk_face_aux)

        Nk_face_plus += delta_05
        Nk_face_minus -= delta_05
        Nk_face_plus[Nk_face_minus <= 0] = Nk_face_plus[Nk_face_minus <= 0] + delta_05[Nk_face_minus<=0]
        Nk_face_minus[Nk_face_minus <= 0] = Nk_face_aux[Nk_face_minus <= 0]
        Nks = np.concatenate((Nk_face_plus,Nk_face_minus),axis=-1)
        Nks = np.concatenate(np.split(Nks, ctes.n_components),axis=2)[0,...]
        Nks = np.concatenate(np.dsplit(Nks, 4),axis=1)[...,0]
        Nk_face_plus, Nk_face_minus = np.hsplit(Nks,2)

        Nkm_aux = Nkm[np.newaxis,:,:] * Nk_aux_matrix[...,0]
        Nkm_plus = np.copy(Nkm_aux)
        Nkm_minus = np.copy(Nkm_aux)
        Nkm_plus  += delta_05m[...,0]
        Nkm_minus -= delta_05m[...,0]
        Nkm_plus[Nkm_minus<=0] = Nkm_plus[Nkm_minus <= 0] + delta_05[Nkm_minus<=0,0]
        Nkm_minus[Nkm_minus<=0] = Nkm_aux[Nkm_minus <= 0]
        Nkms = np.concatenate((Nkm_plus[...,np.newaxis], Nkm_minus[...,np.newaxis]),axis=-1)
        Nkms = np.concatenate(np.split(Nkms, ctes.n_components),axis=2)[0,...]
        Nkms = np.concatenate(np.dsplit(Nkms, 2),axis=1)[...,0]
        Nks = np.concatenate((Nks, Nkms),axis=-1)
        Nkm_plus, Nkm_minus = np.hsplit(Nkms,2)

        ft_Nks = np.tile(ftotal[:,ponteiro],ctes.n_components*6)
        Vp_faceL = np.tile(fprop.Vp[self.v0[ponteiro,0]], ctes.n_components)
        Vp_faceR = np.tile(fprop.Vp[self.v0[ponteiro,1]], ctes.n_components)
        Vp_plus = np.concatenate((Vp_faceL, Vp_faceR), axis=0)
        Vp_faces = np.tile(Vp_plus, 2)
        Vpm = fprop.Vp[self.v0[ponteiro]].sum(axis=-1)*0.5
        Vpms = np.tile(Vpm, ctes.n_components*2)
        Vps = np.concatenate((Vp_faces, Vpms), axis=0)
        P_face = np.tile(P_face[ponteiro,0],3*ctes.n_components*2)

        Fks = self.Fk_from_Nk(fprop, M, Nks, P_face, Vps, ft_Nks, ponteiro)
        Fk_faces_plus, Fk_faces_minus, Fkms = np.hsplit(Fks,3)

        dFkdNk = (Fk_faces_plus -Fk_faces_minus)/(Nk_face_plus - Nk_face_minus).sum(axis=0)
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,np.newaxis],2),axis=2)
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,:,np.newaxis],ctes.n_components),axis=3)
        dFkdNk[np.isnan(dFkdNk)] = 0
        dFkdNk = dFkdNk.transpose(2,1,0,3)

        Fkm_plus, Fkm_minus = np.hsplit(Fkms,2)
        dFkdNk_m = (Fkm_plus - Fkm_minus)/ (Nkm_plus - Nkm_minus).sum(axis=0)
        dFkdNk_m = np.concatenate(np.hsplit(dFkdNk_m[:,:,np.newaxis],ctes.n_components),axis = 2)
        dFkdNk_m[np.isnan(dFkdNk_m)] = 0
        dFkdNk_m = dFkdNk_m.transpose(1,0,2)

        eigval1, v = np.linalg.eig(dFkdNk)
        dFkdNk_eigvalue = eigval1.T

        eigvalm, eigvec = np.linalg.eig(dFkdNk_m)
        dFkdNkm_eigvalue = eigvalm.T
        dFdNkm_eigvec = eigvec.transpose(1,2,0)
        return dFkdNk_eigvalue, dFkdNkm_eigvalue, dFdNkm_eigvec

    def wave_velocity_LLF(self, M, fprop, Nk_face, P_face, ftotal, ponteiro):
        #delta = 1e-12 #1e-3 * np.min(abs(Nk_face[Nk_face>1e-12]))
        delta = 1e-3 * abs(Nk_face[:,ponteiro])
        delta[abs(Nk_face[:,ponteiro])<=0] = 1e-13

        Nkm = (Nk_face[:,ponteiro,1] + Nk_face[:,ponteiro,0])/2
        deltam = 1e-3 * abs(Nkm)
        deltam[abs(Nkm)<=0] = 1e-13

        Nkg = Nkm[:,:,np.newaxis] + (Nk_face[:,ponteiro] - Nkm[:,:,np.newaxis])/(3**(1/2))
        deltag = 1e-3 * abs(Nkg)
        deltag[abs(Nkg)<=0] = 1e-13

        Nk_aux_matrix = np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        matrix_deltas = np.identity(ctes.n_components)[:,:,np.newaxis, np.newaxis] * Nk_aux_matrix
        delta_05 = delta * 0.5 * matrix_deltas
        delta_05m = deltam[...,np.newaxis] * 0.5 * matrix_deltas
        delta_05g = deltag * 0.5 * matrix_deltas

        Nkg_aux = Nkg[np.newaxis,...] * Nk_aux_matrix
        Nkg_plus = np.copy(Nkg_aux)
        Nkg_minus = np.copy(Nkg_aux)
        Nkg_plus += delta_05g
        Nkg_minus -= delta_05g
        Nkg_plus[Nkg_minus<=0] = Nkg_plus[Nkg_minus<=0] + delta_05[Nkg_minus<=0]
        Nkg_minus[Nkg_minus<=0] = Nkg_aux[Nkg_minus <= 0]
        #Nkg_plus[Nkg_aux==0] = 0
        Nkgs = np.concatenate((Nkg_plus,Nkg_minus),axis=-1)
        #Nkgs = np.concatenate(np.split(Nkgs, ctes.n_components),axis=2)[0,...]
        #Nkgs = np.concatenate(np.dsplit(Nkgs, 4),axis=1)[...,0]

        Nk_face_aux = Nk_face[np.newaxis,:,ponteiro,:] * Nk_aux_matrix
        Nk_face_plus = np.copy(Nk_face_aux)
        Nk_face_minus = np.copy(Nk_face_aux)
        Nk_face_plus += delta_05
        Nk_face_minus -= delta_05
        Nk_face_plus[Nk_face_minus <= 0] = Nk_face_plus[Nk_face_minus <= 0] + delta_05[Nk_face_minus<=0]
        Nk_face_minus[Nk_face_minus <= 0] = Nk_face_aux[Nk_face_minus <= 0]
        #Nk_face_plus[Nk_face_aux==0] = 0
        Nk_faces = np.concatenate((Nk_face_plus,Nk_face_minus),axis=-1)
        Nks = np.concatenate((Nkgs, Nk_faces), axis=-1)
        Nks = np.concatenate(np.split(Nks, ctes.n_components),axis=2)[0,...]
        Nks = np.concatenate(np.dsplit(Nks, 8),axis=1)[...,0]

        Nkg_plus, Nkg_minus, Nk_face_plus, Nk_face_minus = np.hsplit(Nks,4)

        Nkm_aux = Nkm[np.newaxis,:,:] * Nk_aux_matrix[...,0]
        Nkm_plus = np.copy(Nkm_aux)
        Nkm_minus = np.copy(Nkm_aux)
        Nkm_plus  += delta_05m[...,0]
        Nkm_minus -= delta_05m[...,0]
        Nkm_plus[Nkm_minus<=0] = Nkm_plus[Nkm_minus <= 0] + delta_05[Nkm_minus<=0,0]
        Nkm_minus[Nkm_minus<=0] = Nkm_aux[Nkm_minus <= 0]
        #Nkm_plus[Nkm_aux==0] = 0
        Nkms = np.concatenate((Nkm_plus[...,np.newaxis], Nkm_minus[...,np.newaxis]),axis=-1)
        Nkms = np.concatenate(np.split(Nkms, ctes.n_components),axis=2)[0,...]
        Nkms = np.concatenate(np.dsplit(Nkms, 2),axis=1)[...,0]
        Nkm_plus, Nkm_minus = np.hsplit(Nkms,2)
        Nks = np.concatenate((Nks, Nkms),axis=-1)

        ft_Nks = np.tile(ftotal[:,ponteiro],ctes.n_components*10)
        Vp_faceL = np.tile(fprop.Vp[self.v0[ponteiro,0]], ctes.n_components)
        Vp_faceR = np.tile(fprop.Vp[self.v0[ponteiro,1]], ctes.n_components)
        Vp_plus = np.concatenate((Vp_faceL, Vp_faceR), axis=0)
        Vp_faces = np.tile(Vp_plus, 4)
        Vpm = fprop.Vp[self.v0[ponteiro]].sum(axis=-1)*0.5
        Vpms = np.tile(Vpm, ctes.n_components*2)
        Vps = np.concatenate((Vp_faces, Vpms), axis=0)
        P_face = np.tile(P_face[ponteiro,0],5*ctes.n_components*2)

        Fks = self.Fk_from_Nk(fprop, M, Nks, P_face, Vps, ft_Nks, ponteiro)

        Fk_Nkg_plus, Fk_Nkg_minus, Fk_faces_plus, Fk_faces_minus, Fkms = np.hsplit(Fks,5)

        Fkm_plus, Fkm_minus = np.hsplit(Fkms,2)
        dFkdNk_m = (Fkm_plus - Fkm_minus)/ (Nkm_plus - Nkm_minus).sum(axis=0)
        dFkdNk_m = np.concatenate(np.hsplit(dFkdNk_m[:,:,np.newaxis],ctes.n_components),axis = 2)
        dFkdNk_m[np.isnan(dFkdNk_m)] = 0
        dFkdNk_m = dFkdNk_m.transpose(1,0,2)

        dFkdNk = ((Fk_faces_plus - Fk_faces_minus)/(Nk_face_plus - Nk_face_minus).sum(axis=0))
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,np.newaxis],2),axis=2)
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,:,np.newaxis],ctes.n_components),axis=3)
        dFkdNk[np.isnan(dFkdNk)] = 0
        dFkdNk = dFkdNk.transpose(2,1,0,3)

        dFkdNk_gauss = (Fk_Nkg_plus - Fk_Nkg_minus)/(Nkg_plus - Nkg_minus).sum(axis=0)
        dFkdNk_gauss = np.concatenate(np.hsplit(dFkdNk_gauss[:,:,np.newaxis],2),axis=2)
        dFkdNk_gauss = np.concatenate(np.hsplit(dFkdNk_gauss[:,:,:,np.newaxis],ctes.n_components),axis=3)
        dFkdNk_gauss[np.isnan(dFkdNk_gauss)] = 0
        if any(np.isnan(dFkdNk_gauss).flatten()):import pdb; pdb.set_trace()
        dFkdNk_gauss = dFkdNk_gauss.transpose(2,1,0,3)

        eigval1, v = np.linalg.eig(dFkdNk)
        dFkdNk_eigvalue = eigval1.T
        eigval2, v = np.linalg.eig(dFkdNk_gauss)
        dFkdNk_gauss_eigvalue = eigval2.transpose(2,1,0)
        eigval3, eigvec_m = np.linalg.eig(dFkdNk_m)
        dFkdNk_m_eigvalue = eigval3.T
        alpha = np.concatenate((dFkdNk_eigvalue, dFkdNk_gauss_eigvalue), axis=-1)
        alpha = np.concatenate((alpha, dFkdNk_m_eigvalue[:,:,np.newaxis]), axis=-1)
        #alpha[Nk_face.sum(axis=-1)==0] = 0
        #import pdb; pdb.set_trace()
        return alpha, eigvec_m.transpose(1,2,0)

    def wave_velocity_MDW(self, M, fprop, Nk_face, P_face, ftotal, Fk_face, ponteiro):
        Nkm = (Nk_face[:,ponteiro,1] + Nk_face[:,ponteiro,0])/2

        Nkg = Nkm[:,:,np.newaxis] + (Nk_face[:,ponteiro] - Nkm[:,:,np.newaxis])/(3**(1/2))

        Nkgs = np.concatenate((Nkg[:,:,0], Nkg[:,:,1]), axis=1)
        Nks = np.concatenate((Nkgs, Nkm), axis=1)

        ft_Nks = np.tile(ftotal[:,ponteiro],3)

        Vp_face = fprop.Vp[ctes.v0]
        Vp_faces = np.concatenate((Vp_face[ponteiro,0], Vp_face[ponteiro,1]), axis=0)
        Vpm = np.sum(Vp_face[ponteiro], axis=-1) * 0.5
        Vps = np.concatenate((Vp_faces, Vpm), axis=0)

        P_faces = np.concatenate((P_face[ponteiro,0], P_face[ponteiro,1]),axis=0)
        P_facem = np.sum(P_face[ponteiro],axis=-1) * 0.5
        P_faces = np.concatenate((P_faces, P_facem),axis=0)

        Fks = self.Fk_from_Nk(fprop, M, Nks, P_faces, Vps, ft_Nks, ponteiro)

        Fk_NkgL, Fk_NkgR, Fk_Nkm = np.hsplit(Fks,3)

        dNkg_Nkface = Nkg - Nk_face[:,ponteiro]
        dNkm_Nkg = Nkm[:,:,np.newaxis] - Nkg

        e = 0

        alpha_L_GL = np.sum((dNkg_Nkface[:,:,0]) * (Fk_NkgL - Fk_face[:,ponteiro,0]), axis = 0) / \
                    np.sum((dNkg_Nkface[:,:,0]) * (dNkg_Nkface[:,:,0]), axis = 0)
        alpha_L_GL[np.sum(((dNkg_Nkface[:,:,0]) * (dNkg_Nkface[:,:,0])), axis = 0)==e] = 0

        alpha_GL_M = np.sum((dNkm_Nkg[:,:,0]) * (Fk_Nkm - Fk_NkgL), axis=0) / \
                     np.sum((dNkm_Nkg[:,:,0]) * (dNkm_Nkg[:,:,0]), axis=0)
        alpha_GL_M[np.sum(((dNkm_Nkg[:,:,0]) * (dNkm_Nkg[:,:,0])), axis = 0)==e] = 0

        alpha_M_GR = np.sum((-dNkm_Nkg[:,:,1]) * (Fk_NkgR - Fk_Nkm), axis=0) / \
                    np.sum((-dNkm_Nkg[:,:,1]) * (-dNkm_Nkg[:,:,1]), axis=0)
        alpha_M_GR[np.sum(((-dNkm_Nkg[:,:,1]) * (-dNkm_Nkg[:,:,1])), axis=0)==e] = 0

        alpha_GR_R = np.sum((-dNkg_Nkface[:,:,1]) * (Fk_face[:,ponteiro,1] - Fk_NkgR), axis=0) / \
                    np.sum((-dNkg_Nkface[:,:,1]) * (-dNkg_Nkface[:,:,1]), axis=0)
        alpha_GR_R[np.sum(((-dNkg_Nkface[:,:,1]) * (-dNkg_Nkface[:,:,1])), axis=0)==e] = 0

        alpha_MDW = np.concatenate((alpha_L_GL[...,np.newaxis], alpha_GL_M[...,np.newaxis]), axis=-1)
        alpha_MDW = np.concatenate((alpha_MDW, alpha_M_GR[...,np.newaxis]), axis=-1)
        alpha_MDW = np.concatenate((alpha_MDW, alpha_GR_R[...,np.newaxis]), axis=-1)
        #alpha_MDW = np.max(abs(alpha_MDW),axis=-1)

        if any(np.isnan(alpha_MDW).flatten()): import pdb; pdb.set_trace()
        return alpha_MDW

    def wave_velocity_DW(self, M, fprop, Nk_face, P_face, ftotal, Fk_face, ponteiro):
        dNk = (Nk_face[:,ponteiro,1] - Nk_face[:,ponteiro,0])
        dFk = (Fk_face[:,ponteiro,1] - Fk_face[:,ponteiro,0])
        #dFk[dNk>0] = 0
        alpha_DW = np.sum((dNk) * (dFk), axis = 0) / \
                    np.sum((dNk) * (dNk), axis = 0)
        alpha_DW[np.sum((dNk) * (dNk), axis = 0)==0] = 0
        #import pdb; pdb.set_trace()
        return alpha_DW

    def update_flux_LLF(self, Fk_face_LLF_all, Nk_face_LLF, alpha_LLF):
        #alpha2 = np.concatenate((alpha_LLF, alpha_RH[:,:,np.newaxis]),axis=-1)
        '''alpha_RH = (Fk_face_LLF_all[:,:,1] - Fk_face_LLF_all[:,:,0]) / \
            (Nk_face_LLF[:,:,1] - Nk_face_LLF[:,:,0])
        alpha_RH[np.isnan(alpha_RH)] = 0
        alpha_RH_cond = np.max(abs(Nk_face_LLF[:,:,1] - Nk_face_LLF[:,:,0])>1,axis=0)
        alpha_LLF[alpha_RH_cond] = np.max(abs(alpha_RH[:,alpha_RH_cond])[:,np.newaxis],axis=0).T
        '''
        Fk_face_LLF = 0.5*(Fk_face_LLF_all.sum(axis=-1) - np.max(abs(alpha_LLF),axis=-1) * \
                    (Nk_face_LLF[:,:,1] - Nk_face_LLF[:,:,0]))
        return Fk_face_LLF

    def update_flux_MDW(self, Fk_face, Nk_face, alpha_MDW):
        alpha = (abs(alpha_MDW))
        Fk_face_MDW = 0.5*(Fk_face.sum(axis=-1) - alpha * \
                    (Nk_face[:,:,1] - Nk_face[:,:,0]))
        return Fk_face_MDW

    def update_flux_DW(self, Fk_face, Nk_face, alpha_DW):
        Fk_face_DW = 0.5*(Fk_face.sum(axis=-1) - abs(alpha_DW) * \
                    (Nk_face[:,:,1] - Nk_face[:,:,0]))
        return Fk_face_DW

    def update_flux_ROE(self, Fk_face, Nk_face, dFkdNk_eigval_m, dFkdNk_eigvec_m):
        A = np.identity(ctes.n_components)[:,:,np.newaxis]
        Gamma = abs(dFkdNk_eigval_m[:,np.newaxis,:] * A)
        R = dFkdNk_eigvec_m
        R_reshaped = R.transpose(2,0,1)
        Gamma_reshaped = Gamma.transpose(2,0,1)
        alpha_ROE = R_reshaped@Gamma_reshaped@np.linalg.inv(R_reshaped)
        dNk = (Nk_face[:,:,1] - Nk_face[:,:,0])
        product = (alpha_ROE@dNk.T[:,:,np.newaxis])[...,0]
        Fk_face_ROE = 0.5*(Fk_face.sum(axis=-1) - product.T)
        return Fk_face_ROE

    def Harten_entropy_corr(self, alpha, alpha_LR):
        e0 = 2
        alpha_ = np.copy(alpha)
        alpha_L = np.max(abs(alpha_LR[...,0]),axis=0)
        alpha_R = np.max(abs(alpha_LR[...,1]),axis=0)
        e = e0 * (alpha_R - alpha_L)
        alpha_[abs(alpha) < e] = ((alpha*alpha + e*e)/2*e)[abs(alpha) < e]
        return alpha_

class FirstOrder:
    def __init__(self):
        pass

    def FOU(self, M, fprop, ftotal, mobilities_internal_faces):
        UPW = Flux()
        Fk_vols_total = UPW.update_flux(M, fprop, ftotal,
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
            ftotal, ponteiro)
        #wave_velocity[:,ponteiro] = 1e-10
        #Fk_face = RS.get_Fk_face(fprop, M, Nk_face, fprop.P[ctes.v0], ftotal)
        #wave_velocity_RH = (Fk_face[...,1] - Fk_face[...,0])/(Nk_face[...,1] - Nk_face[...,0])
        #e = 1e-5
        #wave_velocity[abs(Nk_face[...,1] - Nk_face[...,0])>e] = wave_velocity_RH[abs(Nk_face[...,1] - Nk_face[...,0])>e]
        #import pdb; pdb.set_trace()
        #wave_velocity = Fk_vols_total/fprop.Nk #np.max(abs(wave_velocity),axis=0)
        return Fk_vols_total, wave_velocity

    def LLF(self, M, fprop, ftotal, P_old, Nk_old):

        Nk_face = Nk_old[:,ctes.v0]
        P_face = P_old[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)
        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)

        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, P_face, ftotal)
        Fk_internal_faces, wave_velocity = RS.LLF(M, fprop, Nk_face, P_face,
            ftotal, Fk_face, np.ones_like(ftotal[0],dtype=bool))
        '''wave_velocity_LR = wave_velocity_5[...,:2]
        upw = (wave_velocity_LR[...,0] * wave_velocity_LR[...,1])<=0
        upw[(wave_velocity_LR[...,0]==0) * (wave_velocity_LR[...,1]==0)] = False
        upw_cond = np.sum(upw, axis=0, dtype=bool)'''
        # contour faces upwind clássico (organize this later)
        Fk_internal_faces[:,[0,1]] = Fk_face[:,[0,1],0]

        Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)
        return Fk_vols_total, wave_velocity

    def MDW(self, M, fprop, ftotal, P_old):

        Nk_face = fprop.Nk[:,ctes.v0]
        P_face = P_old[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)

        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)

        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, P_face, ftotal)
        ponteiro = np.zeros(ctes.n_internal_faces,dtype=bool)
        Fk_internal_faces, wave_velocity = RS.MDW(M, fprop, Nk_face, P_face,
            ftotal, Fk_face, ~ponteiro)

        # contour faces upwind clássico (organize this later)
        Fk_internal_faces[...,0] = Fk_face[:,0,0]
        #Fk_internal_faces[...,2] = Fk_face[:,2,0]
        #Fk_internal_faces[...,1] = Fk_face[:,1,0]

        Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)
        return Fk_vols_total, wave_velocity

    def DW(self, M, fprop, ftotal, P_old):

        Nk_face = fprop.Nk[:,ctes.v0]
        P_face = P_old[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)

        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)

        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, P_face, ftotal)
        ponteiro = np.zeros(ctes.n_internal_faces,dtype=bool)
        Fk_internal_faces, wave_velocity = RS.DW(M, fprop, Nk_face, P_face,
            ftotal, Fk_face)

        # contour faces upwind clássico (organize this later)
        Fk_internal_faces[...,0] = Fk_face[:,0,0]
        #Fk_internal_faces[...,1] = Fk_face[:,1,0]

        Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)
        return Fk_vols_total, wave_velocity

    def ROE(self, M, fprop, ftotal, P_old):
        Nk_face = fprop.Nk[:,ctes.v0]
        P_face = P_old[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis], P_face[:,np.newaxis]),axis=1)
        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)

        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, P_face, ftotal)
        Fk_internal_faces, wave_velocity = RS.ROE_MDW(M, fprop, Nk_face, P_face,
            ftotal, Fk_face)

        # contour faces upwind clássico (organize this later)
        Fk_internal_faces[...,0] = Fk_face[:,0,0]
        #Fk_internal_faces[...,1] = Fk_face[:,1,0]

        Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)
        return Fk_vols_total, wave_velocity

class MUSCL:

    """ Class created for the second order MUSCL implementation for the \
    calculation of the advective terms """

    def run(self, M, fprop, wells, P_old, ftot, Pot_hid):

        ''' Global function that calls others '''
        self.P_face = np.sum(P_old[ctes.v0], axis=1) * 0.5
        self.P_face = np.concatenate((self.P_face[:,np.newaxis], self.P_face[:,np.newaxis]),axis=1)
        dNk_vols = self.volume_gradient_reconstruction(M, fprop, wells)
        #self.P_face[ctes_MUSCL.faces_contour] = P_old[ctes.v0[ctes_MUSCL.faces_contour,0]]
        dNk_face, dNk_face_neig = self.get_faces_gradient(M, fprop, dNk_vols)

        Phi = self.Van_Leer_slope_limiter(dNk_face, dNk_face_neig)

        Nk_face, z_face = self.get_extrapolated_compositions(fprop, Phi, dNk_face_neig)

        #G = self.update_gravity_term() # for now, it has no gravity
        alpha, Fk_vols_total = self.update_flux(M, wells, fprop, Nk_face,
            ftot, Pot_hid)
        #alpha = fprop.Fk_vols_total/fprop.Nk
        #import pdb; pdb.set_trace()
        return alpha, Fk_vols_total

    def volume_gradient_reconstruction(self, M, fprop, wells):
        Nk_neig =  fprop.Nk[:,np.newaxis,:] * ctes_MUSCL.allneig_and_vol[np.newaxis,:,:]
        Nk = Nk_neig.transpose(0,2,1)

        dNk = Nk_neig - Nk

        dNk_by_axes = np.repeat(dNk[:,np.newaxis,:,:],3,axis=1)
        dNk_by_axes *= ctes_MUSCL.versor_ds
        dNk_vols = dNk_by_axes.sum(axis = 3)

        dNkds_vols = np.copy(dNk_vols)
        dNkds_vols[:,ctes_MUSCL.ds_vols!=0] = dNk_vols[:,ctes_MUSCL.ds_vols != 0] / \
            ctes_MUSCL.ds_vols[ctes_MUSCL.ds_vols != 0][np.newaxis,:]

        dNkds_vols[:,ctes_MUSCL.all_neig_by_axes==1] = 0 #*dNkds_vols[:,:,all_neig.sum(axis=1)==1]
        #dNkds_vols[:,:,ctes.v0[ctes_MUSCL.faces_contour].flatten()] = 0 # zero in the contour volumes
        return dNkds_vols

    def get_faces_gradient(self, M, fprop, dNkds_vols):
        #ds_face = (M.data['centroid_volumes'][ctes.v0[:,1],:] -  M.data['centroid_volumes'][ctes.v0[:,0],:])
        #if any(ctes_MUSCL.ds_face.flatten()<0): import pdb; pdb.set_trace()
        #versor_ds_face = np.sum(ctes_MUSCL.versor_ds,axis=0)
        #ds_face_abs = abs(ds_face)

        dNk_face = (fprop.Nk[:,ctes.v0[:,1]] - fprop.Nk[:,ctes.v0[:,0]]) #* ctes_MUSCL.versor_ds_face[np.newaxis,:]
        dNk_face_vols = 2. * (dNkds_vols[:,:,ctes.v0] *ctes_MUSCL.ds_face.T[np.newaxis,:,:,np.newaxis]).sum(axis=1)
        dNk_face_neig = dNk_face_vols - dNk_face[:,:,np.newaxis]
        #dNk_face2 = dNk_face[...,np.newaxis]  * np.ones_like(dNk_face_neig)
        #dNk_face_vols[abs(dNk_face_neig)<1e-25] = dNk_face2[abs(dNk_face_neig)<1e-25]
        return dNk_face, dNk_face_neig

    def Van_Leer_slope_limiter(self, dNk_face, dNk_face_neig):
        np.seterr(divide='ignore', invalid='ignore')
        r_face = dNk_face[:,:,np.newaxis] / dNk_face_neig
        r_face[dNk_face_neig==0] = 0
        phi = (r_face + abs(r_face)) / (r_face + 1)
        phi[r_face<0] = 0 #so botei pra caso r==-1
        Phi = phi
        Phi[:,:,1] = -Phi[:,:,1]
        return Phi

    def Van_Albada1_slope_limiter(self, dNk_face, dNk_face_neig):
        np.seterr(divide='ignore', invalid='ignore')
        r_face = dNk_face[:,:,np.newaxis] / dNk_face_neig
        r_face[dNk_face_neig==0] = 0
        phi = (r_face**2 + abs(r_face)) / (r_face**2 + 1)
        phi[r_face<0]=0 #so botei pra caso r==-1
        Phi = phi
        Phi[:,:,1] = -Phi[:,:,1]
        return Phi

    def get_extrapolated_compositions(self, fprop, Phi, dNk_face_neig):
        Phi[abs(dNk_face_neig)<1e-30] = 0
        Nk_face = fprop.Nk[:,ctes.v0] + Phi / 2 * dNk_face_neig
        #Nk_face[(abs(Nk_face)<1e-30)] = fprop.Nk[:,ctes.v0][(abs(Nk_face)<1e-30)]
        if any(Nk_face.flatten()<0): import pdb; pdb.set_trace()
        z_face = Nk_face[0:ctes.Nc] / np.sum(Nk_face[0:ctes.Nc], axis = 0)
        return Nk_face, z_face

    def update_gravity_term(self):
        G = ctes.g * self.rho_j_face * ctes.z[ctes.v0]
        return G

    '''def flux_calculation_conditions_Serna(self, alpha, d2FkdNk):
        #ponteiro_LLF = np.ones(ctes.n_internal_faces,dtype=bool)
        #import pdb; pdb.set_trace()
        ponteiro_LLF = np.ones((ctes.n_components,ctes.n_internal_faces),dtype=bool)
        ponteiro_LLF[alpha[:,:,0] * alpha[:,:,1] <= 0] = False
        ponteiro_LLF[d2FkdNk[:,:,0] * d2FkdNk[:,:,0] <= 0] = False
        ponteiro_LLF = ponteiro_LLF.sum(axis=0,dtype=bool)
        return ponteiro_LLF'''

    def update_flux_upwind(self, Pot_hid, Fk_face_upwind_all, ponteiro):
        Fk_face_upwind = np.empty_like(Fk_face_upwind_all[:,:,0])

        Pot_hidj = Pot_hid[0,ctes.v0[:,0]][ponteiro] #- G[0,:,:,0]
        Pot_hidj_up = Pot_hid[0,ctes.v0[:,1]][ponteiro] #- G[0,:,:,1]

        Fk_face_upwind[:,Pot_hidj_up <= Pot_hidj] = \
            Fk_face_upwind_all[:,Pot_hidj_up <= Pot_hidj, 0]
        Fk_face_upwind[:,Pot_hidj_up > Pot_hidj] = \
            Fk_face_upwind_all[:,Pot_hidj_up > Pot_hidj, 1]
        return Fk_face_upwind

    def update_flux(self, M, wells, fprop, Nk_face, ftotal, Pot_hid):
        Fk_internal_faces = np.zeros((ctes.n_components,ctes.n_internal_faces))
        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)
        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, self.P_face, ftotal)
        ponteiro = np.zeros(ctes.n_internal_faces,dtype=bool)
        'Corrigir o cálculo do fluxo por upw - aplicar antese Fk_face vem um só'
        if ctes.RS['LLF']:
            alpha_wv = np.zeros((ctes.n_internal_faces, 5))
            Fk_internal_faces[:,~ponteiro], alpha_wv[~ponteiro,:] = RS.LLF(M, \
                fprop, Nk_face, self.P_face, ftotal, Fk_face, ~ponteiro)

        elif ctes.RS['MDW']:
            alpha_wv = np.zeros((ctes.n_internal_faces))
            Fk_internal_faces[:,~ponteiro], alpha_wv[~ponteiro] = RS.MDW(M, \
                fprop, Nk_face, self.P_face, ftotal, Fk_face, ~ponteiro)

        elif ctes.RS['DW']:
            alpha_wv = np.zeros((ctes.n_internal_faces))
            Fk_internal_faces[:,~ponteiro], alpha_wv[~ponteiro] = RS.DW(M, \
                fprop, Nk_face, self.P_face, ftotal, Fk_face, ~ponteiro)

        elif ctes.RS['ROE']:
            alpha_wv = np.zeros((ctes.n_internal_faces))
            Fk_internal_faces[:,~ponteiro], alpha_wv[~ponteiro] = RS.ROE_LLF(M, \
                fprop, Nk_face, self.P_face, ftotal, Fk_face)
        else:
            alpha_wv = np.zeros((ctes.n_components, ctes.n_internal_faces))
            ponteiro_alpha = np.zeros_like(ponteiro)
            #dNkmax_small = np.max(abs(Nk_face[:,:,0]-Nk_face[:,:,1]),axis=0)<1e-15
            #ponteiro_alpha[dNkmax_small] = False
            #alpha_wv[:,~ponteiro_alpha], m = RS.medium_wave_velocity(M, fprop, Nk_face, self.P_face, \
            #    ftotal, ~ponteiro_alpha)
            alpha_wv[:,~ponteiro_alpha] = 1e-100
            Fk_internal_faces[:,~ponteiro] = self.update_flux_upwind(Pot_hid, \
                Fk_face, ~ponteiro)
        #import pdb; pdb.set_trace()

        ponteiro[ctes_MUSCL.faces_contour] = True
        Fk_internal_faces[:,ponteiro] = self.update_flux_upwind(fprop.P[np.newaxis,:], \
            Fk_face[:,ponteiro], ponteiro)

        '-------- Perform volume balance to obtain flux through volumes -------'
        Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)
        if any(np.isnan(Fk_vols_total).flatten()): import pdb; pdb.set_trace()
        #Fk_vols_total[:ctes.Nc][fprop.z==0] = 0
        if any(Fk_vols_total[:ctes.Nc][fprop.z==0]<0): import pdb; pdb.set_trace()
        return alpha_wv, Fk_vols_total

class FR:

    def __init__(self):
        'Enviroment for the FR/CPR method - vai ser 1D por enquanto'
        'OBS: O Fk que vai entrar aqui provavelmente, se não interpretei errado \
        tem que ser em mol*m/s, ou seja, a pretransmissibilidade n pode ter dx \
        dividindo'
        '1. Obtain Nk at the SP'
        '2. Compute Fk with the mobility, xkj, rho and csi approximation from \
        volume-weigthed arithmetic average'
        '3. Transform Fk for the SP by using RTo'
        '4. Use a Riemann Solver for computing flux at the interfaces, and it will \
        be used for the continuous flux approximation (This can be done previously'
        '5. Approximate Fk and Nk in the reference domain by using Lagranges \
        polynomial. Where Fk = Fk^D + Fk^C, where Fk^D and Fk^C are also obtained \
        by using Lagranges polynomial with its value from the SP.'
        '6. Obtain Nk for the next time step by the third order RungeKutta'
        '7. Project the Nk solution at the Pspace back to the CV using Gaussian \
        integration'
        'Legend: n_points - number of SP per control volume'

    def run(self, M, fprop, wells, Ft_internal_faces, Nk_SP_old, P_old, q, delta_t, t):
        Nk_SP = np.copy(Nk_SP_old)

        q_SP = q[:,:,np.newaxis] * np.ones_like(Nk_SP)
        self.P_faces = np.sum(P_old[ctes_FR.v0],axis=-1)*0.5
        self.P_SP = self.get_pressure_SP(wells, P_old)
        self.Nk_SP_old = Nk_SP_old
        self.delta_t = delta_t

        self.Ft_SP = self.total_flux_SP(fprop, wells, Ft_internal_faces)

        dFk_SP, wave_velocity = self.dFk_SP_from_Pspace(M, fprop, wells, Ft_internal_faces, np.copy(Nk_SP), P_old)
        Nk_SP = RK3.update_composition_RK3_1(np.copy(Nk_SP_old), q_SP, dFk_SP, delta_t)
        Nk_SP = self.MLP_slope_limiter(M, fprop, Nk_SP, wells)

        '''dFk_SP, wave_velocity = self.dFk_SP_from_Pspace(M, fprop, wells, Ft_internal_faces, np.copy(Nk_SP), P_old)
        Nk_SP = RK3.update_composition_RK3_2(np.copy(Nk_SP_old), q_SP, np.copy(Nk_SP), dFk_SP, delta_t)
        Nk_SP = self.MLP_slope_limiter(M, fprop, Nk_SP, wells)

        dFk_SP, wave_velocity = self.dFk_SP_from_Pspace(M, fprop, wells, Ft_internal_faces, np.copy(Nk_SP), P_old)
        Nk_SP = RK3.update_composition_RK3_3(np.copy(Nk_SP_old), q_SP, np.copy(Nk_SP), dFk_SP, delta_t)
        Nk_SP = self.MLP_slope_limiter(M, fprop, Nk_SP, wells)'''

        Fk_vols_total = np.min(abs(dFk_SP),axis=2)
        Nk = 1 / sum(ctes_FR.weights) * np.sum(ctes_FR.weights * Nk_SP,axis=2)

        z = Nk[0:ctes.Nc,:] / np.sum(Nk[0:ctes.Nc,:], axis = 0)
        if (any((z<0).flatten())): import pdb; pdb.set_trace()
        return wave_velocity, Nk, z, Nk_SP, Fk_vols_total

    def dFk_SP_from_Pspace(self, M, fprop, wells, Ft_internal_faces, Nk_SP, P_old):

        Fk_SP = self.component_flux_SP(fprop, M, Nk_SP)
        #Fk_SP = self.get_Fk_Ft_SP(fprop, M, Nk_SP)

        Fk_faces, Fk_vols_RS_neig, wave_velocity = self.Riemann_Solver(M, fprop, wells, Nk_SP,
            Fk_SP, Ft_internal_faces)

        Fk_D = np.sum(Fk_SP[:,:,:,np.newaxis] * ctes_FR.L[np.newaxis,np.newaxis,:], axis=2)
        dFk_D = np.sum(Fk_SP[:,:,:,np.newaxis] * ctes_FR.dL[np.newaxis,np.newaxis,:], axis=2)
        dFk_C = self.dFlux_Continuous(Fk_SP, Fk_vols_RS_neig)
        dFk_Pspace = (dFk_C + dFk_D)

        #Fk_vols_RS_neig[:,wells['all_wells'],0] = -Fk_vols_RS_neig[:,wells['all_wells'],0]
        if ctes.load_w and not data_loaded['compositional_data']['water_data']['mobility']:
            dFk_Pspace[-1,:] = 0

        #dFk_Pspace[:,wells['all_wells'],1:] = 0
        #dFk_Pspace[:,wells['all_wells'],0] = (Fk_vols_RS_neig[:,wells['all_wells']]).sum(axis=2)/2
        dFk_SP = dFk_Pspace @ ctes_FR.x_points
        dFk_SP = - 2 * dFk_SP #this way only works for uniform mesh
        #up! transforming from local space to original global space (this could be done to the g and L functions
        #only, however, I rather do like this, so it's done just once)
        #if any(dFk_SP[:,0,:].flatten()>0): import pdb; pdb.set_trace()
        #dF = (dFk_SP[:,[0,-1],0] - dFk_SP[:,[0,-1],1]) +  (dFk_SP[:,[0,-1],1] - dFk_SP[:,[0,-1],2])
        #if any(abs(dF).flatten()>1e-15): import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()
        '''self.dFk_C = -2*dFk_C@ ctes_FR.x_points
        self.dFk_D = -2*dFk_D@ ctes_FR.x_points
        self.Fk_SP = Fk_SP'''
        return dFk_SP, wave_velocity

    def total_flux_SP(self, fprop, wells, Ft_internal_faces):
        'RTo'
        phi = np.empty((len(ctes_FR.points),2))
        phi[:,0] = 1 / 4 * (1 + ctes_FR.points)
        phi[:,1] = 1 / 4 * (1 - ctes_FR.points)

        Ft_face_phi = (Ft_internal_faces[:,:,np.newaxis,np.newaxis] * phi[np.newaxis,np.newaxis,:])

        'Look for a faster way to do that'
        Ft_SP_reshaped = np.empty((1,ctes.n_volumes,ctes_FR.n_points))
        contours = np.array([0,ctes_FR.n_points-1])
        for i in range(ctes_FR.n_points):
            lines = np.array([np.zeros_like(ctes_FR.v0[:,0]), np.zeros_like(ctes_FR.v0[:,1])]).astype(int).flatten()
            cols = np.array([ctes_FR.v0[:,0], ctes_FR.v0[:,1]]).flatten()
            data = np.array([Ft_face_phi[:,:,i,0], Ft_face_phi[:,:,i,1]]).flatten()
            Ft_SP_reshaped[:,:,i] = sp.csc_matrix((data, (lines, cols)), shape = (1, ctes.n_volumes)).toarray()
        Ft_SP = 2 * Ft_SP_reshaped #np.concatenate(np.dsplit(Ft_SP_reshaped, ctes_FR.n_points), axis = 2)

        '''if ctes_FR.n_points==3:
            P_f = fprop.P[ctes.v0].sum(axis=-1)/2
            Pot_hid = P_f[ctes_FR.vols_vec]
            Pot_hidj = Pot_hid[:,0]
            Pot_hidj_up = Pot_hid[:,-1]
            z = ctes.z
            z_up = ctes.z
            K_vols = ctes.pretransmissibility_internal_faces[ctes_FR.vols_vec[:,0]]
            #mob_1 = fprop.mobilities_internal_faces[...,ctes_FR.vols_vec].sum(axis=-1)/2
            xx = - np.sum(fprop.mobilities
                * K_vols * ((Pot_hidj_up - Pot_hidj) -
                ctes.g * fprop.rho_j * (z_up - z)), axis = 1)
            Ft_SP[:,1:-1,1] = xx[:,1:-1]'''

        #Ft_SP[0,wells['all_wells'],:] = ((Ft_internal_faces[0,ctes_FR.vols_vec][wells['all_wells']]).sum(axis=-1)/2)[:,np.newaxis]
        #Ft_SP[:,wells['all_wells']] = (Ft_SP_reshaped[:,wells['all_wells']]).sum(axis=-1)[:,:,np.newaxis]
        return Ft_SP

    def get_pressure_SP(self, wells, Pold):
        #P_SP = np.empty((ctes.n_volumes,ctes_FR.n_points))
        #P_SP[:,[0,-1]] = self.P_faces[ctes_FR.vols_vec]
        x_0 = np.copy(ctes_FR.points)
        x_1 = np.copy(ctes_FR.points)
        x_0[x_0>0] = 0
        x_1[x_1<0] = 0
        P_SP = Pold[:, np.newaxis]  + (Pold - self.P_faces[ctes_FR.vols_vec[:,0]])[:,np.newaxis] * x_0 - \
            (Pold - self.P_faces[ctes_FR.vols_vec[:,1]])[:,np.newaxis] * x_1
        #P_SP[wells['all_wells'],:] = Pold[wells['all_wells']][...,np.newaxis] - \
        #    (Pold[wells['all_wells']]-P_SP[wells['all_wells'],0])[...,np.newaxis] \
        #    * ctes_FR.points
        #import pdb; pdb.set_trace()
        #P_SP[:,:] = Pold[:,np.newaxis]
        #P_SP[wells['all_wells'],:] = P_SP[wells['all_wells'],0][...,np.newaxis]
        return P_SP

    def component_flux_SP_inputs(self, fprop):
        Fk_SP_inputs = dict()
        Fk_SP_inputs['ponteiro'] = np.ones(ctes.n_volumes,dtype=bool)
        Fk_SP_inputs['v0'] = np.arange(ctes.n_volumes)[:,np.newaxis] * np.ones((ctes.n_volumes,ctes_FR.n_points))
        Fk_SP_inputs['Ft_SP_flatt'] = np.concatenate(np.dsplit(self.Ft_SP, ctes_FR.n_points),axis=1)[:,:,0]
        Fk_SP_inputs['pretr'] = ctes.pretransmissibility_internal_faces[ctes_FR.vols_vec][:,0]
        Fk_SP_inputs['Vp_SP'] = np.tile(fprop.Vp, ctes_FR.n_points)
        Fk_SP_inputs['P_SP_flatt'] = np.concatenate(np.hsplit(self.P_SP, ctes_FR.n_points),axis=0)[:,0]
        return Fk_SP_inputs

    def component_flux_SP(self, fprop, M, Nk_SP):
        Fk_SP_inputs = self.component_flux_SP_inputs(fprop)
        Nk_SP_flatt = np.concatenate(np.dsplit(Nk_SP, ctes_FR.n_points),axis=1)[:,:,0]
        Fk_SP = RiemannSolvers(Fk_SP_inputs['v0'].astype(int), Fk_SP_inputs['pretr']).Fk_from_Nk(fprop,
            M, Nk_SP_flatt, Fk_SP_inputs['P_SP_flatt'], Fk_SP_inputs['Vp_SP'],
            Fk_SP_inputs['Ft_SP_flatt'], Fk_SP_inputs['ponteiro'])
        Fk_SP = np.concatenate(np.hsplit(Fk_SP[:,:,np.newaxis], ctes_FR.n_points), axis = 2)
        return Fk_SP

    def get_Fk_Ft_SP(self, fprop, M, Nk_SP):
        v = ctes_FR.n_points
        Fk_SP_inputs = self.component_flux_SP_inputs(fprop)
        RS = RiemannSolvers(Fk_SP_inputs['v0'].astype(int), Fk_SP_inputs['pretr'])
        z_SP = Nk_SP[0:ctes.Nc] / np.sum(Nk_SP[0:ctes.Nc], axis = 0)
        z_SP_flatt = np.concatenate(np.dsplit(z_SP, ctes_FR.n_points),axis=1)[:,:,0]
        Nk_SP_flatt = np.concatenate(np.dsplit(Nk_SP, ctes_FR.n_points),axis=1)[:,:,0]

        mobilities_SP_flatt, rho_j_SP, Csi_j_SP, xkj_SP = RS.get_extrapolated_properties(fprop, M, Nk_SP_flatt, z_SP_flatt,
            Fk_SP_inputs['P_SP_flatt'], Fk_SP_inputs['Vp_SP'], ctes_FR.n_points, Fk_SP_inputs['ponteiro'])
        mobilities_SP = np.concatenate(np.split(mobilities_SP_flatt[...,np.newaxis], ctes_FR.n_points, axis=-2),axis=3)

        P_f = fprop.P[ctes.v0].sum(axis=-1)/2
        dP_SP = np.empty((ctes.n_volumes,ctes_FR.n_points))
        dP_SP[:,[0,-1]] = (fprop.P[ctes_FR.v0[:,1]] - fprop.P[ctes_FR.v0[:,0]])[ctes_FR.vols_vec]
        if ctes_FR.n_points==3:
            dP_SP[:,1] = P_f[ctes_FR.vols_vec[:,1]] - P_f[ctes_FR.vols_vec[:,0]]
            dP_SP[0,1] = dP_SP[0,-1]
            dP_SP[-1,1] = dP_SP[-1,-1]
        K_vols = ctes.pretransmissibility_internal_faces[ctes_FR.vols_vec[:,0]][np.newaxis,np.newaxis,:,np.newaxis]

        Ft_SP = - np.sum(mobilities_SP * K_vols * dP_SP[np.newaxis,np.newaxis,:], axis = 1)
        Ft_SP_flatt = np.concatenate(np.dsplit(self.Ft_SP, ctes_FR.n_points),axis=1)[:,:,0]

        f = Flux()

        Pcap_reshaped = np.concatenate(np.dsplit(np.tile(fprop.Pcap[:,Fk_SP_inputs['v0'].astype(int)][:,Fk_SP_inputs['ponteiro']],v),v),axis=1)
        z_reshaped = np.concatenate(np.hsplit(np.tile(ctes.z[Fk_SP_inputs['v0'].astype(int)][Fk_SP_inputs['ponteiro']],v),v),axis=0)
        Fj = f.update_Fj_internal_faces(Ft_SP_flatt, rho_j_SP, mobilities_SP_flatt, Pcap_reshaped,
            z_reshaped, np.tile(Fk_SP_inputs['pretr'][Fk_SP_inputs['ponteiro']],v))
        Fk = f.update_Fk_internal_faces(xkj_SP, Csi_j_SP, Fj)
        Fk_SP = np.concatenate(np.hsplit(Fk[:,:,np.newaxis], ctes_FR.n_points), axis = 2)
        return Fk_SP

    def Riemann_Solver(self, M, fprop, wells, Nk_SP, Fk_SP, Ft_internal_faces):
        Nk_faces = np.empty((ctes.n_components, ctes.n_internal_faces, 2))
        Nk_faces[:,:,1] = Nk_SP[:,ctes_FR.v0[:,1],0] #Nk faces a esquerda dos volumes
        Nk_faces[:,:,0] = Nk_SP[:,ctes_FR.v0[:,0],-1] #Nk nas faces a direita

        P_face = np.concatenate((self.P_faces[:,np.newaxis],self.P_faces[:,np.newaxis]),axis=1)

        Fk_faces = np.empty_like(Nk_faces)
        Fk_faces[:,:,1] = Fk_SP[:,ctes_FR.v0[:,1],0]
        Fk_faces[:,:,0] = Fk_SP[:,ctes_FR.v0[:,0],-1]

        '''Nk_face_contour = np.empty((ctes.n_components,1,2))
        Nk_face_contour[:,0,1] = Nk_SP[:,0,0]
        Nk_face_contour[:,0,0] = Nk_SP[:,-1,-1]
        Nk_face_contour[-1,0,0] = 1*fprop.Csi_j[0,-1,0]*fprop.Vp[0]

        Fk_faces_contour = RiemannSolvers(np.array([0,ctes.n_volumes-1])[np.newaxis]).get_Fk_face(fprop, M,
            Nk_face_contour, P_face[np.newaxis,0], Ft_internal_faces[:,0][:,np.newaxis])
        Fk_face_contour_RS, alpha_wv =  RiemannSolvers(np.array([0,ctes.n_volumes-1])[np.newaxis]).LLF(M, fprop, Nk_face_contour, P_face[np.newaxis,0],
            Ft_internal_faces[:,0][:,np.newaxis], Fk_faces_contour, np.zeros(1,dtype=bool))'''

        Fk_face_RS, alpha_wv =  RiemannSolvers(ctes_FR.v0, ctes.pretransmissibility_internal_faces).LLF(M, fprop, Nk_faces, P_face,
            Ft_internal_faces, Fk_faces, np.ones_like(Ft_internal_faces[0],dtype=bool))
        #ponteiro = np.zeros_like(Ft_internal_faces[0], dtype=bool)
        #ponteiro[[0,1]] = True
        #Fk_face_RS[:,ponteiro] = MUSCL().update_flux_upwind(fprop.P[np.newaxis,:], Fk_faces[:,ponteiro], ponteiro)
        Fk_face_RS[:,0] = Fk_faces[:,0,0]
        Fk_face_RS[:,1] = Fk_faces[:,1,0]

        'Obtaining Flux at each CV side - by finding faces that compounds the CV \
        this only works for 1D problems'
        self.vols_vec = -np.ones((ctes.n_volumes,2),dtype=int)
        lines = np.arange(ctes.n_internal_faces)
        self.vols_vec[ctes_FR.v0[:,0],1] = lines
        self.vols_vec[ctes_FR.v0[:,1],0] = lines
        Fk_vols_RS_neig = Fk_face_RS[:,ctes_FR.vols_vec]
        Fk_vols_RS_neig[:,self.vols_vec<0] = 0 #Fk_face_contour_RS
        #Fk_vols_RS_neig[:,vols_vec<0] = ((Fk_vols_RS_neig[:,wells['all_wells']]).sum(axis=2))
        return Fk_faces, Fk_vols_RS_neig, alpha_wv

    def dFlux_Continuous(self, Fk_SP, Fk_vols_RS_neig):
        Fk_D_l = Fk_SP[...,0]#Fk_D @ x_left
        Fk_D_r = Fk_SP[...,-1]#Fk_D @ x_right
        dFk_C = (Fk_vols_RS_neig[:,:,0] - Fk_D_l)[:,:,np.newaxis] * ctes_FR.dgLB[np.newaxis,:] + \
                (Fk_vols_RS_neig[:,:,1] - Fk_D_r)[:,:,np.newaxis] * ctes_FR.dgRB[np.newaxis,:]
        return dFk_C

    def MLP_slope_limiter(self, M, fprop, Nk_SP_in, wells):

        inds = np.array([0,-1])

        Nk = (np.linalg.inv(ctes_FR.V)[np.newaxis,] @ Nk_SP_in[:,:,:, np.newaxis])[:,:,:,0]
        Nk_SP = self.projections(Nk, ctes_FR.n_points-1)

        #self.machine_error = np.max(abs(Nk_SP - Nk_SP_in))
        #if machine_error<np.finfo(np.float64).eps: machine_error = np.finfo(np.float64).eps
        self.machine_error = 1e-10 #323 #1e-150 #1e-323 #1e-700 #1e-150

        inds = np.array([0,-1])

        'Projected n=0'
        Nk_P0 = self.projections(Nk, 0)
        Nk_P0_vertex = Nk_P0[:,:,inds]

        'Projected n=1'
        Nk_P1 = self.projections(Nk, 1)
        Nk_P1_vertex = Nk_P1[:,:,inds]

        'Neigboring vertex points values'
        Nk_faces = np.empty((ctes.n_components,ctes.n_internal_faces,2))

        Nk_faces[:,:,1] = Nk_SP[:,ctes_FR.v0[:,1], 0]
        Nk_faces[:,:,0] = Nk_SP[:,ctes_FR.v0[:,0], -1]

        Nk_neig = np.empty((ctes.n_components,ctes.n_internal_faces,2))

        Nk_neig[:,:,0] = Nk_P0_vertex[:,ctes_FR.v0[:,0],0]
        Nk_neig[:,:,1] = Nk_P0_vertex[:,ctes_FR.v0[:,1],0]

        Phi_P1 = self.P1_limiter(M, Nk, Nk_SP, Nk_P0_vertex, Nk_P1_vertex, Nk_neig, \
            Nk_faces)

        phi_P2 = np.zeros_like(Phi_P1)
        phi_Pn = np.zeros_like(Phi_P1)
        Nk_P2 = Nk_SP

        if ctes_FR.n_points > 2:
            '---------------Hierarchical MLP limiting procedure----------------'

            #axis=1 due to reshaping of the argument [~phi_Pn.astype(bool)]
            phi_Pn = self.troubled_cell_marker(M, fprop, Nk_P1_vertex, Nk_P0_vertex,\
                Nk_SP, Nk_neig, Nk_faces)
            phi_P2 = phi_Pn

            if ctes_FR.n_points==4:

                'Projected n=2'
                Nk2 = np.zeros_like(Nk)
                Nk2[:,:,0:3] = Nk[:,:,0:3]
                Nk_P2 = self.projections(Nk, 2)
                Nk_P2_vertex = Nk_P2[:,:,inds]

                'Projected P2 into n=1'
                Nk_P12 = self.projections(Nk2, 1)
                Nk_P1_vertex2 = Nk_P12[:,:,inds]

                'Projected P2 into n=0'
                Nk_P02 = self.projections(Nk2, 0)
                Nk_P0_vertex2 = Nk_P02[:,:,inds]

                'Neigboring vertex points values'
                Nk_faces2 = np.empty((ctes.n_components,ctes.n_internal_faces,2))

                Nk_faces2[:,:,1] = Nk_P2[:,ctes_FR.v0[:,1], 0]
                Nk_faces2[:,:,0] = Nk_P2[:,ctes_FR.v0[:,0], -1]

                Nk_neig2 = np.empty((ctes.n_components,ctes.n_internal_faces,2))

                Nk_neig2[:,:,0] = Nk_P0_vertex2[:,ctes_FR.v0[:,0],0]
                Nk_neig2[:,:,1] = Nk_P0_vertex2[:,ctes_FR.v0[:,1],0]

                phi_P2 = self.troubled_cell_marker(M, fprop, Nk_P1_vertex2, \
                    Nk_P0_vertex2, Nk_P2, Nk_neig2, Nk_faces2)

            phi_P2[phi_Pn==1] = 1
            Phi_P1[phi_P2==1] = 1
            #phi_P2[np.sum(Nk_P2<0,axis=-1, dtype=bool)*np.sum(Nk_SP<0,axis=-1, dtype=bool)] = 0
            ##Phi_P1[np.sum(Nk_P1<0,axis=-1, dtype=bool)*(phi_P2==0)] = 0
        '''Phi_P1 = np.min(Phi_P1, axis=0)[np.newaxis,:]
        phi_P2 = np.min(phi_P2, axis=0)[np.newaxis,:]
        phi_Pn = np.min(phi_Pn, axis=0)[np.newaxis,:]'''


        Nk_SPlim = (Nk_P0 + Phi_P1[:,:,np.newaxis] * (Nk_P1 - Nk_P0) + \
            phi_P2[:,:,np.newaxis] * ((Nk_P2 - Nk_P1) + phi_Pn[:,:,np.newaxis] * (Nk_SP - Nk_P2)))

        Nk_SPlim[(Nk_SPlim<0)*(Nk_SPlim>-1e-323)] = 0

        if any(abs(Nk_SPlim[Nk_SPlim<0])):
            import pdb; pdb.set_trace()

        '''C_high_order_check1 = (Nk_SPlim[:,:,inds] <= np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]) #<= machine_error
        C_high_order_check2 = (Nk_SPlim[:,:,inds] >= np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]) #>= -machine_error
        C_high_order_check1[abs((Nk_SPlim[:,:,inds] - np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])) <= 1e-15] = True
        C_high_order_check2[abs((Nk_SPlim[:,:,inds] - np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])) <= 1e-15] = True

        high_order_check = C_high_order_check1 * C_high_order_check2
        #if any(np.min(1 * high_order_troubled_cells,axis=2)[phi_Pn==1]==0):import pdb; pdb.set_trace()
        phi_Pn_check = np.min(1 * high_order_check,axis = 2)#[phi_Pn==1]
        if len(phi_Pn_check[phi_Pn_check==0])>0: import pdb; pdb.set_trace()'''
        #import pdb; pdb.set_trace()
        return Nk_SPlim

    def projections(self, Nk, m):
        Nkm = np.zeros_like(Nk)
        Nkm[:,:,0:m+1] = Nk[:,:,0:m+1]
        Nk_Pm = (ctes_FR.V[np.newaxis,] @ Nkm[:,:,:,np.newaxis])[:,:,:,0]
        return Nk_Pm

    def P1_limiter(self, M, Nk, Nk_Pm, Nk_P0_vertex, Nk_P1_vertex, Nk_neig, Nk_faces):
        Linear_term = Nk_P1_vertex - Nk_P0_vertex

        'Neigboring vertex points values'

        #Phi_r = self.MLP_u1_mod(Nk_neig, Nk_P0_vertex, Linear_term)
        Phi_r = self.MLP_u1(Nk_neig, Nk_P0_vertex, Linear_term)
        #Phi_r = self.MLP_vk(M, Nk_neig, Nk_P0_vertex, Linear_term)

        Phi_P1 = np.ones_like(Phi_r)
        Phi_P1[abs(Nk_P1_vertex - Nk_P0_vertex) >= self.machine_error] = Phi_r[abs(Nk_P1_vertex - Nk_P0_vertex) >= self.machine_error]
        Phi_P1 = np.min(Phi_P1, axis = 2)
        #Phi_P1[:,-1] = 0
        #Phi_P1[:,0] = 0
        return Phi_P1

    def MLP_u1_mod(self, Nk_neig, Nk_P0_vertex, Linear_term):
        """ You and Kim paper page 27 of pdf"""
        Nk_avg_vertex = Nk_P0_vertex[:,ctes_FR.v0,0].sum(axis=-1)/2
        Nk_avg_vertex_vols = Nk_avg_vertex[:,ctes_FR.vols_vec]
        #Nk_avg_vertex_vols[:,self.vols_vec<0] = Nk_P0_vertex[:,self.vols_vec<0]
        f_min = Nk_avg_vertex_vols + np.heaviside(Nk_avg_vertex_vols - Nk_P0_vertex, np.ones_like(Nk_P0_vertex)) * \
                (np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec] - Nk_avg_vertex_vols)
        f_max = Nk_avg_vertex_vols + np.heaviside(Nk_avg_vertex_vols - Nk_P0_vertex, np.ones_like(Nk_P0_vertex)) * \
                (np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec] - Nk_avg_vertex_vols)

        r1 = (f_min - Nk_P0_vertex) / Linear_term
        r2 = (f_max - Nk_P0_vertex) / Linear_term
        rs = np.concatenate((r1[...,np.newaxis], r2[...,np.newaxis]),axis = -1)
        rs[Linear_term == 0] = 1
        r = np.max(rs, axis = -1)
        Phi_r_u1 = r
        Phi_r_u1[Phi_r_u1>1] = 1
        Phi_r_u1[Phi_r_u1<0] = 0
        return Phi_r_u1

    def MLP_u1(self, Nk_neig, Nk_P0_vertex, Linear_term):
        f_min = np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]
        f_max = np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]
        r1 = (f_min - Nk_P0_vertex) / Linear_term
        r2 = (f_max - Nk_P0_vertex) / Linear_term
        rs = np.concatenate((r1[...,np.newaxis], r2[...,np.newaxis]),axis = -1)
        rs[Linear_term == 0] = 1
        r = np.max(rs, axis = -1)
        Phi_r_u1 = r
        Phi_r_u1[Phi_r_u1>1] = 1
        Phi_r_u1[Phi_r_u1<0] = 0
        return Phi_r_u1

    def MLP_vk(self, M, Nk_neig, Nk_P0_vertex, Linear_term):
        delta_minus = Linear_term
        Nk_neig_max = np.max(Nk_neig,axis=2)[:,ctes_FR.vols_vec]
        Nk_neig_max[Nk_neig_max < Nk_P0_vertex] = Nk_P0_vertex[Nk_neig_max < Nk_P0_vertex]

        Nk_neig_min = np.min(Nk_neig,axis=2)[:,ctes_FR.vols_vec]
        Nk_neig_min[Nk_neig_min > Nk_P0_vertex] = Nk_P0_vertex[Nk_neig_min > Nk_P0_vertex]

        delta_plus = Nk_neig_min - Nk_P0_vertex
        delta_plus2 = Nk_neig_max - Nk_P0_vertex
        delta_plus[delta_minus>0] = delta_plus2[delta_minus>0]

        x_vols = M.data['centroid_volumes'][0,0]
        dx_vols = x_vols * 2
        K1 = 6 #1e-10 #BL use K1=2 or K1=6 for P3
        #K2 = 1e-14
        e2 = (K1 * dx_vols)**3

        #dNk_max_min = (np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec] - np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        #theta = dNk_max_min/(K2*dx_vols**1.5)
        #e2 = K1/(1+theta) * dNk_max_min**2

        Phi_r_u2 = 1/delta_minus * ((delta_plus * delta_plus + e2) * delta_minus
                + 2 * delta_minus * delta_minus * delta_plus) / \
                (delta_plus * delta_plus + 2 * delta_minus * delta_minus
                + delta_minus * delta_plus + e2)

        Phi_r_u2[delta_minus==0] = 1

        Phi_r_u2[Phi_r_u2 > 1] = 1
        Phi_r_u2[Phi_r_u2 < 0] = 0

        return Phi_r_u2

    def Phi_r_u2(self, M, r):
        Phi_r_u2 = (r**2 + 2*r + self.machine_error) / (r**2 + r + 2 + self.machine_error)
        return Phi_r_u2

    def P1_projected_cond(self, Nk_P1_vertex, Nk_neig):
        'P1-projected MLP condition - troubled-cell marker'

        C_troubled1 = (Nk_P1_vertex <= np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        C_troubled2 = (Nk_P1_vertex >= np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        #C_troubled1[abs((Nk_P1_vertex - np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])) <= self.machine_error] = True
        #C_troubled2[abs((Nk_P1_vertex - np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])) <= self.machine_error] = True
        troubled_cells = (C_troubled1 * C_troubled2)
        return troubled_cells

    def smooth_extrema_cond(self, fprop, Nk_P0_vertex, Nk_P1_vertex, Nk_Pmvertex, Nk_neig):
        'Smooth extrema detector'
        High_order_term = Nk_Pmvertex - Nk_P1_vertex
        Linear_term = Nk_P1_vertex - Nk_P0_vertex

        Nk_less_max = (Nk_Pmvertex < np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        Nk_big_min = (Nk_Pmvertex > np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        #Nk_less_max[abs(Nk_Pmvertex - np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]) <= self.machine_error] = True
        #Nk_big_min[abs(Nk_Pmvertex - np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]) <= self.machine_error] = True
        C1 = ((Linear_term) > 0) * ((High_order_term) < 0) * (Nk_big_min)
        C2 = ((Linear_term) < 0) * (High_order_term > 0) * (Nk_less_max)
        return C1, C2

    def Gibbs_Wilbraham_oscilations(self, M, Nk_faces):
        jump1 = Nk_faces[:,:,1] - Nk_faces[:,:,0]
        c_vols = M.data['centroid_volumes'][:,1]
        c_faces = c_vols[ctes.v0]
        hf = c_vols[0]*2 #abs(c_faces[:,1] - c_faces[:,0])

        smooth_bound_indicator = abs(jump1)/(abs(Nk_faces[:,:,1] + Nk_faces[:,:,0])/2) #em 1D apenas !!
        smooth_bound_indicator[np.isnan(smooth_bound_indicator)] = 0
        theta_f = smooth_bound_indicator/(hf**((ctes_FR.n_points)/2))
        trbl_bound_detector = np.empty_like(theta_f)

        trbl_bound_detector[theta_f<1] = 0 #normal bound
        trbl_bound_detector[theta_f>=1] = 1 #troubled bound
        trbl_bound_detector[:,0] = 1 #contour faces
        #trbl_bound_detector[:,1] = 1 #contour faces
        trbl_bound_detector_vols = np.empty_like(ctes_FR.vols_vec)
        trbl_bound_detector_vols = trbl_bound_detector[:,ctes_FR.vols_vec]
        #trbl_bound_detector_vols[:,[0,-1]] = 1 #Contour volumes
        Nf = np.sum(trbl_bound_detector_vols, axis = -1)
        return Nf

    def troubled_cell_marker(self, M, fprop, Nk_P1_vertex, Nk_P0_vertex, Nk_Pm, Nk_neig, Nk_faces):

        #Nk_faces_avg = 1/2 * np.sum(ctes_FR.weights[np.newaxis,np.newaxis,:,np.newaxis] * Nk_faces, axis=2)

        inds = np.array([0,-1])
        Nk_Pmvertex = Nk_Pm[:,:,inds]

        P1_proj_cond = self.P1_projected_cond(Nk_P1_vertex, Nk_neig)

        'Augmented MLP condition - troubled-cell marker'

        '''C_troubled1 = (Nk_Pmvertex <= np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        C_troubled2 = (Nk_Pmvertex >= np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        C_troubled1[abs((Nk_Pmvertex - np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])) <= machine_error] = True
        C_troubled2[abs((Nk_Pmvertex - np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])) <= machine_error] = True
        troubled_cells = (C_troubled1 * C_troubled2)
        phi_Pm =  np.min(1*(troubled_cells), axis = 2)'''

        C1, C2 = self.smooth_extrema_cond(fprop, Nk_P0_vertex, Nk_P1_vertex, Nk_Pmvertex, Nk_neig)

        #cC3 = np.concatenate((1e-3*abs(Nk_P0_vertex), fprop.Vp[np.newaxis,:,np.newaxis]*\
        #    np.ones_like(Nk_P0_vertex)),axis=-1)
        #C3 = (abs(Nk_Pmvertex - Nk_P0_vertex) <= np.max(cC3,axis=-1)[:,:,np.newaxis])

        #Nf = self.Gibbs_Wilbraham_oscilations(M, Nk_faces)

        '''A4-1 step'''
        P1_proj_cond_cells = np.min(1*(P1_proj_cond),axis=-1)

        'Type I trouble'

        '''aux_I = P1_proj_cond_cells[P1_proj_cond_cells==1]
        aux_Nf_I = Nf[P1_proj_cond_cells==1]
        aux_I[aux_Nf_I>=1] = 0
        P1_proj_cond_cells[P1_proj_cond_cells==1] = aux_I'''

        '''A4-2 step'''
        smooth_extrema = C1+C2#+C3
        smooth_extrema_cells = np.min(1*smooth_extrema,axis=-1)

        'Type II trouble'

        '''aux_II = smooth_extrema_cells[smooth_extrema_cells==1]
        aux_Nf_II = Nf[smooth_extrema_cells==1]
        aux_II[aux_Nf_II>=2] = 0
        smooth_extrema_cells[smooth_extrema_cells==1] = aux_II'''

        phi_Pm = 1*((P1_proj_cond_cells + smooth_extrema_cells).astype(bool))
        #phi_Pm = 1*((P1_proj_cond_cells).astype(bool))

        if any((phi_Pm[:,:,np.newaxis]*Nk_Pmvertex).flatten()<0): import pdb; pdb.set_trace()
        #phi_Pm = np.ones_like(phi_Pm)

        if any(phi_Pm.flatten()>1): import pdb; pdb.set_trace()
        return phi_Pm
