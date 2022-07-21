from packs.utils import constants as ctes
from packs.directories import data_loaded
from .flux_volumes import Flux
from ..IMPEC.properties_calculation import PropertiesCalc
import numpy as np


if ctes.miscible_w:
    from packs.compositional.satability_check_3ph import StabilityCheck
else:
    if data_loaded['compositional_data']['component_data']['constant_K']:
        from packs.compositional.Kflash import StabilityCheck
    else:
        from packs.compositional.stability_check import StabilityCheck

class RiemannSolvers:
    def __init__(self, v0, pretransmissibility):
        self.v0 = v0
        self.pretransmissibility = pretransmissibility

    def UPW(self, M, fprop, Nk_face, P_face, Ft_internal, Fk_face, ponteiro_UPW):
        #"upwind direction"
        Pot_hid = fprop.P + fprop.Pcap #- G[0,:,:]
        Pot_hidj = Pot_hid[0,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[0,ctes.v0[:,1]]
        Fk_internal_faces = np.copy(Fk_face[...,0])
        Fk_internal_faces[:,Pot_hidj>Pot_hidj_up] = np.copy(Fk_face[:,Pot_hidj>Pot_hidj_up,1])

        ponteiro = np.ones(ctes.n_internal_faces,dtype=bool)
        #dNkmax_small = np.max(abs(Nk_face[:,:,0]-Nk_face[:,:,1]),axis=0)<1e-20
        #ponteiro[dNkmax_small] = False
        wave_velocity = np.zeros((ctes.n_components, ctes.n_internal_faces))
        if any(ponteiro): #try to remove this if
            wave_velocity[:,ponteiro],m = self.medium_wave_velocity(M, fprop, Nk_face, P_face, \
            Ft_internal, ponteiro)
        return Fk_internal_faces, wave_velocity

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
        alpha_5 = np.empty((ctes.n_components, len(ftotal[0]), 5))
        if any(ponteiro):
            alpha_5[:,ponteiro,:], eigvec_m = self.wave_velocity_LLF(M, fprop, Nk_face,
                                P_face, ftotal, ponteiro)

        alpha_LR, alpha_m, eigvec_m = self.LRM_wave_velocity(M, fprop, Nk_face,
                                P_face, ftotal, ponteiro)
        alpha_5 = np.concatenate((alpha_LR,alpha_m[...,np.newaxis]),axis=-1)

        #alpha_5[:,~ponteiro] = 0
        #alpha_5[dNk_small] = 0
        #alpha_5, eigvec_m = self.LR_wave_velocity(M, fprop, Nk_face,
        #                        P_face, ftotal, ponteiro)
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

        alphaLR,e = self.LR_wave_velocity(M, fprop, Nk_face, P_face, \
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

    def ROE(self, M, fprop, Nk_face, P_face, ftotal, Fk_face, ponteiro):
        #ponteiro = np.ones_like(ftotal[0],dtype=bool)
        Vpm = fprop.Vp[ctes.v0].sum(axis=-1)/2
        alpha_5, eigvec_m = self.wave_velocity_LLF(M, fprop, Nk_face, P_face, \
            ftotal, ponteiro)
        alpha_m = np.max(abs(alpha_5),axis=-1)
        alpha_LLF = np.max(abs(alpha_5),axis=0)
        alpha_m = self.Harten_entropy_corr(alpha_m, alpha_5[...,[0,1]])
        ponteiro_LLF = self.umbilic_points(alpha_m)
        Fk_internal_faces = np.empty_like(Fk_face[...,0])
        ponteiro_LLF[:] = False

        '''Fk_internal_faces[...,ponteiro_LLF] = self.update_flux_LLF(Fk_face[:,ponteiro_LLF],
            Nk_face[:,ponteiro_LLF], alpha_LLF[ponteiro_LLF,:])'''

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
                L_face, V_face, A_face, xkj_face, Csi_j_face, rho_j_face = \
                    StabilityCheck(P_face, fprop.T).run_init(P_face, z_face, \
                    ksi_W = Csi_j_face[:,ctes.n_phases-1], rho_W = rho_j_face[:,ctes.n_phases-1])
            else:
                L_face = np.ones(len(P_face)); V_face = np.zeros(len(P_face))

                rho_j_face[0,0,:] = np.tile(fprop.rho_j[0,0,self.v0[ponteiro,0]], v)
                rho_j_face[0,1,:] = np.tile(fprop.rho_j[0,1,self.v0[ponteiro,0]], v) #constante, independe da pressao
                #self.reshape_constant_property(fprop.rho_j[0,0:2,:], ponteiro, v)
                Csi_j_face[0,0,:] = np.tile(fprop.Csi_j[0,0,self.v0[ponteiro,0]], v)
                Csi_j_face[0,1,:] = np.tile(fprop.Csi_j[0,1,self.v0[ponteiro,0]], v)
                #self.reshape_constant_property(fprop.Csi_j[0,0:2,:], ponteiro, v)
        else: L_face = []; V_face = []

        if ctes.load_w and not ctes.miscible_w:
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

    def get_Fk_face(self, fprop, M, Nk_face, P_face, Vp_face, ftotal):
        ft_Nks = np.tile(ftotal,2)
        Nks = np.concatenate((Nk_face[:,:,0], Nk_face[:,:,1]), axis=1)
        #Vp_face = np.concatenate((fprop.Vp[ctes.v0[:,0]], fprop.Vp[ctes.v0[:,1]]), axis=0)
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
        return dFkdNk_eigvalue, v

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
