import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from ..stability_check import StabilityCheck
from ..properties_calculation import PropertiesCalc
from .composition_solver import RK3
import scipy.sparse as sp
from sympy import Symbol, diff
import sympy.utilities.lambdify as lambdify
from scipy import integrate
import matplotlib.pyplot as plt
from packs.compositional import prep_FR as ctes_FR
import math

'Todo esse código vai ainda ser ajeitado! As funções de Flux devem ser funçoes para \
calculo independente do fluxo. Funcoes de rotina. E algumas funçes do MUSCL devem também \
seguir essa logica (sao funçoes que nao pertencem ao metodo em si mas sao auxiliares a ele \
e acabam sendo necessarias a outros metodos tambem)'

class Flux:
    """ Class created for computing flux accondingly to the First Order Upwind \
    Method - actually, it's only Flux because of the choice of mobilities and \
    other properties to be taken as function of the potencial gradient. But here \
    flux is computed at the interfaces, and then, volumes, only that."""

    def update_flux(self, fprop, Ft_internal_faces, rho_j_internal_faces, mobilities_internal_faces):
        ''' Main function that calls others '''
        self.Nk = fprop.Nk
        Fj_internal_faces = self.update_Fj_internal_faces(Ft_internal_faces,
            rho_j_internal_faces, mobilities_internal_faces, fprop.Pcap[:,ctes.v0],
            ctes.z[ctes.v0], ctes.pretransmissibility_internal_faces)
        Fk_internal_faces = self.update_Fk_internal_faces(
            fprop.xkj_internal_faces, fprop.Csi_j_internal_faces, Fj_internal_faces)
        fprop.Fk_vols_total = self.update_flux_volumes(Fk_internal_faces)

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
        # M.flux_faces[M.faces.internal] = Ft_internal_faces * M.faces.normal[M.faces.internal].T

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

    def wave_velocity_upw(self, M, fprop, P_old, Ft_internal_faces):
        arg = (fprop.Nk[:,ctes.v0[:,1]] - fprop.Nk[:,ctes.v0[:,0]])==0
        arg_vols = np.sum(arg,axis=0,dtype=bool)
        alpha = np.empty(arg.shape)

        Fj_internal_faces_0 = self.update_Fj_internal_faces(Ft_internal_faces,
            fprop.rho_j[:,:,ctes.v0[:,0]], fprop.mobilities[:,:,ctes.v0[:,0]], fprop.Pcap[:,ctes.v0],
            ctes.z[ctes.v0], ctes.pretransmissibility_internal_faces)
        Fk_internal_faces_0 = self.update_Fk_internal_faces(fprop.xkj[:,:,ctes.v0[:,0]],
            fprop.Csi_j[:,:,ctes.v0[:,0]], Fj_internal_faces_0)

        Fj_internal_faces_1 = self.update_Fj_internal_faces(Ft_internal_faces,
            fprop.rho_j[:,:,ctes.v0[:,1]], fprop.mobilities[:,:,ctes.v0[:,1]], fprop.Pcap[:,ctes.v0],
            ctes.z[ctes.v0], ctes.pretransmissibility_internal_faces)
        Fk_internal_faces_1 = self.update_Fk_internal_faces(fprop.xkj[:,:,ctes.v0[:,1]],
        fprop.Csi_j[:,:,ctes.v0[:,1]], Fj_internal_faces_1)
        alpha[~arg] = ((Fk_internal_faces_1 - Fk_internal_faces_0)[~arg] / \
                    (fprop.Nk[:,ctes.v0[:,1]] - fprop.Nk[:,ctes.v0[:,0]])[~arg])

        alpha[arg] = 1
        Nk_face = fprop.Nk[:,ctes.v0]
        P_face = np.sum(P_old[ctes.v0], axis=1) * 0.5
        ponteiro = np.ones_like(P_face, dtype=bool)
        alpha = fprop.Fk_vols_total
        #ponteiro[arg_vols] = True
        #alpha2 = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces).\
        #         LR_wave_velocity(M, fprop, Nk_face, P_face, Ft_internal_faces, ponteiro)

        #alpha[arg] = np.max(abs(alpha2),axis=2)[arg]
        return alpha

class RiemannSolvers:
    def __init__(self, v0, pretransmissibility):
        self.v0 = v0
        self.pretransmissibility = pretransmissibility

    def LLF(self, M, fprop, Nk_face, P_face, ftotal, Fk_face, ponteiro):
        alpha = self.wave_velocity(M, fprop, Nk_face, P_face, ftotal, np.copy(~ponteiro))
        Fk_internal_faces, alpha_LLF = self.update_flux_LLF(Fk_face[:,~ponteiro,:],
            Nk_face[:,~ponteiro,:], alpha)
        return Fk_internal_faces, alpha_LLF

    def reshape_constant_property(self, y, ponteiro, v):
        try:
            y_aux0 = y[self.v0][ponteiro,0]
            y_aux1 = y[self.v0][ponteiro,1]
        except:
            y_aux0 = y[:,self.v0][:,ponteiro,0]
            y_aux1 = y[:,self.v0][:,ponteiro,1]

        y_reshaped1 = np.tile(y_aux0,(int(v/2) * \
            (np.sign(v - ctes.n_components)**2) + int(v)*(1 - np.sign(v -
            ctes.n_components)**2)))
        y_reshaped2 = np.tile(y_aux1,(int(v/2) * \
            (np.sign(v - ctes.n_components)**2) + int(v/(ctes.n_components+1))*(1 -
            np.sign(v-ctes.n_components))))

        y_reshaped = np.concatenate((y_reshaped1, y_reshaped2), axis=-1)

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

    def get_extrapolated_properties(self, fprop, M, Nk_face, z_face, P_face, Vp, v, ponteiro):
        xkj_face = np.empty((ctes.n_components, ctes.n_phases, len(P_face)))
        Csi_j_face = np.empty((1, ctes.n_phases, len(P_face)))
        rho_j_face = np.empty((1, ctes.n_phases, len(P_face)))

        ''' Flash calculations and properties calculations at each side of the \
        interface '''
        if ctes.compressible_k:
            L_face, V_face, xkj_face[0:ctes.Nc,0,...], xkj_face[0:ctes.Nc,1,...], \
            Csi_j_face[0,0,...], Csi_j_face[0,1,...], rho_j_face[0,0,...], \
            rho_j_face[0,1,...] = StabilityCheck(P_face, fprop.T).run_init(P_face, z_face)

        else:
            L_face = np.ones(len(P_face)); V_face = np.zeros(len(P_face))
            xkj_face[0:ctes.Nc,0:2,:] = 1
            rho_j_face[0,0:2,:] = self.reshape_constant_property(fprop.rho_j[0,0:2,:], ponteiro, v)
            Csi_j_face[0,0:2,:] = self.reshape_constant_property(fprop.Csi_j[0,0:2,:], ponteiro, v)


        if ctes.load_w:
            xkj_face[-1,-1,...] = 1
            xkj_face[-1,0:-1,...] = 0
            xkj_face[0:ctes.Nc,-1,...] = 0

            if data_loaded['compositional_data']['water_data']['mobility']:
                Csi_W0_face = self.reshape_constant_property(fprop.Csi_W0, ponteiro, v)

                Sw_face, Csi_j_face[0,-1,...], rho_j_face[0,-1,...] = \
                PropertiesCalc().update_water_saturation(fprop, Nk_face[-1,...],
                P_face, Vp, Csi_W0_face)

            else:
                Sw_face = self.reshape_constant_property(fprop.Sw, ponteiro, v)
                rho_j_face[0,-1] = self.reshape_constant_property(fprop.rho_j[0,-1], ponteiro, v)
                Csi_j_face[0,-1] = self.reshape_constant_property(fprop.Csi_j[0,-1], ponteiro, v)
        else:
            Sw_face = self.reshape_constant_property(fprop.Sw, ponteiro, v)
            rho_j_face[0,-1] = self.reshape_constant_property(fprop.rho_j[0,-1], ponteiro, v)
            Csi_j_face[0,-1] = self.reshape_constant_property(fprop.Csi_j[0,-1], ponteiro, v)

        So_face, Sg_face =  PropertiesCalc().update_saturations(Sw_face,
            Csi_j_face, L_face, V_face)

        mobilities_face = PropertiesCalc().update_mobilities(fprop, So_face,
        Sg_face, Sw_face, Csi_j_face, xkj_face)

        return mobilities_face, rho_j_face, Csi_j_face, xkj_face

    def Fk_from_Nk(self, fprop, M, Nk, P_face, ftotal, ponteiro):
        ''' Function to compute component flux based on a given composition (Nk) '''

        v = int(len(Nk[0,:])/len(ponteiro[ponteiro]))
        z = Nk[0:ctes.Nc] / np.sum(Nk[0:ctes.Nc], axis = 0)

        Vp_reshaped = self.reshape_constant_property(fprop.Vp, ponteiro, v)

        mobilities, rho_j, Csi_j, xkj = self.get_extrapolated_properties(fprop, M, Nk, z,
        np.tile(P_face[ponteiro],v), Vp_reshaped, v, ponteiro)

        f = Flux()

        Pcap_reshaped = np.concatenate(np.dsplit(np.tile(fprop.Pcap[:,self.v0][:,ponteiro],v),v),axis=1)
        z_reshaped = np.concatenate(np.hsplit(np.tile(ctes.z[self.v0][ponteiro],v),v),axis=0)
        Fj = f.update_Fj_internal_faces(ftotal, rho_j, mobilities, Pcap_reshaped,
        z_reshaped, np.tile(self.pretransmissibility[ponteiro],v))

        Fk = f.update_Fk_internal_faces(xkj, Csi_j, Fj)
        return Fk

    def get_Fk_face(self, fprop, M, Nk_face, P_face, ftotal):
        ''' Function that computes the flux in each face side (Left and Right)'''

        Fk_face = self.Fk_from_Nk(fprop, M, np.concatenate((Nk_face[:,:,0],Nk_face[:,:,1]),axis=1), P_face,
            np.tile(ftotal[:,ponteiro],2), np.ones(len(ftotal[0,:]), dtype=bool))

        Fk_face = np.concatenate(np.hsplit(Fk_face[:,:,np.newaxis],2),axis = 2)
        return Fk_face

    def LR_wave_velocity(self, M, fprop, Nk_face, P_face, ftotal, ponteiro):
        delta = 0.001

        Nk_face_plus = Nk_face[np.newaxis,:,ponteiro,:] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        Nk_face_minus = Nk_face[np.newaxis,:,ponteiro,:] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        matrix_deltas = np.identity(ctes.n_components)[:,:,np.newaxis, np.newaxis] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]),2])
        Nk_face_plus += delta * 0.5 * matrix_deltas
        Nk_face_minus -= delta * 0.5 * matrix_deltas
        Nk_face_plus = np.concatenate(np.split(Nk_face_plus, ctes.n_components),axis=2)[0,...]
        Nk_face_minus = np.concatenate(np.split(Nk_face_minus, ctes.n_components),axis=2)[0,...]
        Nk_face_plus = np.concatenate(np.dsplit(Nk_face_plus, 2),axis=1)[:,:,0]
        Nk_face_minus = np.concatenate(np.dsplit(Nk_face_minus, 2),axis=1)[:,:,0]
        dFkdNk = ((self.Fk_from_Nk(fprop, M, Nk_face_plus, P_face, np.tile(ftotal[:,ponteiro],ctes.n_components*2), ponteiro) -
        self.Fk_from_Nk(fprop, M, Nk_face_minus, P_face, np.tile(ftotal[:,ponteiro],ctes.n_components*2), ponteiro))/(Nk_face_plus - Nk_face_minus).sum(axis=0))
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,np.newaxis],2),axis=2)
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,:,np.newaxis],ctes.n_components),axis=3)
        dFkdNk = dFkdNk.transpose(2,1,0,3)

        eigval1, v = np.linalg.eig(dFkdNk)
        dFkdNk_eigvalue = eigval1.T
        return dFkdNk_eigvalue


    def wave_velocity(self, M, fprop, Nk_face, P_face, ftotal, ponteiro):
        delta = 0.001

        Nkm = (Nk_face[:,ponteiro,1] + Nk_face[:,ponteiro,0])/2

        Nkg = Nkm[:,:,np.newaxis] + (Nk_face[:,ponteiro] - Nkm[:,:,np.newaxis])/(3**(1/2))

        Nkg_plus = Nkg[np.newaxis,...] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        Nkg_minus = Nkg[np.newaxis,...] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        matrix_deltas = np.identity(ctes.n_components)[:,:,np.newaxis, np.newaxis] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]),2])
        Nkg_plus += delta * 0.5 * matrix_deltas
        Nkg_minus -= delta * 0.5 * matrix_deltas
        Nkg_plus = np.concatenate(np.split(Nkg_plus, ctes.n_components),axis=2)[0,...]
        Nkg_minus = np.concatenate(np.split(Nkg_minus, ctes.n_components),axis=2)[0,...]
        Nkg_plus = np.concatenate(np.dsplit(Nkg_plus, 2),axis=1)[:,:,0]
        Nkg_minus = np.concatenate(np.dsplit(Nkg_minus, 2),axis=1)[:,:,0]
        dFkdNk_gauss = ((self.Fk_from_Nk(fprop, M, Nkg_plus, P_face, np.tile(ftotal[:,ponteiro],ctes.n_components*2),ponteiro) -
        self.Fk_from_Nk(fprop, M, Nkg_minus, P_face, np.tile(ftotal[:,ponteiro],ctes.n_components*2), ponteiro))/(Nkg_plus - Nkg_minus).sum(axis=0))
        dFkdNk_gauss = np.concatenate(np.hsplit(dFkdNk_gauss[:,:,np.newaxis],2),axis=2)
        dFkdNk_gauss = np.concatenate(np.hsplit(dFkdNk_gauss[:,:,:,np.newaxis],ctes.n_components),axis=3)
        dFkdNk_gauss = dFkdNk_gauss.transpose(2,1,0,3)

        Nkm_plus = Nkm[np.newaxis,:,:] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro])])
        Nkm_minus = Nkm[np.newaxis,:,:] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro])])
        matrix_deltas = np.identity(ctes.n_components)[:,:,np.newaxis] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro])])
        Nkm_plus  += delta * 0.5 * matrix_deltas
        Nkm_minus -= delta * 0.5 * matrix_deltas
        Nkm_plus = np.concatenate(np.split(Nkm_plus, ctes.n_components),axis=2)[0,...]
        Nkm_minus = np.concatenate(np.split(Nkm_minus, ctes.n_components),axis=2)[0,...]
        dFkdNk_m = ((self.Fk_from_Nk(fprop, M, Nkm_plus, P_face, np.tile(ftotal[:,ponteiro],ctes.n_components), ponteiro) -
        self.Fk_from_Nk(fprop, M, Nkm_minus, P_face, np.tile(ftotal[:,ponteiro],ctes.n_components), ponteiro))/ (Nkm_plus - Nkm_minus).sum(axis=0))
        dFkdNk_m = np.concatenate(np.hsplit(dFkdNk_m[:,:,np.newaxis],ctes.n_components),axis = 2)
        dFkdNk_m = dFkdNk_m.transpose(1,0,2)

        Nk_face_plus = Nk_face[np.newaxis,:,ponteiro,:] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        Nk_face_minus = Nk_face[np.newaxis,:,ponteiro,:] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        matrix_deltas = np.identity(ctes.n_components)[:,:,np.newaxis, np.newaxis] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]),2])
        Nk_face_plus += delta * 0.5 * matrix_deltas
        Nk_face_minus -= delta * 0.5 * matrix_deltas
        Nk_face_plus = np.concatenate(np.split(Nk_face_plus, ctes.n_components),axis=2)[0,...]
        Nk_face_minus = np.concatenate(np.split(Nk_face_minus, ctes.n_components),axis=2)[0,...]
        Nk_face_plus = np.concatenate(np.dsplit(Nk_face_plus, 2),axis=1)[:,:,0]
        Nk_face_minus = np.concatenate(np.dsplit(Nk_face_minus, 2),axis=1)[:,:,0]
        dFkdNk = ((self.Fk_from_Nk(fprop, M, Nk_face_plus, P_face, np.tile(ftotal[:,ponteiro],ctes.n_components*2), ponteiro) -
        self.Fk_from_Nk(fprop, M, Nk_face_minus, P_face, np.tile(ftotal[:,ponteiro],ctes.n_components*2), ponteiro))/(Nk_face_plus - Nk_face_minus).sum(axis=0))
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,np.newaxis],2),axis=2)
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,:,np.newaxis],ctes.n_components),axis=3)
        dFkdNk = dFkdNk.transpose(2,1,0,3)

        eigval1, v = np.linalg.eig(dFkdNk)
        dFkdNk_eigvalue = eigval1.T
        eigval2, v = np.linalg.eig(dFkdNk_gauss)
        dFkdNk_gauss_eigvalue = eigval2.transpose(2,1,0)
        eigval3, v = np.linalg.eig(dFkdNk_m)
        dFkdNk_m_eigvalue = eigval3.T

        alpha = np.concatenate((dFkdNk_eigvalue[:,ponteiro], dFkdNk_gauss_eigvalue), axis=-1)
        alpha = np.concatenate((alpha, dFkdNk_m_eigvalue[:,:,np.newaxis]), axis=-1)
        return alpha

    def update_flux_LLF(self, Fk_face_LLF_all, Nk_face_LLF, alpha_LLF):
        #alpha2 = np.concatenate((alpha_LLF, alpha_RH[:,:,np.newaxis]),axis=-1)
        alpha_RH = (Fk_face_LLF_all[:,:,1] - Fk_face_LLF_all[:,:,0]) / \
            (Nk_face_LLF[:,:,1] - Nk_face_LLF[:,:,0])
        alpha_LLF[abs(Nk_face_LLF[:,:,1] - Nk_face_LLF[:,:,0])>1e-3] = alpha_RH[abs(Nk_face_LLF[:,:,1] -
            Nk_face_LLF[:,:,0])>1e-3][:,np.newaxis]
        alpha = np.max(abs(alpha_LLF),axis = 0)

        Fk_face_LLF = 0.5*(Fk_face_LLF_all.sum(axis=-1) - np.max(abs(alpha),axis=-1) * \
                    (Nk_face_LLF[:,:,1] - Nk_face_LLF[:,:,0]))
        return Fk_face_LLF, alpha

class MUSCL:

    """ Class created for the second order MUSCL implementation for the \
    calculation of the advective terms """

    def run(self, M, fprop, wells, P_old, ftot, t, delta_t):
        self.t = t+delta_t
        ''' Global function that calls others '''
        self.P_face = np.sum(P_old[ctes.v0], axis=1) * 0.5
        dNk_vols = self.volume_gradient_reconstruction(M, fprop, wells)
        dNk_face, dNk_face_neig = self.get_faces_gradient(M, fprop, dNk_vols)
        order = 2
        if order == 2:
            Phi = self.Van_Leer_slope_limiter(dNk_face, dNk_face_neig)

            Nk_face, z_face = self.get_extrapolated_compositions(fprop, Phi, dNk_face_neig)
        else:
            Phi = np.zeros_like(dNk_face_neig)
            Nk_face, z_face = self.get_extrapolated_compositions(fprop, Phi, dNk_face_neig)
        #G = self.update_gravity_term() # for now, it has no gravity
        alpha = self.update_flux(M, fprop, Nk_face, ftot)
        return alpha

    def volume_gradient_reconstruction(self, M, fprop, wells):
        neig_vols = M.volumes.bridge_adjacencies(M.volumes.all,2,3)
        matriz = np.zeros((ctes.n_volumes,ctes.n_volumes))

        lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        data = np.array([np.ones(len(ctes.v0[:, 0])), np.ones(len(ctes.v0[:, 0])),
                        np.zeros(len(ctes.v0[:, 0])), np.zeros(len(ctes.v0[:, 0]))]).flatten()
        all_neig = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes)).toarray()
        all_neig = all_neig.astype(int)
        all_neig2 = all_neig + np.identity(ctes.n_volumes)
        allneig2 = all_neig2.astype(int)

        Nk_neig =  fprop.Nk[:,np.newaxis,:] * all_neig2[np.newaxis,:,:]
        Nk = Nk_neig.transpose(0,2,1)
        pos_neig = M.data['centroid_volumes'].T[:,np.newaxis,:] * all_neig2[np.newaxis,:,:]

        pos = pos_neig.transpose(0,2,1)

        ds = pos_neig - pos
        ds_norm = np.linalg.norm(ds, axis=0)
        versor_ds = np.empty(ds.shape)
        versor_ds[:,ds_norm==0] = 0
        versor_ds[:,ds_norm!=0] = ds[:,ds_norm!=0] / ds_norm[ds_norm!=0]
        dNk = Nk_neig - Nk

        dNk_by_axes = np.repeat(dNk[:,np.newaxis,:,:],3,axis=1)
        dNk_by_axes = dNk_by_axes * versor_ds
        dNk_vols = dNk_by_axes.sum(axis = 3)

        ds_vols = ds * versor_ds
        ds_vols = ds_vols.sum(axis = 2)
        dNkds_vols = np.copy(dNk_vols)
        dNkds_vols[:,ds_vols!=0] = dNk_vols[:,ds_vols != 0] / ds_vols[ds_vols != 0][np.newaxis,:]
        self.all_neig = all_neig.sum(axis=1)
        self.faces_contour = self.identify_contour_faces()
        #dNkds_vols[:,:,all_neig.sum(axis=1)==1] = 0 #*dNkds_vols[:,:,all_neig.sum(axis=1)==1]

        dNkds_vols[:,:,ctes.v0[self.faces_contour].flatten()] = 0 # zero in the contour volumes
        return dNkds_vols

    def get_faces_gradient(self, M, fprop, dNkds_vols):
        dNk_face =  fprop.Nk[:,ctes.v0[:,1]] - fprop.Nk[:,ctes.v0[:,0]]
        ds_face = M.data['centroid_volumes'][ctes.v0[:,1],:] -  M.data['centroid_volumes'][ctes.v0[:,0],:]
        dNk_face_vols = 2. * (dNkds_vols[:,:,ctes.v0] * ds_face.T[np.newaxis,:,:,np.newaxis]).sum(axis=1)
        dNk_face_neig = dNk_face_vols - dNk_face[:,:,np.newaxis]
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
        import pdb; pdb.set_trace()
        return Phi

    def get_extrapolated_compositions(self, fprop, Phi, dNk_face_neig):
        Nk_face = fprop.Nk[:,ctes.v0] + Phi / 2 * dNk_face_neig
        z_face = Nk_face[0:ctes.Nc] / np.sum(Nk_face[0:ctes.Nc], axis = 0)
        return Nk_face, z_face

    def MLP_slope_limiter(self, M, fprop, dNk_face_neig, dNk_face):
        np.seterr(divide='ignore', invalid='ignore')
        '------------------------Seccond Order term----------------------------'
        delta_inf = dNk_face_neig/2
        delta_inf[:,:,1] = - delta_inf[:,:,1]
        Phi_MLP = np.ones_like(dNk_face_neig)
        r = dNk_face[:,:,np.newaxis]/(dNk_face_neig)
        r[r<0] = 0
        Phi_MLP[abs(dNk_face_neig)>=1e-10] = r[abs(dNk_face_neig)>=1e-10]
        Phi_MLP[r > 1] = 1
        e2 = 0#(5 * (M.data['area'][M.faces.internal])**3)[np.newaxis,:,np.newaxis]
        #Phi_MLP = 1 / delta_inf * (((delta_sup**2 + e2) * delta_inf + 2 * delta_inf**2 * delta_sup)
        #    / (delta_sup**2 + 2 * delta_inf**2 + delta_sup * delta_inf + e2))
        #Phi_MLP[delta_inf==0] = 0
        #Nk_face = fprop.Nk[:,ctes.v0] + Phi_MLP*delta_inf
        #Phi_MLP[Nk_face < Nk_aux_min] = 0
        #Phi_MLP[Nk_face > Nk_aux_max] = 0

        '-------------------------Third Order term-----------------------------'
        Theta_MLP2 = np.ones_like(Phi_MLP)
        d2Nk_vols = dNk_face[:,:,np.newaxis] - dNk_face_neig
        aux_d2Nk = np.copy(d2Nk_vols)
        d2Nk_vols[:,:,1] = - d2Nk_vols[:,:,1]

        Nk_face = fprop.Nk[:,ctes.v0] + delta_inf #+ 1/8 * aux_d2Nk
        cond1 = (Nk_face > np.max(fprop.Nk[:,ctes.v0],axis=2)[:,:,np.newaxis])
        cond2 = (Nk_face < np.min(fprop.Nk[:,ctes.v0],axis=2)[:,:,np.newaxis])
        Theta_MLP2[cond1 + cond2] = 0

        Nk_face = fprop.Nk[:,ctes.v0] + delta_inf + 1/8 * aux_d2Nk
        cond3 = (Nk_face > np.min(fprop.Nk[:,ctes.v0],axis=2)[:,:,np.newaxis]) * \
            (dNk_face_neig > 0) * (d2Nk_vols < 0)
        cond4 = (Nk_face < np.max(fprop.Nk[:,ctes.v0],axis=2)[:,:,np.newaxis]) * \
            (dNk_face_neig < 0) * (d2Nk_vols > 0)

        Theta_MLP2_aux = Theta_MLP2[Theta_MLP2==0]
        Theta_MLP2_aux[(cond3 + cond4)[Theta_MLP2==0]] = 1
        #cond5 = (np.abs((Nk_face - fprop.Nk[:,ctes.v0])/fprop.Nk[:,ctes.v0]) <= 0.001)
        #Theta_MLP2_aux[cond5[Theta_MLP2==0]] = 1
        Theta_MLP2[Theta_MLP2==0] = Theta_MLP2_aux

        #Theta_MLP2[Phi_MLP==0] = 0
        #Nk_face = fprop.Nk[:,ctes.v0] + delta_inf + 1/8 * aux_d2Nk
        #Theta_MLP2[Nk_face > np.max(fprop.Nk[:,ctes.v0],axis=2)[:,:,np.newaxis]] = 0
        #Theta_MLP2[Nk_face < np.min(fprop.Nk[:,ctes.v0],axis=2)[:,:,np.newaxis]] = 0
        #import pdb; pdb.set_trace()

        'Creating mapping between volumes and internal faces, to get minimum \
        value of the limiter projection at the face. It will be a value for \
        each control volume, and not for each face, as it was beeing done '
        ctes_FR.vols_vec = -np.ones((ctes.n_volumes,2),dtype=int)
        lines = np.arange(ctes.n_internal_faces)
        ctes_FR.vols_vec[ctes.v0[:,0],0] = lines
        ctes_FR.vols_vec[ctes.v0[:,1],1] = lines
        contour_vols = np.argwhere(ctes_FR.vols_vec<0)[:,0]
        ctes_FR.vols_vec[ctes_FR.vols_vec < 0] = ctes_FR.vols_vec[contour_vols,:][ctes_FR.vols_vec[contour_vols]>=0]

        Phi_MLPv = np.empty((ctes.n_components,ctes.n_volumes,2))
        Phi_MLPv[:,:,0] = Phi_MLP[:,ctes_FR.vols_vec[:,0],0]
        Phi_MLPv[:,:,1] = Phi_MLP[:,ctes_FR.vols_vec[:,1],1]
        Phi_MLPv = np.min(Phi_MLPv,axis=2)
        Phi_MLP = Phi_MLPv[:,ctes.v0]

        Theta_MLP2v = np.empty((ctes.n_components,ctes.n_volumes,2))
        Theta_MLP2v[:,:,0] = Theta_MLP2[:,ctes_FR.vols_vec[:,0],0]
        Theta_MLP2v[:,:,1] = Theta_MLP2[:,ctes_FR.vols_vec[:,1],1]
        Theta_MLP2v = np.min(Theta_MLP2v,axis=2)
        Theta_MLP2 = Theta_MLP2v[:,ctes.v0]
        Phi_MLP[:,:,1] = - Phi_MLP[:,:,1]
        Theta_MLP2[:,:,1] = - Theta_MLP2[:,:,1]
        return Phi_MLP, Theta_MLP2

    def get_extrapolated_compositions_MLP(self, fprop, Phi_MLP, Theta_MLP2, dNk_face_neig, dNk_face):
        'Activating limiters at the faces of the contour control volumes '
        #r_face = dNk_face[:,:,np.newaxis] / dNk_face_neig
        d2Nk_vols = dNk_face[:,:,np.newaxis] - dNk_face_neig
        d2Nk_vols[:,:,1] = -d2Nk_vols[:,:,1]
        Nk_face = fprop.Nk[:,ctes.v0] + Phi_MLP/2 * dNk_face_neig + Theta_MLP2/8 * d2Nk_vols
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

    def flux_calculation_conditions(self, alpha):
        ponteiro_LLF = np.ones((ctes.n_components,ctes.n_internal_faces),dtype=bool)
        ponteiro_LLF[alpha[:,:,0] * alpha[:,:,1] <= 0] = False
        difs = np.empty((ctes.n_components, ctes.n_components,ctes.n_internal_faces, 2))
        ind = np.arange(ctes.n_components).astype(int)
        for k in range(ctes.n_components):
            difs[k] = abs(alpha[k,:] - alpha[ind,:])
            difs[k,k] = 1e5
        cond = np.min(difs,axis = 1)
        ponteiro_LLF = ~ponteiro_LLF.sum(axis=0,dtype=bool)
        ponteiro_LLF2 = np.ones((ctes.n_components,ctes.n_internal_faces),dtype=bool)
        ponteiro_LLF2[cond[:,:,0]<0.01*np.max(abs(alpha[:,:,0]),axis=0)] = False
        ponteiro_LLF2[cond[:,:,1]<0.01*np.max(abs(alpha[:,:,1]),axis=0)] = False
        ponteiro_LLF2 = ponteiro_LLF2.sum(axis=0,dtype=bool)
        ponteiro_LLF += ~ponteiro_LLF2
        ponteiro_LLF = ~ponteiro_LLF
        return ponteiro_LLF

    def identify_contour_faces(self):
        vols_contour = np.argwhere(self.all_neig==1).flatten()
        faces_contour = np.empty_like(vols_contour)

        for i in range(len(vols_contour)):
            try: faces_contour[i] = np.argwhere(ctes.v0[:,0] == vols_contour[i]).flatten()
            except: faces_contour[i] = np.argwhere(ctes.v0[:,1] == vols_contour[i]).flatten()
        return faces_contour

    def update_flux(self, M, fprop, Nk_face, ftotal):
        Fk_internal_faces = np.empty((ctes.n_components,ctes.n_internal_faces))
        alpha_wv = np.empty((ctes.n_internal_faces, 5))
        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)
        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, self.P_face, ftotal)

        ponteiro = np.zeros(ctes.n_internal_faces,dtype=bool)

        #LR_eigval, d2FkdNk_eigval = self.get_LR_eigenvalues_Serna(M, fprop, Nk_face)
        #LR_eigval = self.get_LR_eigenvalues(M, fprop, Nk_face,  np.ones(ctes.n_internal_faces,dtype=bool))
        #ponteiro = self.flux_calculation_conditions_Serna(LR_eigval, d2FkdNk_eigval)
        #ponteiro = self.flux_calculation_conditions(LR_eigval)

        Fk_internal_faces[:,~ponteiro], alpha_wv[~ponteiro,:] = RS.LLF(M, fprop, Nk_face, self.P_face,
            ftotal, Fk_face, ponteiro)


        '''if any(ponteiro):
            alpha_wv[:,ponteiro,:] = 0
            #alpha_wv[:,ponteiro,0:2] = LR_eigval[:,ponteiro]
            Fk_internal_faces[:,ponteiro] = self.update_flux_upwind(fprop,
                                                Fk_face[:,ponteiro,:], np.copy(ponteiro))'''

        '-------- Perform volume balance to obtain flux through volumes -------'
        fprop.Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)
        return alpha_wv

    '''def get_LR_eigenvalues_Serna(self, M, fprop, Nk_face):
        dFkdNk = np.empty((ctes.n_internal_faces, ctes.n_components, ctes.n_components))
        dFkdNk_eigvalue = np.empty((ctes.n_components,ctes.n_internal_faces, 2))
        d2FkdNk = np.empty_like(dFkdNk)
        d2FkdNk_eigvalue = np.empty((ctes.n_components,ctes.n_internal_faces, 2))
        delta = 0.001
        for i in range(2):
            for k in range(0,ctes.n_components):
                Nk_face_plus = np.copy(Nk_face[:,:,i])
                Nk_face_minus = np.copy(Nk_face[:,:,i])
                Nk_face_plus[k] += delta*0.5
                Nk_face_minus[k] -= delta*0.5
                dFkdNk[:,:,k] = ((self.Fk_from_Nk(fprop, M, Nk_face_plus, np.ones(ctes.n_internal_faces, dtype=bool)) -
                                    self.Fk_from_Nk(fprop, M, Nk_face_minus, np.ones(ctes.n_internal_faces, dtype=bool)))\
                                    /(Nk_face_plus[k]-Nk_face_minus[k])).T
                d2FkdNk[:,:,k] = ((self.dFk_dNk(fprop, M, Nk_face_plus, delta, k, np.ones(ctes.n_internal_faces, dtype=bool)) -
                                    self.dFk_dNk(fprop, M, Nk_face_minus, delta, k, np.ones(ctes.n_internal_faces, dtype=bool)))\
                                    /(Nk_face_plus[k]-Nk_face_minus[k])).T
            eigval1, v = np.linalg.eig(dFkdNk)
            dFkdNk_eigvalue[:,:,i] = eigval1.T
            eigval2, v = np.linalg.eig(d2FkdNk)
            d2FkdNk_eigvalue[:,:,i] = eigval2.T
        return dFkdNk_eigvalue, d2FkdNk_eigvalue'''

    def get_LR_eigenvalues(self, M, fprop, Nk_face, ponteiro):

        delta = 0.001
        Nk_face_plus = Nk_face[np.newaxis,:,ponteiro,:] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        Nk_face_minus = Nk_face[np.newaxis,:,ponteiro,:] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]), 2])
        matrix_deltas = np.identity(ctes.n_components)[:,:,np.newaxis, np.newaxis] * np.ones([ctes.n_components, ctes.n_components, len(ponteiro[ponteiro]),2])
        Nk_face_plus += delta * 0.5 * matrix_deltas
        Nk_face_minus -= delta * 0.5 * matrix_deltas
        Nk_face_plus = np.concatenate(np.split(Nk_face_plus, ctes.n_components),axis=2)[0,...]
        Nk_face_minus = np.concatenate(np.split(Nk_face_minus, ctes.n_components),axis=2)[0,...]
        Nk_face_plus = np.concatenate(np.dsplit(Nk_face_plus, 2),axis=1)[:,:,0]
        Nk_face_minus = np.concatenate(np.dsplit(Nk_face_minus, 2),axis=1)[:,:,0]
        dFkdNk = ((self.Fk_from_Nk(fprop, M, Nk_face_plus, ponteiro) -
        self.Fk_from_Nk(fprop, M, Nk_face_minus, ponteiro))/(Nk_face_plus - Nk_face_minus).sum(axis=0))
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,np.newaxis],2),axis=2)
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,:,np.newaxis],ctes.n_components),axis=3)
        dFkdNk = dFkdNk.transpose(2,1,0,3)

        eigval1, v = np.linalg.eig(dFkdNk)
        dFkdNk_eigvalue = eigval1.T

        return dFkdNk_eigvalue


    def update_flux_upwind(self, fprop, Fk_face_upwind_all, ponteiro):
        Fk_face_upwind = np.empty_like(Fk_face_upwind_all[:,:,0])

        Pot_hid = fprop.P #+ fprop.Pcap
        Pot_hidj = Pot_hid[ctes.v0[:,0]][ponteiro] #- G[0,:,:,0]
        Pot_hidj_up = Pot_hid[ctes.v0[:,1]][ponteiro] #- G[0,:,:,1]

        Fk_face_upwind[Fk_face_upwind_all[:,:,1] <= Fk_face_upwind_all[:,:,0]] = \
            Fk_face_upwind_all[Fk_face_upwind_all[:,:,1] <= Fk_face_upwind_all[:,:,0], 0]
        Fk_face_upwind[Fk_face_upwind_all[:,:,1] > Fk_face_upwind_all[:,:,0]] = \
            Fk_face_upwind_all[Fk_face_upwind_all[:,:,1] > Fk_face_upwind_all[:,:,0], 1]
        #import pdb; pdb.set_trace()
        #Fk_face_upwind = np.max(Fk_face_upwind_all,axis=-1)
        return Fk_face_upwind

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
        self.t = t
        Nk_SP = np.copy(Nk_SP_old)

        q_SP = q[:,:,np.newaxis] * np.ones_like(Nk_SP)


        dFk_SP, wave_velocity = self.dFk_SP_from_Pspace(M, fprop, wells, Ft_internal_faces, np.copy(Nk_SP), q_SP, P_old)
        Nk_SP = RK3.update_composition_RK3_1(np.copy(Nk_SP_old), q_SP, dFk_SP, delta_t)
        Nk_SP = self.MLP_slope_limiter(M, fprop, Nk_SP)

        dFk_SP, wave_velocity = self.dFk_SP_from_Pspace(M, fprop, wells, Ft_internal_faces, np.copy(Nk_SP), q_SP, P_old)
        Nk_SP = RK3.update_composition_RK3_2(np.copy(Nk_SP_old), q_SP, np.copy(Nk_SP), dFk_SP, delta_t)
        Nk_SP = self.MLP_slope_limiter(M, fprop, Nk_SP)

        dFk_SP, wave_velocity = self.dFk_SP_from_Pspace(M, fprop, wells, Ft_internal_faces, Nk_SP, q_SP, P_old)
        Nk_SP = RK3.update_composition_RK3_3(np.copy(Nk_SP_old), q_SP, np.copy(Nk_SP), dFk_SP, delta_t)
        Nk_SP = self.MLP_slope_limiter(M, fprop, Nk_SP)
        fprop.Fk_vols_total = np.min(abs(dFk_SP),axis=2)

        Nk = 1 / sum(ctes_FR.weights) * np.sum(ctes_FR.weights * Nk_SP,axis=2)

        z = Nk[0:ctes.Nc,:] / np.sum(Nk[0:ctes.Nc,:], axis = 0)

        if len(Nk_SP[Nk_SP<0])>0:
            import pdb; pdb.set_trace()
        return wave_velocity, Nk, z, Nk_SP

    def dFk_SP_from_Pspace(self, M, fprop, wells, Ft_internal_faces, Nk_SP, q_SP, P_old):

        Ft_SP = self.total_flux_SP(fprop, wells, M, Ft_internal_faces)
        Fk_SP = self.component_flux_SP(fprop, M, Nk_SP, P_old, Ft_SP)

        Fk_faces, Fk_vols_RS_neig, wave_velocity = self.Riemann_Solver(M, fprop, Nk_SP, Fk_SP, P_old, Ft_internal_faces)
        Fk_D = np.sum(Fk_SP[:,:,:,np.newaxis] * ctes_FR.L[np.newaxis,np.newaxis,:], axis=2)
        dFk_D = np.sum(Fk_SP[:,:,:,np.newaxis] * ctes_FR.dL[np.newaxis,np.newaxis,:], axis=2)
        dFk_C = self.dFlux_Continuous(Fk_D, Fk_vols_RS_neig)
        dFk_Pspace = (dFk_C + dFk_D)

        #Fk_vols_RS_neig[:,wells['all_wells'],0] = -Fk_vols_RS_neig[:,wells['all_wells'],0]
        if not data_loaded['compositional_data']['water_data']['mobility']:
            dFk_Pspace[-1,:] = 0

        #dFk_Pspace[:,wells['all_wells'],1:] = 0
        #dFk_Pspace[:,wells['all_wells'],0] = (Fk_vols_RS_neig[:,wells['all_wells']]).sum(axis=2)/2
        dFk_SP = dFk_Pspace @ ctes_FR.x_points
        dFk_SP = - 2 * dFk_SP #this way only works for uniform mesh

        #up! transforming from local space to original global space (this could be done to the g and L functions
        #only, however, I rather do like this, so it's done just once)

        #dFk_SP [:,wells['all_wells'],[0,-1]] = -Fk_vols_RS_neig[:,wells['all_wells'], [0,1]]#.sum(axis=2)

        #dFk_SP = np.empty_like(Fk_SP[:,ctes.vols_no_wells])
        #for i in range(ctes_FR.n_points):
        #    dFk_SP[:,:,i] = np.array(dFk_func(ctes_FR.points[i]))
        return dFk_SP, wave_velocity

    def total_flux_SP(self, fprop, wells, M, Ft_internal_faces):
        'RTo'
        phi = np.empty((len(ctes_FR.points),2))
        phi[:,0] = 1 / 4 * (1 + ctes_FR.points)
        phi[:,1] = 1 / 4 * (1 - ctes_FR.points)

        Ft_face_phi = (Ft_internal_faces[:,:,np.newaxis,np.newaxis] * phi[np.newaxis,np.newaxis,:])
        #Fk_face_phi_reshaped = np.concatenate(np.split(Fk_face_phi, ctes_FR.n_points,axis=3),axis=1)[:,:,:,0]

        'Look for a faster way to do that'
        Ft_SP_reshaped = np.empty((1,ctes.n_volumes,ctes_FR.n_points))
        contours = np.array([0,ctes_FR.n_points-1])
        for i in range(ctes_FR.n_points):
            lines = np.array([np.zeros_like(ctes_FR.v0[:,0]), np.zeros_like(ctes_FR.v0[:,1])]).astype(int).flatten()
            cols = np.array([ctes_FR.v0[:,0], ctes_FR.v0[:,1]]).flatten()
            data = np.array([Ft_face_phi[:,:,i,0], Ft_face_phi[:,:,i,1]]).flatten()
            Ft_SP_reshaped[:,:,i] = sp.csc_matrix((data, (lines, cols)), shape = (1, ctes.n_volumes)).toarray()
        Ft_SP = 2 * np.concatenate(np.dsplit(Ft_SP_reshaped, ctes_FR.n_points), axis = 2)
        #Fk_SP[:,wells['all_wells']] = Fk_face[:,wells['all_wells'],1][:,:,np.newaxis]/2
        return Ft_SP

    def component_flux_SP(self, fprop, M, Nk_SP, P, Ft_SP):
        P_SP = P[:,np.newaxis] * np.ones_like(Ft_SP[0])
        ponteiro = np.ones(ctes.n_volumes,dtype=bool)
        v0 = np.arange(ctes.n_volumes)[:,np.newaxis] * np.ones((ctes.n_volumes,ctes_FR.n_points))
        Nk_SP_flatt = np.concatenate(np.dsplit(Nk_SP, ctes_FR.n_points),axis=1)[:,:,0]
        Ft_SP_flatt = np.concatenate(np.dsplit(Ft_SP, ctes_FR.n_points),axis=1)[:,:,0]
        pretr = ctes.pretransmissibility_internal_faces[ctes_FR.vols_vec][:,0]

        Fk_SP = RiemannSolvers(v0.astype(int), pretr).Fk_from_Nk(fprop, M, Nk_SP_flatt, P, Ft_SP_flatt, ponteiro)
        Fk_SP = np.concatenate(np.hsplit(Fk_SP[:,:,np.newaxis], ctes_FR.n_points), axis = 2)
        return Fk_SP

    def Riemann_Solver(self, M, fprop, Nk_SP, Fk_SP, P_old, Ft_internal_faces):

        Nk_faces = np.empty((ctes.n_components, ctes.n_internal_faces, 2))

        Nk_faces[:,:,1] = Nk_SP[:,ctes_FR.v0[:,1],0] #Nk faces a esquerda dos volumes
        Nk_faces[:,:,0] = Nk_SP[:,ctes_FR.v0[:,0],-1] #Nk nas faces a direita
        P_face = np.sum(P_old[ctes_FR.v0],axis=1) * 0.5

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
            Ft_internal_faces, Fk_faces, np.zeros(ctes.n_internal_faces,dtype=bool))

        'Obtaining Flux at each CV side - by finding faces that compounds the CV \
        this only works for 1D problems'
        vols_vec = -np.ones((ctes.n_volumes,2),dtype=int)
        lines = np.arange(ctes.n_internal_faces)
        vols_vec[ctes_FR.v0[:,0],1] = lines
        vols_vec[ctes_FR.v0[:,1],0] = lines

        Fk_vols_RS_neig = Fk_face_RS[:,ctes_FR.vols_vec]

        Fk_vols_RS_neig[:,vols_vec<0] = 0 #Fk_face_contour_RS
        return Fk_faces, Fk_vols_RS_neig, alpha_wv

    def dFlux_Continuous(self, Fk_D, Fk_vols_RS_neig):
        #norm_point_coords =  np.array([(x,y,z) for x in coords for y in coords for z in coords])
        x_left = np.array([(-1)**i for i in range(ctes_FR.n_points)])
        x_right = np.ones(ctes_FR.n_points)
        Fk_D_l = Fk_D @ x_left
        Fk_D_r = Fk_D @ x_right
        dFk_C = (Fk_vols_RS_neig[:,:,0] - Fk_D_l)[:,:,np.newaxis] * ctes_FR.dgLB[np.newaxis,:] + \
                (Fk_vols_RS_neig[:,:,1] - Fk_D_r)[:,:,np.newaxis] * ctes_FR.dgRB[np.newaxis,:]

        return dFk_C

    def MLP_slope_limiter(self, M, fprop, Nk_SP_in):

        inds = np.array([0,-1])

        Nk = (np.linalg.inv(ctes_FR.V)[np.newaxis,] @ Nk_SP_in[:,:,:, np.newaxis])[:,:,:,0]
        Nk_SP = self.projections(Nk, ctes_FR.n_points-1)

        #machine_error = np.max(abs(Nk_SP - Nk_SP_in))
        #if machine_error<np.finfo(np.float64).eps: machine_error = np.finfo(np.float64).eps
        machine_error = 1e-8

        inds = np.array([0,-1])

        'Projected n=0'
        Nk_P0 = self.projections(Nk, 0)
        Nk_P0_vertex = Nk_P0[:,:,inds]

        'Projected n=1'
        Nk_P1 = self.projections(Nk, 1)
        Nk_P1_vertex = Nk_P1[:,:,inds]

        Phi_P1 = self.P1_limiter(M, Nk, Nk_SP, Nk_P0_vertex, Nk_P1_vertex, machine_error)

        phi_P2 = np.zeros_like(Phi_P1)
        phi_Pn = np.zeros_like(Phi_P1)
        Nk_P2 = Nk_SP

        '''cC3 = np.concatenate((1e-3*Nk_P0_vertex, fprop.Vp[np.newaxis,:,np.newaxis]*\
            np.ones_like(Nk_P0_vertex)),axis=-1)
        C3 = (abs(Nk_SPvertex - Nk_P0_vertex) <= np.max(abs(cC3),axis=-1)[:,:,np.newaxis])
        smooth_extrema = np.min(C3, axis=2)
        Phi_P1[Phi_P1<1] = 1*smooth_extrema[Phi_P1<1]'''

        if ctes_FR.n_points > 2:
            '---------------Hierarchical MLP limiting procedure----------------'

            #axis=1 due to reshaping of the argument [~phi_Pn.astype(bool)]
            phi_Pn = self.troubled_cell_marker(M, fprop, Nk_P1_vertex, Nk_P0_vertex, Nk_SP, machine_error)
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

                phi_P2 = self.troubled_cell_marker(M, fprop, Nk_P1_vertex2, Nk_P0_vertex2, Nk_P2, machine_error)

            phi_P2[phi_Pn==1] = 1
            Phi_P1[phi_P2==1] = 1

        Nk_SPlim = (Nk_P0 + Phi_P1[:,:,np.newaxis] * (Nk_P1 - Nk_P0) + \
            phi_P2[:,:,np.newaxis] * ((Nk_P2 - Nk_P1) + phi_Pn[:,:,np.newaxis] * (Nk_SP - Nk_P2)))

        '''C_high_order_check1 = (Nk_SPlim[:,:,inds] <= np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]) #<= machine_error
        C_high_order_check2 = (Nk_SPlim[:,:,inds] >= np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]) #>= -machine_error
        C_high_order_check1[abs((Nk_SPlim[:,:,inds] - np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])) <= machine_error] = True
        C_high_order_check2[abs((Nk_SPlim[:,:,inds] - np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])) <= machine_error] = True

        high_order_check = C_high_order_check1 * C_high_order_check2
        #if any(np.min(1 * high_order_troubled_cells,axis=2)[phi_Pn==1]==0):import pdb; pdb.set_trace()
        phi_Pn_check = np.min(1 * high_order_check,axis = 2)#[phi_Pn==1]
        if len(phi_Pn_check[phi_Pn_check==0])>0: import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()'''
        #if self.t>0.8: import pdb; pdb.set_trace()

        return Nk_SPlim

    def projections(self, Nk, m):
        Nkm = np.zeros_like(Nk)
        Nkm[:,:,0:m+1] = Nk[:,:,0:m+1]
        Nk_Pm = (ctes_FR.V[np.newaxis,] @ Nkm[:,:,:,np.newaxis])[:,:,:,0]
        return Nk_Pm

    def P1_limiter(self, M, Nk, Nk_Pm, Nk_P0_vertex, Nk_P1_vertex, machine_error):
        Linear_term = Nk_P1_vertex - Nk_P0_vertex

        'Neigboring vertex points values'
        Nk_faces = np.empty((ctes.n_components,ctes.n_internal_faces,2))
        Nk_neig = np.copy(Nk_faces)
        Nk_faces[:,:,1] = Nk_Pm[:,:,0][:,ctes_FR.v0[:,1]]
        Nk_faces[:,:,0] = Nk_Pm[:,:,-1][:,ctes_FR.v0[:,0]]
        Nk_neig[:,:,0] = Nk_P0_vertex[:,ctes_FR.v0[:,0],0]
        Nk_neig[:,:,1] = Nk_P0_vertex[:,ctes_FR.v0[:,1],0]

        #Phi_r = self.MLP_u1_mod(Nk_neig, Nk_P0_vertex, Linear_term)
        #Phi_r = self.MLP_u1(Nk_neig, Nk_P0_vertex, Linear_term)
        Phi_r = self.MLP_u2_mod(M, Nk_neig, Nk_P0_vertex, Linear_term)

        Phi_P1 = np.ones_like(Phi_r)
        Phi_P1[abs(Nk_P1_vertex - Nk_P0_vertex) >= machine_error] = Phi_r[abs(Nk_P1_vertex - Nk_P0_vertex) >= machine_error]
        Phi_P1 = np.min(Phi_P1, axis = 2)
        #Phi_P1[:,-1] = 1
        #Phi_P1[:,0] = 1

        return Phi_P1

    def MLP_u1_mod(self, Nk_neig, Nk_P0_vertex, Linear_term):
        Nk_avg_vertex = Nk_P0_vertex[:,ctes_FR.v0,0].sum(axis=-1)/2
        Nk_avg_vertex_vols = Nk_avg_vertex[:,ctes_FR.vols_vec]
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
        #Phi_r_u1[Phi_r_u1<0] = 0
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

    def MLP_u2_mod(self, M, Nk_neig, Nk_P0_vertex, Linear_term):
        delta_minus = np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec] - Nk_P0_vertex
        delta_minus2 = np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec] - Nk_P0_vertex
        delta_minus[delta_minus==0] = delta_minus2[delta_minus==0]
        delta_plus = Linear_term
        x_vols = M.data['centroid_volumes'][0,0]
        dx_vols = x_vols*2
        K = 8
        e2 = (K * dx_vols**3)
        Phi_r_u2 = 1/delta_minus * ((delta_plus**2 + e2)*delta_minus + 2*delta_minus**2*delta_plus) / \
                    (delta_plus**2 + 2*delta_minus**2 + delta_minus*delta_plus + e2)
        Phi_r_u2[Phi_r_u2>1] = 1
        Phi_r_u2[Phi_r_u2<0] = 0
        return Phi_r_u2

    def Phi_r_u2(self, M, r, machine_error):
        Phi_r_u2 = (r**2 + 2*r + machine_error) / (r**2 + r + 2 + machine_error)
        return Phi_r_u2

    def troubled_cell_marker(self, M, fprop, Nk_P1_vertex, Nk_P0_vertex, Nk_Pm, machine_error):
        Linear_term = Nk_P1_vertex - Nk_P0_vertex

        'Neigboring vertex points values'
        Nk_faces = np.empty((ctes.n_components,ctes.n_internal_faces,2))
        Nk_neig = np.empty((ctes.n_components,ctes.n_internal_faces,2))

        Nk_faces[:,:,1] = Nk_Pm[:,ctes_FR.v0[:,1], 0]
        Nk_faces[:,:,0] = Nk_Pm[:,ctes_FR.v0[:,0], -1]
        #Nk_faces_avg = 1/2 * np.sum(ctes_FR.weights[np.newaxis,np.newaxis,:,np.newaxis] * Nk_faces, axis=2)
        Nk_avg = Nk_P0_vertex[:,:,0]

        Nk_neig[:,:,0] = Nk_avg[:,ctes_FR.v0[:,0]]
        Nk_neig[:,:,1] = Nk_avg[:,ctes_FR.v0[:,1]]
        inds = np.array([0,-1])
        Nk_Pmvertex = Nk_Pm[:,:,inds]

        'P1-projected MLP condition - troubled-cell marker'

        C_troubled1 = (Nk_P1_vertex <= np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        C_troubled2 = (Nk_P1_vertex >= np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        C_troubled1[abs((Nk_P1_vertex - np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])) <= machine_error] = True
        C_troubled2[abs((Nk_P1_vertex - np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])) <= machine_error] = True
        troubled_cells = (C_troubled1 * C_troubled2)

        phi_Pm =  np.min(1 * troubled_cells, axis = 2)

        'Smooth extrema detector'
        '''High_order_term = Nk_Pmvertex - Nk_P1_vertex
        Nk_less_max = (Nk_Pmvertex < np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        Nk_big_min = (Nk_Pmvertex > np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        #Nk_less_max[abs(Nk_Pmvertex - np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]) < machine_error] = True
        #Nk_big_min[abs(Nk_Pmvertex - np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]) < machine_error] = True
        C1 = (Linear_term > 0) * (High_order_term < 0) * \
            (Nk_big_min)
        C2 = (Linear_term < 0) * (High_order_term > 0) * \
            (Nk_less_max)
        cC3 = np.concatenate((1e-3*abs(Nk_P0_vertex), fprop.Vp[np.newaxis,:,np.newaxis]*\
            np.ones_like(Nk_P0_vertex)),axis=-1)
        C3 = (abs(Nk_Pmvertex - Nk_P0_vertex) <= np.max(cC3,axis=-1)[:,:,np.newaxis])
        smooth_extrema = (C1 + C2) + C3
        smooth_extrema_cell = np.min(1*smooth_extrema,axis=-1)

        #jump = Nk_faces[:,:,:,1] - Nk_faces[:,:,:,0]

        jump = Nk_faces[:,:,1] - Nk_faces[:,:,0]
        c_vols = M.data['centroid_volumes'][:,1]
        c_faces = c_vols[ctes.v0]
        hf = c_vols[0]*2 #abs(c_faces[:,1] - c_faces[:,0])

        smooth_bound_indicator = abs(jump)/(abs(Nk_faces[:,:,1] + Nk_faces[:,:,0])/2) #em 1D apenas !!
        theta_f = smooth_bound_indicator/(hf**((ctes_FR.n_points)/2))
        trbl_bound_detector = np.empty_like(theta_f)

        trbl_bound_detector[theta_f<1] = 0 #normal bound
        trbl_bound_detector[theta_f>=1] = 1 #troubled bound
        trbl_bound_detector_vols = np.empty_like(ctes_FR.vols_vec)
        trbl_bound_detector_vols = trbl_bound_detector[:,ctes_FR.vols_vec]
        Nf = np.sum(trbl_bound_detector_vols, axis = -1)


        'Type I trouble'

        aux_I = phi_Pm[phi_Pm==1]
        aux_Nf_I = Nf[phi_Pm==1]
        aux_I[aux_Nf_I>=1] = 0
        phi_Pm[phi_Pm==1] = aux_I

        'Type II trouble'

        aux_II = smooth_extrema_cell[smooth_extrema_cell==1]
        aux_Nf_II = Nf[smooth_extrema_cell==1]
        aux_II[aux_Nf_II>=1] = 0
        smooth_extrema_cell[smooth_extrema_cell==1] = aux_II
        phi_Pm[~phi_Pm.astype(bool)] = smooth_extrema_cell[~phi_Pm.astype(bool)]'''
        return phi_Pm
