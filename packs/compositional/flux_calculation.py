import numpy as np
from ..directories import data_loaded
from ..utils import constants as ctes
from .stability_check import StabilityCheck
from .properties_calculation import PropertiesCalc
import scipy.sparse as sp
from multiprocessing import Pool

class FOUM:
    """ Class created for computing flux accondingly to the First Order Upwind \
    Method """

    def update_flux(self, fprop, Ft_internal_faces, rho_j_internal_faces, mobilities_internal_faces):
        ''' Main function that calls others '''
        self.update_Fj_internal_faces(Ft_internal_faces,
        rho_j_internal_faces, mobilities_internal_faces, fprop.Pcap[:,ctes.v0],
        ctes.z[ctes.v0], ctes.pretransmissibility_internal_faces)

        Fk_internal_faces = self.update_Fk_internal_faces(
        fprop.xkj_internal_faces,fprop.Csi_j_internal_faces)

        self.update_flux_volumes(fprop, Fk_internal_faces)

    def update_Fj_internal_faces(self, Ft_internal_faces, rho_j_internal_faces,
        mobilities_internal_faces, Pcap_face, z_face,
        pretransmissibility_internal_faces):
        ''' Function to calculate phase flux '''

        frj = mobilities_internal_faces[0,...] / \
        np.sum(mobilities_internal_faces[0,...], axis = 0)

        self.Fj_internal_faces = frj[np.newaxis,...] * (Ft_internal_faces +
            pretransmissibility_internal_faces * (np.sum(mobilities_internal_faces *
            (Pcap_face[:,:,1] - Pcap_face[:,:,0] - ctes.g * rho_j_internal_faces *
            (z_face[:,1] - z_face[:,0])), axis=1) - np.sum(mobilities_internal_faces,
            axis=1) * (Pcap_face[:,:,1] - Pcap_face[:,:,0] - ctes.g *
            rho_j_internal_faces * (z_face[:,1] - z_face[:,0]))))
        # M.flux_faces[M.faces.internal] = Ft_internal_faces * M.faces.normal[M.faces.internal].T

    def update_Fk_internal_faces(self, xkj_internal_faces, Csi_j_internal_faces):
        ''' Function to compute component flux '''

        Fk_internal_faces = np.sum(xkj_internal_faces * Csi_j_internal_faces *
        self.Fj_internal_faces, axis = 1)
        return Fk_internal_faces

    def update_flux_volumes(self, fprop, Fk_internal_faces):
        ''' Function to compute component flux balance through the control \
        volume interfaces'''

        cx = np.arange(ctes.n_components)
        lines = np.array([np.repeat(cx,len(ctes.v0[:,0])), np.repeat(cx,len(ctes.v0[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(ctes.v0[:,0],ctes.n_components), np.tile(ctes.v0[:,1], ctes.n_components)]).flatten()
        data = np.array([-Fk_internal_faces, Fk_internal_faces]).flatten()
        fprop.Fk_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_components, ctes.n_volumes)).toarray()


class MUSCL:

    """ Class created for the second order MUSCL implementation for the \
    calculation of the advective terms """

    def run(self, M, fprop, wells, P_old, ftot, order=1):
        ''' Global function that calls others '''
        self.Ft_internal_faces = ftot
        self.P_face = np.sum(P_old[ctes.v0], axis=1) * 0.5
        dNk_vols = self.volume_gradient_reconstruction(M, fprop, wells)
        dNk_face, dNk_face_neig = self.get_faces_gradient(M, fprop, dNk_vols)

        if order == 2:
            Phi = self.Van_Leer_slope_limiter(dNk_face, dNk_face_neig)
            Nk_face, z_face = self.get_extrapolated_compositions(fprop, Phi, dNk_face_neig)
        elif order == 3:
            Phi_MLP, Theta_MLP2 = self.MLP_slope_limiter(M, fprop, np.copy(dNk_face_neig), dNk_face)
            Nk_face, z_face = self.get_extrapolated_compositions_MLP(fprop, Phi_MLP, Theta_MLP2, dNk_face_neig, dNk_face)
        else:
            Phi = np.zeros_like(dNk_face_neig)
            Nk_face, z_face = self.get_extrapolated_compositions(fprop, Phi, dNk_face_neig)
        #G = self.update_gravity_term() # for now, it has no gravity
        alpha = self.update_flux(M, fprop, Nk_face)
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
        r = (dNk_face[:,:,np.newaxis]/(dNk_face_neig))
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
        vols_vec = -np.ones((ctes.n_volumes,2),dtype=int)
        lines = np.arange(ctes.n_internal_faces)
        vols_vec[ctes.v0[:,0],0] = lines
        vols_vec[ctes.v0[:,1],1] = lines
        contour_vols = np.argwhere(vols_vec<0)[:,0]
        vols_vec[vols_vec < 0] = vols_vec[contour_vols,:][vols_vec[contour_vols]>0]

        Phi_MLPv = np.empty((ctes.n_components,ctes.n_volumes,2))
        Phi_MLPv[:,:,0] = Phi_MLP[:,vols_vec[:,0],0]
        Phi_MLPv[:,:,1] = Phi_MLP[:,vols_vec[:,1],1]
        Phi_MLPv = np.min(Phi_MLPv,axis=2)
        Phi_MLP = Phi_MLPv[:,ctes.v0]

        Theta_MLP2v = np.empty((ctes.n_components,ctes.n_volumes,2))
        Theta_MLP2v[:,:,0] = Theta_MLP2[:,vols_vec[:,0],0]
        Theta_MLP2v[:,:,1] = Theta_MLP2[:,vols_vec[:,1],1]
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

    def update_flux(self, M, fprop, Nk_face):
        Fk_internal_faces = np.empty((ctes.n_components,ctes.n_internal_faces))
        alpha_wv = np.empty((ctes.n_components, ctes.n_internal_faces,5))
        Fk_face = self.get_Fk_face(fprop, M, Nk_face)

        ponteiro = np.zeros(ctes.n_internal_faces,dtype=bool)

        #LR_eigval, d2FkdNk_eigval = self.get_LR_eigenvalues_Serna(M, fprop, Nk_face)
        #LR_eigval = self.get_LR_eigenvalues(M, fprop, Nk_face,  np.ones(ctes.n_internal_faces,dtype=bool))
        #ponteiro = self.flux_calculation_conditions_Serna(LR_eigval, d2FkdNk_eigval)
        #ponteiro = self.flux_calculation_conditions(LR_eigval)

        alpha = self.wave_velocity_LLF(M, fprop, Nk_face, np.copy(~ponteiro))
        Fk_internal_faces[:,~ponteiro], alpha_wv[:,~ponteiro,:] = \
        self.update_flux_LLF(Fk_face[:,~ponteiro,:], Nk_face[:,~ponteiro,:],
        alpha)

        '''if any(ponteiro):
            alpha_wv[:,ponteiro,:] = 0
            #alpha_wv[:,ponteiro,0:2] = LR_eigval[:,ponteiro]
            Fk_internal_faces[:,ponteiro] = self.update_flux_upwind(fprop,
                                                Fk_face[:,ponteiro,:], np.copy(ponteiro))'''

        '-------- Perform volume balance to obtain flux through volumes -------'
        FOUM().update_flux_volumes(fprop, Fk_internal_faces)
        return alpha_wv

    def get_extrapolated_properties(self, fprop, M, Nk_face, z_face, P_face, Vp, v):
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
            L_face = np.ones(len(z_face[0,:])); V_face = np.zeros(len(z_face[0,:]))
            xkj_face[0:ctes.Nc,0:2,:,i] = 1
            rho_j_face = fprop.rho_j[:,:,ctes.v0]
            Csi_j_face = fprop.Csi_j[:,:,ctes.v0]

        if ctes.load_w:
            xkj_face[-1,-1,...] = 1
            xkj_face[-1,0:-1,...] = 0
            xkj_face[0:ctes.Nc,-1,...] = 0

            Sw_face, Csi_j_face[0,-1,...], rho_j_face[0,-1,...] = \
            PropertiesCalc().update_water_saturation(fprop, Nk_face[-1,...],
            P_face, Vp)

        So_face, Sg_face =  PropertiesCalc().update_saturations(Sw_face,
        Csi_j_face, L_face, V_face)

        mobilities_face = PropertiesCalc().update_mobilities(fprop, So_face,
        Sg_face, Sw_face, Csi_j_face, xkj_face)
        return mobilities_face, rho_j_face, Csi_j_face, xkj_face

    def Fk_Nk(self, fprop, M, Nk, ponteiro):
        ''' Function to compute component flux based on a given composition (Nk) '''
        v = int(len(Nk[0,:])/len(ponteiro[ponteiro]))
        z = Nk[0:ctes.Nc] / np.sum(Nk[0:ctes.Nc], axis = 0)

        Vp_reshaped1 = np.tile(fprop.Vp[ctes.v0][ponteiro,0],(int(v/2) * \
        (np.sign(v - ctes.n_components)**2) + int(v)*(1 - np.sign(v -
        ctes.n_components)**2)))
        Vp_reshaped2 = np.tile(fprop.Vp[ctes.v0][ponteiro,1],(int(v/2) * \
        (np.sign(v - ctes.n_components)**2) + int(v/(ctes.n_components+1))*(1 -
        np.sign(v-ctes.n_components))))

        Vp_reshaped = np.concatenate((Vp_reshaped1, Vp_reshaped2))

        mobilities, rho_j, Csi_j, xkj = self.get_extrapolated_properties(fprop, M, Nk, z,
        np.tile(self.P_face[ponteiro],v), Vp_reshaped, v)
        ftotal = self.Ft_internal_faces[:,ponteiro]

        f = FOUM()

        Pcap_reshaped = np.concatenate(np.dsplit(np.tile(fprop.Pcap[:,ctes.v0][:,ponteiro],v),v),axis=1)
        z_reshaped = np.concatenate(np.hsplit(np.tile(ctes.z[ctes.v0][ponteiro],v),v),axis=0)
        f.update_Fj_internal_faces(np.tile(ftotal,v), rho_j, mobilities, Pcap_reshaped,
        z_reshaped, np.tile(ctes.pretransmissibility_internal_faces[ponteiro],v))
        Fk = f.update_Fk_internal_faces(xkj, Csi_j)
        return Fk

    def dFk_dNk(self, fprop, M, Nk_face, delta, k, ponteiro):
        Nk_face_plus = np.copy(Nk_face)
        Nk_face_minus = np.copy(Nk_face)
        Nk_face_plus[k] += delta * 0.5
        Nk_face_minus[k] -= delta * 0.5
        dFkdNk = (self.Fk_Nk(fprop, M, Nk_face_plus, np.ones(ctes.n_internal_faces, dtype=bool)) -
        self.Fk_Nk(fprop, M, Nk_face_minus, np.ones(ctes.n_internal_faces, dtype=bool)))\
        /(Nk_face_plus[k]-Nk_face_minus[k])
        return dFkdNk

    def get_Fk_face(self, fprop, M, Nk_face):
        ''' Function that computes the flux in each face side (Left and Right)'''
        Fk_face = self.Fk_Nk(fprop, M, np.concatenate((Nk_face[:,:,0],Nk_face[:,:,1]),axis=1),
        np.ones(ctes.n_internal_faces, dtype=bool))
        Fk_face = np.concatenate(np.hsplit(Fk_face[:,:,np.newaxis],2),axis = 2)
        return Fk_face

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
                dFkdNk[:,:,k] = ((self.Fk_Nk(fprop, M, Nk_face_plus, np.ones(ctes.n_internal_faces, dtype=bool)) -
                                    self.Fk_Nk(fprop, M, Nk_face_minus, np.ones(ctes.n_internal_faces, dtype=bool)))\
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
        dFkdNk = ((self.Fk_Nk(fprop, M, Nk_face_plus, ponteiro) -
        self.Fk_Nk(fprop, M, Nk_face_minus, ponteiro))/(Nk_face_plus - Nk_face_minus).sum(axis=0))
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,np.newaxis],2),axis=2)
        dFkdNk = np.concatenate(np.hsplit(dFkdNk[:,:,:,np.newaxis],ctes.n_components),axis=3)
        dFkdNk = dFkdNk.transpose(2,1,0,3)

        eigval1, v = np.linalg.eig(dFkdNk)
        dFkdNk_eigvalue = eigval1.T

        return dFkdNk_eigvalue

    def wave_velocity_LLF(self, M, fprop, Nk_face, ponteiro):
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
        dFkdNk_gauss = ((self.Fk_Nk(fprop, M, Nkg_plus, ponteiro) -
        self.Fk_Nk(fprop, M, Nkg_minus, ponteiro))/(Nkg_plus - Nkg_minus).sum(axis=0))
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
        dFkdNk_m = ((self.Fk_Nk(fprop, M, Nkm_plus, ponteiro) -
        self.Fk_Nk(fprop, M, Nkm_minus, ponteiro))/ (Nkm_plus - Nkm_minus).sum(axis=0))
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
        dFkdNk = ((self.Fk_Nk(fprop, M, Nk_face_plus, ponteiro) -
        self.Fk_Nk(fprop, M, Nk_face_minus, ponteiro))/(Nk_face_plus - Nk_face_minus).sum(axis=0))
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

    def wave_velocity_DW(self, M, fprop, Nk_face):
        np.seterr(divide='ignore', invalid='ignore')
        delta = 0.0001
        Nkm = (Nk_face[:,:,1] + Nk_face[:,:,0])/2
        Fk_face = np.empty((ctes.n_components,ctes.n_internal_faces, 2))

        dFkdNk_RLG = np.empty_like(Fk_face)
        dFkdNk_MG = np.empty_like(Fk_face)
        dFkdNk_m = np.empty((ctes.n_internal_faces, ctes.n_components, ctes.n_components))

        for i in range(2):
            f = FOUM()
            f.update_Fj_internal_faces(self.Ft_internal_faces,
                        self.rho_j_face[:,:,:,i], self.mobilities_face[:,:,:,i],
                        fprop.Pcap[:,ctes.v0], ctes.z[ctes.v0], ctes.pretransmissibility_internal_faces)
            Fk_face[:,:,i] = f.update_Fk_internal_faces(self.xkj_face[:,:,:,i],
                             self.Csi_j_face[:,:,:,i])
            if i==0:
                for k in range(ctes.n_components):

                    Nkm_plus = np.copy(Nkm)
                    Nkm_minus = np.copy(Nkm)
                    Nkm_plus[k,:] += delta/2
                    Nkm_minus[k,:] -= delta/2
                    dFkdNk_m[:,:,k] = (self.Fk_Nk(fprop, M, Nkm_plus) -
                                self.Fk_Nk(fprop, M, Nkm_minus)).T/delta
            '''Nkg = Nkm + (Nk_face[:,:,i] - Nkm)/3**(1/2)
            dFkdNk_RLG[:,:,i] = (self.Fk_Nk(fprop, M, Nk_face[:,:,i], P_old) -
                                self.Fk_Nk(fprop, M, Nkg, P_old))/(Nk_face[:,:,i] - Nkg)

            dFkdNk_MG[:,:,i] = (self.Fk_Nk(fprop, M, Nkm, P_old) -
                                self.Fk_Nk(fprop, M, Nkg, P_old))/(Nkm - Nkg)
            dFkdNk_RLG[(Nk_face[:,:,i] - Nkg)==0] = 0;
            dFkdNk_MG[(Nkm - Nkg)==0] = 0''' #THIS WAS FOR MDW SCHEME

        eigval, v = np.linalg.eig(dFkdNk_m)
        dFkdNk_m = eigval.T
        alpha = (Fk_face[:,:,1] - Fk_face[:,:,0]) / (Nk_face[:,:,1] - Nk_face[:,:,0])
        alpha[(Nk_face[:,:,1] - Nk_face[:,:,0])<=0.01] = dFkdNk_m[(Nk_face[:,:,1] - Nk_face[:,:,0])<=0.01]
        return Fk_face, alpha

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

    def update_flux_DW(self, Fk_face_DW_all, Nk_face_DW, alpha_DW):
        Fk_face_DW = 0.5*(Fk_face_DW_all.sum(axis=-1) - abs(alpha_DW) * \
                    (Nk_face_DW[:,:,1] - Nk_face_DW[:,:,0]))
        alpha = np.max(abs(alpha_DW),axis = 0)
        return Fk_face_DW, alpha

    def update_flux_LLF(self, Fk_face_LLF_all, Nk_face_LLF, alpha_LLF):
        alpha = np.max(abs(alpha_LLF),axis = 0)
        Fk_face_LLF = 0.5*(Fk_face_LLF_all.sum(axis=-1) - np.max(abs(alpha),axis=-1) * \
                    (Nk_face_LLF[:,:,1] - Nk_face_LLF[:,:,0]))
        return Fk_face_LLF, alpha

class FR:
    def __init__(self):
        'Enviroment for the FR/CPR method'
