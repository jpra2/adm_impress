import numpy as np
from .properties_calculation import PropertiesCalc
from ..IMPEC.flux_calculation import Flux, RiemannSolvers
from ...directories import data_loaded
from packs.utils import constants as ctes
import scipy.sparse as sp
import time


class saturation:
    def __init__(self,M):
        c_int = M.faces.center(M.faces.internal)
        c_vols = M.volumes.center(M.volumes.all)
        pos = (c_int[:,np.newaxis,:] - c_vols[ctes.v0]).sum(axis=2)
        self.v0 = np.copy(ctes.v0)
        self.v0[:,0] = ctes.v0[pos>0]
        self.v0[:,1] = ctes.v0[pos<0]

    def TVD_auxiliary_matrix(self, coef0, coef1, coef2):
        lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        data = np.array([np.ones(len(ctes.v0[:, 0])) * coef2, np.ones(len(ctes.v0[:, 0])) * coef0,
                        np.zeros(len(ctes.v0[:, 0])), np.zeros(len(ctes.v0[:, 0]))]).flatten()
        all_neig = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes)).toarray()
        all_neig = all_neig.astype(float)
        all_neig = all_neig + np.identity(ctes.n_volumes) * coef1
        allneig_weights = all_neig.astype(float)
        return allneig_weights

    def TVD_auxiliary_matrixes(self, M):

        neig_vols = M.volumes.bridge_adjacencies(M.volumes.all,2,3)
        matriz = np.zeros((ctes.n_volumes,ctes.n_volumes))


        allneig_weights_positive = self.TVD_auxiliary_matrix(-1/6,5/6,2/6)
        allneig_weights_negative = self.TVD_auxiliary_matrix(2/6,5/6,-1/6)

        pos_neig = M.data['centroid_volumes'].T[:,np.newaxis,:] * np.sign(abs(allneig_weights_positive))[np.newaxis,:,:]

        pos = pos_neig.transpose(0,2,1)

        ds = pos_neig - pos
        ds_norm = np.linalg.norm(ds, axis=0)
        self.versor_ds = np.empty(ds.shape)

        self.versor_ds[:,ds_norm==0] = 0
        self.versor_ds[:,ds_norm!=0] = ds[:,ds_norm!=0] / ds_norm[ds_norm!=0]
        aux = self.versor_ds[abs(ds).sum(axis=-1)>0]
        aux[ds[abs(ds).sum(axis=-1)>0]==0] = 1
        self.versor_ds[abs(ds).sum(axis=-1)>0] = aux #só um auxiliar
        self.faces_dir = abs(M.faces.normal(M.faces.internal))
        return allneig_weights_positive, allneig_weights_negative

    def identify_contour_faces(self):
        lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        data = np.array([np.ones(len(ctes.v0[:, 0])), np.ones(len(ctes.v0[:, 0])),
                        np.zeros(len(ctes.v0[:, 0])), np.zeros(len(ctes.v0[:, 0]))]).flatten()
        all_neig = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes)).toarray()
        all_neig = all_neig.astype(int)
        all_neig = all_neig.sum(axis=1)
        vols_contour = np.argwhere(all_neig==1).flatten()
        faces_contour = np.empty_like(vols_contour)

        for i in range(len(vols_contour)):
            try: faces_contour[i] = np.argwhere(ctes.v0[:,0] == vols_contour[i]).flatten()
            except: faces_contour[i] = np.argwhere(ctes.v0[:,1] == vols_contour[i]).flatten()
        return faces_contour

    def TVD_properties_matrix_by_axis(self, prop_vols, allneig_weights):
        prop_vols_neig = prop_vols[:,:,np.newaxis,:] * allneig_weights
        prop_vols_neig_by_axis = np.repeat(prop_vols_neig[:,np.newaxis,:,:],3,axis=1)
        prop_vols_neig_by_axis = prop_vols_neig_by_axis * abs(self.versor_ds)[np.newaxis,:,np.newaxis,...] #axis==1 : x,y,z
        return prop_vols_neig_by_axis

    def third_order_TVD(self, M, wells, prop_vols, mobilities, Pot_hidj):
        #vec_coefs = np.array([-1/6,5/6,2/6])
        ph = mobilities.sum(axis=-1)

        blocks_upw = np.sum(mobilities[ph!=0,:]<1e-16,axis=0, dtype=bool)#will be upwinded
        faces_upw = blocks_upw[ctes.v0].sum(axis=-1,dtype=bool)
        allneig_weights_positive, allneig_weights_negative = self.TVD_auxiliary_matrixes(M)

        Pot_hidj_faces = Pot_hidj[:,ctes.v0[:,0]] - Pot_hidj[:,ctes.v0[:,1]]
        prop_vols_neig_pos_by_axis = self.TVD_properties_matrix_by_axis(prop_vols, allneig_weights_positive)
        prop_vols_neig_neg_by_axis = self.TVD_properties_matrix_by_axis(prop_vols, allneig_weights_negative)
        prop_faces_neig_by_axis_back = prop_vols_neig_pos_by_axis[:,:,:,self.v0[:,0]] * \
                self.faces_dir.T[np.newaxis,:,np.newaxis,:,np.newaxis]
        prop_faces_neig_by_axis_front = prop_vols_neig_neg_by_axis[:,:,:,self.v0[:,1]] * \
                self.faces_dir.T[np.newaxis,:,np.newaxis,:,np.newaxis]
        prop_faces_neig_back = prop_faces_neig_by_axis_back.sum(axis=1)
        prop_faces_neig_front = prop_faces_neig_by_axis_front.sum(axis=1)

        prop_internal_faces = np.ones_like(prop_faces_neig_back)
        faces_contour = self.identify_contour_faces()
        faces_upw[faces_contour] = True
        prop_internal_faces[:,:,~faces_upw] = prop_faces_neig_back[:,:,~faces_upw] * (Pot_hidj_faces>0)[np.newaxis,:,~faces_upw,np.newaxis] + \
                prop_faces_neig_front[:,:,~faces_upw] * (Pot_hidj_faces<=0)[np.newaxis,:,~faces_upw,np.newaxis]
        prop_internal_faces = prop_internal_faces.sum(axis=-1)
        faces_upw[np.sum((prop_internal_faces<0) + (prop_internal_faces>1), axis=1)[0]] = True
        prop_internal_faces[:,:,faces_upw] = prop_vols[:,:,ctes.v0[faces_upw,0]] * (Pot_hidj_faces>0)[:,faces_upw] + \
            prop_vols[:,:,ctes.v0[faces_upw,1]] * (Pot_hidj_faces<=0)[:,faces_upw]

        return prop_internal_faces

    def update_pore_volume(self, P):
        Vp = PropertiesCalc().update_porous_volume(P)
        return Vp

    def update_relative_perm(self, fprop, Sj):
        So = Sj[0]
        Sg = Sj[1]
        if ctes.n_phases == 3: Sw = Sj[2]
        else: Sw = fprop.Sw
        krs_new = PropertiesCalc().update_relative_permeabilities(fprop, So, Sg, Sw)
        return krs_new

    def properties_faces_upwind(self, Pot_hid, properties):
        #Pot_hid = P_new + fprop.Pcap - G[0,:,:]
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]

        properties_internal_faces = np.zeros([properties.shape[0], ctes.n_phases, ctes.n_internal_faces])
        properties_vols = properties[...,ctes.v0[:,0]]
        properties_vols_up = properties[...,ctes.v0[:,1]]
        properties_internal_faces[:,Pot_hidj_up <= Pot_hidj] = properties_vols[:,Pot_hidj_up <= Pot_hidj]
        properties_internal_faces[:,Pot_hidj_up > Pot_hidj] = properties_vols_up[:,Pot_hidj_up > Pot_hidj]
        return properties_internal_faces

    def update_well_term(self, wells, fprop, qk, mobilities):
        wp = wells['ws_p']
        if len(wp)>0:
            mob_ratio = mobilities[:,:,wp] / np.sum(mobilities[:,:,wp], axis = 1)
            qj = mob_ratio * np.sum(fprop.q_phase, axis=1)
            qk[:,wp] = np.sum(fprop.xkj[:,:,wp] * fprop.Csi_j[:,:,wp] * qj, axis = 1)
        return qk

    def update_well_term_der(self, wells, fprop, dqk_dSj, krs, mobilities_new, Sj):
        wp = wells['ws_p']

        if len(wp)>0:
            if ctes.n_phases==2:
                Sj = np.concatenate((Sj,fprop.Sw[np.newaxis,:]),axis=0)
                krs = np.concatenate((krs,np.zeros((1,1,ctes.n_volumes))), axis=1)

            dkrsdSj = PropertiesCalc().relative_permeability_derivative_call(krs[:,:,wp], Sj[:,wp])
            if not ctes.load_w: dkrsdSj = dkrsdSj[:,:-1,:-1,:]

            dfrj_new = (1/self.mis[:,:,np.newaxis,wp] * dkrsdSj *
                    np.sum(mobilities_new[...,wp], axis=1)[0] -
                    (mobilities_new[:,:,np.newaxis,wp] *
                    np.sum(1/self.mis[:,:,np.newaxis,wp] *
                    dkrsdSj, axis=1)[:,np.newaxis,...])) *\
                    (1/(np.sum(mobilities_new[...,wp], axis=1)[0])**2)
            dqs_dSj = dfrj_new * np.sum(fprop.q_phase, axis=1)
            dqk_dSj[:,:,wp] = np.sum(fprop.xkj[:,:,np.newaxis, wp] * fprop.Csi_j[:,:,np.newaxis, wp]
                * dqs_dSj[:,...], axis = 1)
        return dqk_dSj

    def dfractional_flux(self, dkrsdSj, phase_viscosity_internal_faces, mobilities_internal_faces):
        dfrk_faces_dSj_vols = (1/phase_viscosity_internal_faces[:,:,np.newaxis,:] *
                dkrsdSj * np.sum(mobilities_internal_faces, axis=1)[0] -
                (mobilities_internal_faces[:,:,np.newaxis,:] *
                np.sum(1/phase_viscosity_internal_faces[:,:,np.newaxis,:] *
                dkrsdSj, axis=1)[:,np.newaxis,:,:])) *\
                (1/(np.sum(mobilities_internal_faces, axis=1)[0])**2)
        return dfrk_faces_dSj_vols

    def dFk_faces_upw(self, fprop, Pot_hidj, Sj, krs_internal_faces, phase_viscosity_internal_faces,
        mobilities_internal_faces, Csi_j_internal_faces, xkj_internal_faces, ftotal):

        Pot_hid = Pot_hidj[:,ctes.v0[:,0]][:,np.newaxis,:] * np.ones((ctes.n_phases,ctes.n_phases,ctes.n_internal_faces))
        Pot_hid_up = Pot_hidj[:,ctes.v0[:,1]][:,np.newaxis,:] * np.ones((ctes.n_phases,ctes.n_phases,ctes.n_internal_faces))

        Sj_internal_faces = self.properties_faces_upwind(Pot_hidj, Sj[np.newaxis,:])
        if not ctes.load_w:
            w_int = np.zeros(ctes.n_internal_faces)
            Sj_internal_faces = np.concatenate((Sj_internal_faces,w_int[np.newaxis,np.newaxis,:]),axis=1)
            krs_internal_faces = np.concatenate((krs_internal_faces, w_int[np.newaxis,np.newaxis,:]),axis=1)
            dkrsdSj_upw = PropertiesCalc().relative_permeability_derivative_call(krs_internal_faces,
                Sj_internal_faces[0])
            dkrsdSj_upw = dkrsdSj_upw[:,:-1,:-1,:]
        else:
            dkrsdSj_upw = PropertiesCalc().relative_permeability_derivative_call(krs_internal_faces,
                Sj_internal_faces[0])

        #dkrsdSj_opp = np.zeros_like(dkrsdSj_upw)
        dfrs_faces_dSj_vols = np.zeros((1,ctes.n_phases,ctes.n_phases,ctes.n_internal_faces,2))

        dfrs_faces_dSj_vols_upw = self.dfractional_flux(dkrsdSj_upw, phase_viscosity_internal_faces, mobilities_internal_faces)
        #dfrs_faces_dSj_vols_opp = self.dfractional_flux(dkrsdSj_opp, phase_viscosity_internal_faces, mobilities_internal_faces)

        dfrs_faces_dSj_vols[:,Pot_hid>=Pot_hid_up,0] = dfrs_faces_dSj_vols_upw[:,Pot_hid>=Pot_hid_up]
        dfrs_faces_dSj_vols[:,Pot_hid<Pot_hid_up,1] = dfrs_faces_dSj_vols_upw[:,Pot_hid<Pot_hid_up]
        dFj_internal_faces = dfrs_faces_dSj_vols * (ftotal[:,:,np.newaxis]) #lembrar de quando inserir Pcap, voltar pra eq. original

        dFk_internal_faces = Flux().update_Fk_internal_faces(
            xkj_internal_faces[:,:,np.newaxis,:,np.newaxis],
            Csi_j_internal_faces[:,:,np.newaxis,:,np.newaxis], dFj_internal_faces)
        return dFk_internal_faces

    def dFk_faces_TVD(self, fprop, wells, Pot_hidj, Sj, krs_new, phase_viscosity_internal_faces,
        mobilities_internal_faces, xkj_internal_faces, Csi_j_internal_faces, ftotal):
        #dkrjdSj = diag(dkrsdSj) # ver como obter
        if ctes.n_phases==2:
            Sj = np.concatenate((Sj,fprop.Sw),axis=-1)

        Pot_hidj_faces = Pot_hidj[:,ctes.v0[:,0]] - Pot_hidj[:,ctes.v0[:,1]]
        weights_faces_neig = np.ones((ctes.n_internal_faces, 2, 2))
        weights_faces_neig[:,:,0] *= (np.array([5/6,2/6])[np.newaxis,:])
        weights_faces_neig[:,:,1] *= (np.array([2/6,5/6])[np.newaxis,:])

        dkrsdSj = PropertiesCalc().relative_permeability_derivative_call(krs_new, Sj)
        if not ctes.load_w: dkrsdSj=dkrsdSj[:,:-1,:-1,:]
        dkrsdSj_faces_neig = dkrsdSj[:,:,:,self.v0,np.newaxis] * weights_faces_neig[np.newaxis,np.newaxis,np.newaxis,]
        dkrsdSj_faces_neig_back = dkrsdSj_faces_neig[...,0]
        dkrsdSj_faces_neig_front = dkrsdSj_faces_neig[...,1]

        dkrsdSj_faces_neig = dkrsdSj_faces_neig_back * (Pot_hidj_faces>0)[np.newaxis,:,np.newaxis,:,np.newaxis] + \
                dkrsdSj_faces_neig_front * (Pot_hidj_faces<=0)[np.newaxis,:,np.newaxis,:,np.newaxis]

        dfrk_faces_dSj_vols = (1/phase_viscosity_internal_faces[:,:,np.newaxis,:,np.newaxis] *
                dkrsdSj_faces_neig * np.sum(mobilities_internal_faces, axis=1)[0,...,np.newaxis] -
                (mobilities_internal_faces[:,:,np.newaxis,:,np.newaxis] *
                np.sum(1/phase_viscosity_internal_faces[:,:,np.newaxis,:,np.newaxis] *
                dkrsdSj_faces_neig, axis=1)[:,np.newaxis,:,:])) *\
                (1/(np.sum(mobilities_internal_faces, axis=1)[0,...,np.newaxis])**2)

        dFj_internal_faces = dfrk_faces_dSj_vols * (ftotal[:,:,np.newaxis]) #lembrar de quando inserir Pcap, voltar pra eq. original
        dFk_internal_faces = Flux().update_Fk_internal_faces(
            xkj_internal_faces[:,:,np.newaxis, :,np.newaxis],
            Csi_j_internal_faces[:,:,np.newaxis, :,np.newaxis], dFj_internal_faces)
        return dFk_internal_faces

    def update_dFk_vols_total(self, fprop, dFk_internal_faces):

        ''' Function to calculate phase flux '''

        dFk_vols_total = np.empty((ctes.n_components, ctes.n_phases, ctes.n_volumes))

        for i in range(ctes.n_phases):
            cx = np.arange(ctes.n_components)
            lines = np.array([np.repeat(cx,len(ctes.v0[:,0])), np.repeat(cx,len(ctes.v0[:,1]))]).astype(int).flatten()
            cols = np.array([np.tile(ctes.v0[:,0],ctes.n_components), np.tile(ctes.v0[:,1], ctes.n_components)]).flatten()
            data = np.array([-dFk_internal_faces[:,i,:,0], dFk_internal_faces[:,i,:,1]]).flatten()
            dFk_vols_total[:,i] = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_components, ctes.n_volumes)).toarray()

        return dFk_vols_total

    def implicit_solver(self, M, fprop, wells, Pot_hid, ftotal, dVjdNk, dVjdP, P_old, qk, delta_t):
        P_new = fprop.P
        Pot_hid += P_new - P_old
        Vp_new = self.update_pore_volume(P_new)
        q_total = fprop.q_phase.sum(axis=1) #for the producer well

        ponteiro_j = np.ones_like(dVjdP[0,:2,:],dtype=bool)

        if ctes.n_phases==3:
            Sj_new = np.array([fprop.So, fprop.Sg, fprop.Sw]) #old sats as initial estimatives bruno's paper IMPSAT
        elif ctes.n_phases==2:
            Sj_new = np.array([fprop.So, fprop.Sg])

        #mobilities_new = np.empty_like(fprop.mobilities)
        #cte viscosity
        self.mis = PropertiesCalc().update_phase_viscosities(fprop, fprop.Csi_j, fprop.xkj)

        'Initialization'
        mobilities_new = np.copy(fprop.mobilities)
        mobilities_internal_faces_new = np.copy(fprop.mobilities_internal_faces)
        phase_viscosity_internal_faces = self.properties_faces_upwind(Pot_hid, self.mis)
        dqk = np.zeros((ctes.n_components,ctes.n_phases, ctes.n_volumes))
        #xkj_internal_faces= self.third_order_TVD(M, wells, fprop.xkj, mobilities_new, Pot_hid)

        xkj_internal_faces = self.properties_faces_upwind(Pot_hid, fprop.xkj)
        Csi_j_internal_faces = self.properties_faces_upwind(Pot_hid,fprop.Csi_j)

        qk_new = np.copy(qk)

        Sj_old = np.copy(Sj_new) * 1e100
        j=0

        while np.max(abs(Sj_new-Sj_old))>1e-5:
            j+=1
            Sj_old = np.copy(Sj_new)
            krs_new = self.update_relative_perm(fprop, np.copy(Sj_old))
            mobilities_new = krs_new/self.mis
            #krs_internal_faces_new = self.third_order_TVD(M, wells, krs_new, mobilities_new, Pot_hid)
            krs_internal_faces_new = self.properties_faces_upwind(Pot_hid, krs_new)

            mobilities_internal_faces_new = krs_internal_faces_new / \
                phase_viscosity_internal_faces
            Fk_vols_total_new = Flux().update_flux(M, fprop, ftotal,
                                 fprop.rho_j_internal_faces,
                                 mobilities_internal_faces_new)

            qk_new = self.update_well_term(wells, fprop, np.copy(qk), mobilities_new)

            dFk_internal_faces = self.dFk_faces_upw(fprop, Pot_hid, Sj_new, krs_internal_faces_new,
                phase_viscosity_internal_faces, mobilities_internal_faces_new, Csi_j_internal_faces,
                xkj_internal_faces, ftotal)
            #dFk_internal_faces = self.dFk_faces_TVD(fprop, wells, Pot_hid, Sj_new, krs_new,
            #    phase_viscosity_internal_faces, mobilities_internal_faces_new, xkj_internal_faces,
            #    Csi_j_internal_faces, ftotal)

            dFk_vols_total_new = self.update_dFk_vols_total(fprop, dFk_internal_faces)
            dqk_new = self.update_well_term_der(wells, fprop, np.copy(dqk), krs_new, mobilities_new, Sj_old)

            Rj = (Sj_old[:-1] * Vp_new - fprop.Vj[0,:-1] - dVjdP[0,:-1] * (P_new - P_old) - \
                delta_t * np.sum(dVjdNk[:,:-1] * (Fk_vols_total_new + qk_new)[:,np.newaxis], axis=0))#[:,ctes.vols_no_wells]

            dRj = (Vp_new - delta_t * np.sum(dVjdNk[:,:-1] *
                (dFk_vols_total_new[:,:-1] + dqk_new[:,:-1]), axis=0)) #[:,ctes.vols_no_wells]

            Sj_new[:-1,:] = Sj_old[:-1,:] - Rj/dRj
            #Sj_new[:2][(Sj_new[:2]>Sjmax)] = (Sjmax)[(Sj_new[:2]>Sjmax)]#/2
            #if any(abs(Sj_new.flatten())>1): import pdb; pdb.set_trace()
            Sj_new[-1,:] = 1 - Sj_new[:-1,:].sum(axis=0)
            if j>600: import pdb; pdb.set_trace()
            #import pdb; pdb.set_trace()

        So = Sj_new[0,:]
        Sg = Sj_new[1,:]
        Sj_new[Sj_new<np.finfo(float).eps] = 0

        if ctes.n_phases==3:
            Sw = Sj_new[2,:]
        else: Sw = fprop.Sw
        if any(Sj_new.flatten()<0): import pdb; pdb.set_trace()
        #if (Fk_vols_total_new + qk_new)[-1,-1]!=0: import pdb; pdb.set_trace()
        #fprop.Fk_vols_total = Fk_vols_total_new
        #import pdb; pdb.set_trace()
        if any(np.isnan(Sj_new.flatten())): import pdb; pdb.set_trace()
        '''krs_new = self.update_relative_perm(fprop, Sj_new)
        krs_internal_faces_new = self.properties_faces_upwind(Pot_hid, krs_new)
        phase_viscosity_internal_faces = self.properties_faces_upwind(Pot_hid, self.mis)
        mobilities_new = krs_new/self.mis
        mobilities_internal_faces_new = krs_internal_faces_new/phase_viscosity_internal_faces
        Fk_vols_total_new, wave_velocity = Flux().update_flux(M, fprop, ftotal,
                             fprop.rho_j_internal_faces, mobilities_internal_faces_new)
        Sj = -(- Vj[0] - dVjdP[0] * (P_new - P_old) - delta_t * np.sum(dVjdNk * (Fk_vols_total_new + qk_new)[:,np.newaxis,:], axis=0))/Vp_new
        '''
        ponteiro = np.ones_like(ftotal[0],dtype=bool)
        Vpm = fprop.Vp[ctes.v0]
        Nk_face = fprop.Nk[:,ctes.v0]
        P_face = fprop.P[ctes.v0].sum(axis=-1)/2
        P_face = np.concatenate((P_face[:,np.newaxis],P_face[:,np.newaxis]), axis=1)
        #import pdb; pdb.set_trace()
        wave_velocity, eigvec_m = RiemannSolvers(ctes.v0,ctes.pretransmissibility_internal_faces).\
            medium_wave_velocity(M, fprop, Nk_face, P_face, ftotal, ponteiro)
        #wave_velocity = Fk_vols_total_new/fprop.Nk
        return So, Sg, Sw, Fk_vols_total_new, wave_velocity, qk_new, mobilities_new
