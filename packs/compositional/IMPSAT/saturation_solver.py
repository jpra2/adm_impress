import numpy as np
from .properties_calculation import PropertiesCalc
from ..IMPEC.flux_calculation import Flux
from ...directories import data_loaded
from packs.utils import constants as ctes
import scipy.sparse as sp
import time


class saturation:
    def __init__(self, M):
        c_int = M.faces.center(M.faces.internal)
        c_vols = M.volumes.center(M.volumes.all)
        pos = (c_int[:,np.newaxis,:] - c_vols[ctes.v0]).sum(axis=2)
        self.v0 = np.copy(ctes.v0)
        self.v0[:,0] = ctes.v0[pos>0]
        self.v0[:,1] = ctes.v0[pos<0]

        neig_vols = M.volumes.bridge_adjacencies(M.volumes.all,2,3)
        matriz = np.zeros((ctes.n_volumes,ctes.n_volumes))

        lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        data = np.array([np.ones(len(ctes.v0[:, 0])) * (2/6), np.ones(len(ctes.v0[:, 0])) * (-1/6),
                        np.zeros(len(ctes.v0[:, 0])), np.zeros(len(ctes.v0[:, 0]))]).flatten()
        all_neig = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes)).toarray()
        all_neig = all_neig.astype(float)
        all_neig2 = all_neig + np.identity(ctes.n_volumes) * (5/6)
        self.allneig_weights = all_neig2.astype(float)

        pos_neig = M.data['centroid_volumes'].T[:,np.newaxis,:] * np.sign(abs(all_neig2))[np.newaxis,:,:]

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

    def third_order_TVD(self, wells, prop_vols, Pot_hidj):
        #vec_coefs = np.array([-1/6,5/6,2/6])
        Pot_hidj_faces = Pot_hidj[:,ctes.v0[:,0]] - Pot_hidj[:,ctes.v0[:,1]]
        prop_vols_neig = prop_vols[:,:,np.newaxis,:] * self.allneig_weights
        prop_vols_neig_by_axes = np.repeat(prop_vols_neig[:,np.newaxis,:,:],3,axis=1)
        prop_vols_neig_by_axes = prop_vols_neig_by_axes * abs(self.versor_ds)[np.newaxis,:,np.newaxis,...] #axis==1 : x,y,z
        aux_wells = prop_vols_neig_by_axes[:,:,:,wells['all_wells']]
        aux_wells[aux_wells!=0] = 1/2 * aux_wells[aux_wells!=0]
        prop_vols_neig_by_axes[:,:,:,wells['all_wells']] = aux_wells
        prop_faces_neig_by_axes_back = prop_vols_neig_by_axes[:,:,:,self.v0[:,0]] * \
                self.faces_dir.T[np.newaxis,:,np.newaxis,:,np.newaxis]
        prop_faces_neig_by_axes_front = prop_vols_neig_by_axes[:,:,:,self.v0[:,1]] * \
                self.faces_dir.T[np.newaxis,:,np.newaxis,:,np.newaxis]
        prop_faces_neig_back = prop_faces_neig_by_axes_back.sum(axis=1)
        prop_faces_neig_front = prop_faces_neig_by_axes_front.sum(axis=1)
        prop_internal_faces_TVD = prop_faces_neig_back * (Pot_hidj_faces>0)[np.newaxis,:,:,np.newaxis] + \
                prop_faces_neig_front * (Pot_hidj_faces<0)[np.newaxis,:,:,np.newaxis]
        prop_internal_faces_TVD = prop_internal_faces_TVD.sum(axis=-1)
        return prop_internal_faces_TVD

    def update_pore_volume(self, P):
        Vp = PropertiesCalc().update_porous_volume(P)
        return Vp

    def update_relative_perm(self, fprop, Sj):

        So = Sj[0]
        Sg = Sj[1]
        Sw = Sj[2]

        krs_new = PropertiesCalc().update_relative_permeabilities(fprop, So, Sg, Sw)
        return krs_new

    def properties_faces_upwind(self, Pot_hid, properties):
        #Pot_hid = P_new + fprop.Pcap - G[0,:,:]
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]

        properties_internal_faces = np.zeros([properties.shape[0], ctes.n_phases, ctes.n_internal_faces])
        properties_vols = properties[:,:,ctes.v0[:,0]]
        properties_vols_up = properties[:,:,ctes.v0[:,1]]
        properties_internal_faces[0,Pot_hidj_up <= Pot_hidj] = properties_vols[0,Pot_hidj_up <= Pot_hidj]
        properties_internal_faces[0,Pot_hidj_up > Pot_hidj] = properties_vols_up[0,Pot_hidj_up > Pot_hidj]
        return properties_internal_faces


    def update_well_term(self, wells, fprop, qk, mobilities):
        wp = wells['ws_p']
        if len(wp)>0:
            mob_ratio = mobilities[:,:,wp] / np.sum(mobilities[:,:,wp], axis = 1)
            qj = mob_ratio * np.sum(fprop.q_phase, axis=1)
            qk[:,wp] = np.sum(fprop.xkj[:,:,wp] * fprop.Csi_j[:,:,wp] * qj, axis = 1)
        return qk

    def update_well_term_der(self, wells, fprop, dqk_dSj, dkrsdSj, dfrj_new):
        wp = wells['ws_p']
        if len(wp)>0:
            dqs_dSj = dfrj_new[:,:,:,wp] * np.sum(fprop.q_phase, axis=1)
            dqk_dSj[:,:,wp] = np.sum(fprop.xkj[:,:,np.newaxis, wp] * fprop.Csi_j[:,:,np.newaxis, wp]
                * dqs_dSj[:,...], axis = 1)
        return dqk_dSj


    def dfractional_flux(self, fprop, Pot_hid, saturations, krs, mobilities_new):
        #dkrjdSj = diag(dkrsdSj) # ver como obter
        #saturations_ = np.empty((3,len(saturations[0])))
        #saturations_[0:2,:] = saturations
        #saturations[-1,:] = 1 - np.sum(saturations[0:2],axis=0)

        dkrsdSj = PropertiesCalc().relative_permeability_derivative_call(krs, saturations)
        dkrsdSj[:,:,saturations==0] = 0
        dfrj_new = (1/self.mis[:,:,np.newaxis,:] * dkrsdSj * np.sum(mobilities_new, axis=1)[0] -
                (mobilities_new[:,:,np.newaxis,:] * np.sum(1/self.mis[:,:,np.newaxis,:] * dkrsdSj, axis=1)[:,np.newaxis,:,:])) *\
                (1/(np.sum(mobilities_new[:,:,np.newaxis,np.newaxis,:], axis=1))**2)
        #faces_contour = self.identify_contour_faces()
        #dfrj_new[:,:,:,ctes.v0[faces_contour].flatten()] =  0#- dfrj_new[:,:,:,ctes.v0[faces_contour].flatten()]
        dfrj_internal_faces = self.properties_faces_upwind(Pot_hid, dfrj_new[0])
        #dfrj_internal_faces[:,:,faces_contour] = 0
        return dfrj_new, dfrj_internal_faces, dkrsdSj[0]

    def update_dFk_internal_faces(self, fprop, dfrj_internal_faces_new, Ft_internal_faces, rho_j_internal_faces,
        mobilities_internal_faces, Pcap_face, z_face, pretransmissibility_internal_faces):

        ''' Function to calculate phase flux '''
        dFj_internal_faces = dfrj_internal_faces_new * (Ft_internal_faces) #lembrar de quando inserir Pcap, voltar pra eq. original

        dFk_internal_faces = Flux().update_Fk_internal_faces(
            self.xkj_internal_faces[:,:,np.newaxis, :],
            fprop.Csi_j_internal_faces[:,:,np.newaxis, :], dFj_internal_faces[np.newaxis,...])

        dFk_vols_total = np.empty((ctes.n_components, ctes.n_phases, ctes.n_volumes))
        # da pra vetorizar isso ainda, só não pensei com calma e pra nao fazer besteira, vou deixar pra depois
        for i in range(ctes.n_phases):
            cx = np.arange(ctes.n_components)
            lines = np.array([np.repeat(cx,len(ctes.v0[:,0])), np.repeat(cx,len(ctes.v0[:,1]))]).astype(int).flatten()
            cols = np.array([np.tile(ctes.v0[:,0],ctes.n_components), np.tile(ctes.v0[:,1], ctes.n_components)]).flatten()
            data = np.array([-dFk_internal_faces[:,i,:], dFk_internal_faces[:,i,:]]).flatten()
            dFk_vols_total[:,i] = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_components, ctes.n_volumes)).toarray()

        return dFk_vols_total

    def implicit_solver(self, M, fprop, wells, Pot_hid, Ft_internal_faces, dVjdNk, dVjdP, P_new, P_old, qk, delta_t):
        Pot_hid += P_new - P_old
        Vj = fprop.Nj/fprop.Csi_j
        Vp_new = self.update_pore_volume(P_new)
        #dVjdP = dVjdP[0,0:2,:] #only hydrocarbon phases
        #dVjdNk = dVjdNk[:-1,0:2,:] #only hydrocarbon phases
        q_total = fprop.q_phase.sum(axis=1) #for the producer well

        ponteiro_j = np.ones_like(dVjdP[0,:2,:],dtype=bool)

        Sj_new = np.array([fprop.So, fprop.Sg, fprop.Sw]) #old sats as initial estimatives bruno's paper IMPSAT
        #mobilities_new = np.empty_like(fprop.mobilities)

        self.xkj_internal_faces = self.properties_faces_upwind(Pot_hid, fprop.xkj)
        self.mis = PropertiesCalc().update_phase_viscosities(fprop, fprop.Csi_j, fprop.xkj)

        Sgr = float(data_loaded['compositional_data']['residual_saturations']['Sgr'])
        Swr = float(data_loaded['compositional_data']['residual_saturations']['Swr'])
        Sor = fprop.Sor
        Sjmax = np.empty((2,ctes.n_volumes))
        Sjmin = np.empty((2,ctes.n_volumes))

        Sjmax[0,:] = 1 - Sgr - Swr
        Sjmax[1,:] = 1 - Sor - Swr
        Sjmin[0,:] = Sor
        Sjmin[1,:] = Sgr
        save_Smax = np.copy(Sjmax)
        save_Smin = np.copy(Sjmin)

        'Initialization'
        mobilities_new = np.copy(fprop.mobilities)
        Fk_vols_total_new = fprop.Fk_vols_total
        mobilities_internal_faces_new = fprop.mobilities_internal_faces
        krs_new = self.mis * mobilities_new
        dqk = np.zeros((ctes.n_components,ctes.n_phases, len(fprop.mobilities[0,0,:])))
        qk_new = np.copy(qk)

        Sj_old = np.copy(Sj_new) * 1/2
        j=0

        while np.max(abs(Sj_new-Sj_old)[:,ctes.vols_no_wells])>1e-5:
            j+=1
            Sj_old[:,ctes.vols_no_wells] = np.copy(Sj_new)[:,ctes.vols_no_wells]
            krs_new[:,:,ctes.vols_no_wells] = self.update_relative_perm(fprop, np.copy(Sj_old)[:,ctes.vols_no_wells])

            krs_internal_faces_new = self.properties_faces_upwind(Pot_hid, krs_new)#self.third_order_TVD(wells, krs_new, Pot_hid)
            phase_viscosity_internal_faces = self.properties_faces_upwind(Pot_hid, self.mis)
            mobilities_new = krs_new/self.mis
            mobilities_internal_faces_new = krs_internal_faces_new/phase_viscosity_internal_faces
            Fk_vols_total_new, wave_velocity = Flux().update_flux(M, fprop, P_old, Ft_internal_faces,
                                 fprop.rho_j_internal_faces, mobilities_internal_faces_new)
            qk_new = self.update_well_term(wells, fprop, np.copy(qk), mobilities_new)

            dfrj_new, dfrj_internal_faces, dkrsdSj = self.dfractional_flux(fprop, Pot_hid, Sj_old, krs_new, mobilities_new)

            dFk_vols_total_new = self.update_dFk_internal_faces(fprop, dfrj_internal_faces, Ft_internal_faces, fprop.rho_j_internal_faces,
                mobilities_internal_faces_new, fprop.Pcap[:,ctes.v0], ctes.z[ctes.v0],
                ctes.pretransmissibility_internal_faces)
            dqk_new = self.update_well_term_der(wells, fprop, np.copy(dqk), dkrsdSj, dfrj_new)
            #dqk_new[:,:,wells['all_wells']] = - dFk_vols_total_new[:,:,wells['all_wells']]

            Rj = Sj_old[:2] * Vp_new - Vj[0,:2] - dVjdP[0,:2] * (P_new - P_old) - \
                delta_t * np.sum(dVjdNk[:,:2] * (Fk_vols_total_new + qk_new), axis=0)

            dRj = Vp_new - delta_t * np.sum(dVjdNk[:,:2] *
                (dFk_vols_total_new[:,:2] + dqk_new[:,:2]), axis=0)

            Sj_new[:2] = Sj_old[:2] - Rj/dRj

            #Sj_new[:2][(Sj_new[:2]-Sjmax)>1e-3] = (Sjmax + Sjmin)[(Sj_new[:2]-Sjmax)>1e-3]/2
            if any(abs(Sj_new.flatten())>1): import pdb; pdb.set_trace()
            #if j>00: import pdb; pdb.set_trace()
            Sj_new[-1,:] = 1 - Sj_new[:2].sum(axis=0)
            import pdb; pdb.set_trace()



        if abs(Sj_new[0,-1]-0.8)>0.01: import pdb; pdb.set_trace()

        if any(Sj_new[-1,:]<0): import pdb; pdb.set_trace()
        So = Sj_new[0,:]
        Sg = Sj_new[1,:]
        Sw = Sj_new[2,:]
        #if (Fk_vols_total_new + qk_new)[-1,-1]!=0: import pdb; pdb.set_trace()
        #fprop.Fk_vols_total = Fk_vols_total_new

        krs_new = self.update_relative_perm(fprop, Sj_new)
        krs_internal_faces_new = self.properties_faces_upwind(Pot_hid, krs_new)
        phase_viscosity_internal_faces = self.properties_faces_upwind(Pot_hid, self.mis)
        mobilities_new = krs_new/self.mis
        mobilities_internal_faces_new = krs_internal_faces_new/phase_viscosity_internal_faces
        Fk_vols_total_new, wave_velocity = Flux().update_flux(M, fprop, P_old, Ft_internal_faces,
                             fprop.rho_j_internal_faces, mobilities_internal_faces_new)
        #Sj  = -(- Vj[0] - dVjdP[0] * (P_new - P_old) - delta_t * np.sum(dVjdNk * (Fk_vols_total_new + qk_new)[:,np.newaxis,:], axis=0))/Vp_new

        return So, Sg, Sw, Fk_vols_total_new, wave_velocity, qk_new, mobilities_new
