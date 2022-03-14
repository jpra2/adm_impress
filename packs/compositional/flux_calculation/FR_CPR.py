from .riemann_solvers import RiemannSolvers
from .. import prep_FR as ctes_FR
from ..IMPEC.properties_calculation import PropertiesCalc
from ..IMPEC.composition_solver import RK3, Euler
from .MUSCL import MUSCL
from packs.utils import constants as ctes
from packs.directories import data_loaded
import scipy.sparse as sp
import numpy as np

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

    def run(self, M, fprop, wells, Ft_internal, Nk_SP_old, P_old, delta_t, t):
        Nk_SP = np.copy(Nk_SP_old)

        q_SP = fprop.qk_molar[:,:,np.newaxis] * np.ones_like(Nk_SP)
        #q_SP[:] = 0 #BASTIAN

        self.P_faces = np.sum(P_old[ctes_FR.v0],axis=-1)*0.5
        self.P_SP = self.get_pressure_SP(wells, P_old)
        self.Nk_SP_old = Nk_SP_old
        self.delta_t = delta_t

        self.Ft_SP = self.total_flux_SP(fprop, wells, Ft_internal)

        dFk_SP, wave_velocity = self.dFk_SP_from_Pspace(M, fprop, wells, Ft_internal, np.copy(Nk_SP), P_old)
        Nk_SP = RK3.update_composition_RK3_1(np.copy(Nk_SP_old), q_SP, dFk_SP, delta_t)
        Nk_SP = self.slopeLimN(M, Nk_SP)
        #Nk_SP = self.MLP_slope_limiter(M, fprop, Nk_SP, wells)

        dFk_SP, wave_velocity = self.dFk_SP_from_Pspace(M, fprop, wells, Ft_internal, np.copy(Nk_SP), P_old)
        Nk_SP = RK3.update_composition_RK3_2(np.copy(Nk_SP_old), q_SP, np.copy(Nk_SP), dFk_SP, delta_t)
        Nk_SP = self.slopeLimN(M, Nk_SP)
        #Nk_SP = self.MLP_slope_limiter(M, fprop, Nk_SP, wells)

        dFk_SP, wave_velocity = self.dFk_SP_from_Pspace(M, fprop, wells, Ft_internal, np.copy(Nk_SP), P_old)
        Nk_SP = RK3.update_composition_RK3_3(np.copy(Nk_SP_old), q_SP, np.copy(Nk_SP), dFk_SP, delta_t)
        Nk_SP = self.slopeLimN(M, Nk_SP)
        #Nk_SP = self.MLP_slope_limiter(M, fprop, Nk_SP, wells)

        Fk_vols_total = np.min(abs(dFk_SP),axis=2)
        Nk = 1 / sum(ctes_FR.weights) * np.sum(ctes_FR.weights * Nk_SP,axis=2)
        z = Nk[0:ctes.Nc,:] / np.sum(Nk[0:ctes.Nc,:], axis = 0)

        if (any((z<0).flatten())): import pdb; pdb.set_trace()
        return wave_velocity, Nk, z, Nk_SP, Fk_vols_total

    def dFk_SP_from_Pspace(self, M, fprop, wells, Ft_internal, Nk_SP, P_old):

        Fk_SP = self.component_flux_SP(fprop, M, Nk_SP)
        #Fk_SP = self.get_Fk_Ft_SP(fprop, M, Nk_SP)

        #FOR THE BURGERS PROBLEM
        #Fk_SP = (Nk_SP**2/2)/(1/ctes.n_volumes)

        #FOR THE BASTIAN PROBLEM
        '''Csi_j_SP = fprop.Csi_j[:,:,:,np.newaxis]*np.ones((1,ctes.n_phases,ctes.n_volumes,ctes_FR.n_points))
        Sw_SP = Nk_SP[-1] * (1/Csi_j_SP[0,-1]) / fprop.Vp[0]
        mobilities = np.empty_like(Csi_j_SP)
        xkj = np.zeros((ctes.n_components,ctes.n_phases,ctes.n_volumes,ctes_FR.n_points))
        xkj[0:-1,0:-1,:] = 1
        xkj[-1,-1,:] = 1
        L_face = np.ones_like(Sw_SP)
        V_face = 1 - L_face
        So_SP = np.empty_like(Sw_SP)
        Sg_SP = np.empty_like(Sw_SP)
        So_SP, Sg_SP =  PropertiesCalc().update_saturations(Sw_SP,
            Csi_j_SP, L_face, V_face)
        for i in range(ctes_FR.n_points):
            mobilities[:,:,:,i] = PropertiesCalc().update_mobilities(fprop, So_SP[:,i], Sg_SP[:,i],
                Sw_SP[:,i], Csi_j_SP[:,:,:,i], xkj[:,:,:,i])
        f = mobilities[0,...]/(np.sum(mobilities[0,...],axis=0))
        Fj_SP = f * Ft_internal[0,0] #self.Flux_SP(fprop, wells, M, Fk_faces)
        Fk_SP = np.sum(xkj * Csi_j_SP * Fj_SP, axis = 1)'''

        Fk_faces, Fk_vols_RS_neig, wave_velocity = self.Riemann_Solver(M, fprop, wells, Nk_SP,
            Fk_SP, Ft_internal)

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

        return dFk_SP, wave_velocity

    def total_flux_SP(self, fprop, wells, Ft_internal):
        'RTo'
        phi = np.empty((len(ctes_FR.points),2))
        phi[:,0] = 1 / 4 * (1 + ctes_FR.points)
        phi[:,1] = 1 / 4 * (1 - ctes_FR.points)

        Ft_face_phi = (Ft_internal[:,:,np.newaxis,np.newaxis] * phi[np.newaxis,np.newaxis,:])

        'Look for a faster way to do that'
        Ft_SP_reshaped = np.empty((1,ctes.n_volumes,ctes_FR.n_points))
        contours = np.array([0,ctes_FR.n_points-1])
        for i in range(ctes_FR.n_points):
            lines = np.array([np.zeros_like(ctes_FR.v0[:,0]), np.zeros_like(ctes_FR.v0[:,1])]).astype(int).flatten()
            cols = np.array([ctes_FR.v0[:,0], ctes_FR.v0[:,1]]).flatten()
            data = np.array([Ft_face_phi[:,:,i,0], Ft_face_phi[:,:,i,1]]).flatten()
            Ft_SP_reshaped[:,:,i] = sp.csc_matrix((data, (lines, cols)), shape = (1, ctes.n_volumes)).toarray()
        Ft_SP = 2 * Ft_SP_reshaped #np.concatenate(np.dsplit(Ft_SP_reshaped, ctes_FR.n_points), axis = 2)

        #Ft_SP[0,wells['all_wells'],:] = ((Ft_internal[0,ctes_FR.vols_vec][wells['all_wells']]).sum(axis=-1)/2)[:,np.newaxis]
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
        P_SP[:,:] = Pold[:,np.newaxis]
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

    def Riemann_Solver(self, M, fprop, wells, Nk_SP, Fk_SP, Ft_internal):
        Nk_faces = np.empty((ctes.n_components, ctes.n_internal_faces, 2))
        Nk_faces[:,:,1] = Nk_SP[:,ctes_FR.v0[:,1],0] #Nk faces a esquerda dos volumes
        Nk_faces[:,:,0] = Nk_SP[:,ctes_FR.v0[:,0],-1] #Nk nas faces a direita

        P_face = np.concatenate((self.P_faces[:,np.newaxis],self.P_faces[:,np.newaxis]),axis=1)

        Fk_faces = np.empty_like(Nk_faces)
        Fk_faces[:,:,1] = Fk_SP[:,ctes_FR.v0[:,1],0]
        Fk_faces[:,:,0] = Fk_SP[:,ctes_FR.v0[:,0],-1]

        #FOR THE BURGERS PROBLEM AND BASTIAN
        '''Nk_face_contour = np.empty((ctes.n_components,1,2))
        Nk_face_contour[:,0,1] = Nk_SP[:,0,0]
        Nk_face_contour[:,0,0] = Nk_SP[:,-1,-1]
        #Nk_face_contour[-1,0,0] = 1*fprop.Csi_j[0,-1,0]*fprop.Vp[0] #BASTIAN

        Vp_face_contour = fprop.Vp[:2]
        RS_contour = RiemannSolvers(np.array([0,ctes.n_volumes-1])[np.newaxis], np.array([ctes.pretransmissibility_internal_faces[0]]))
        Fk_faces_contour = RS_contour.get_Fk_face(fprop, M, Nk_face_contour, \
            P_face[np.newaxis,0], Vp_face_contour, Ft_internal[:,0][:,np.newaxis])
        Fk_face_contour_RS, alpha_wv =  RS_contour.LLF(M, fprop, Nk_face_contour, P_face[np.newaxis,0],
            Ft_internal[:,0][:,np.newaxis], Fk_faces_contour, np.ones(1,dtype=bool))

        RS = RiemannSolvers(ctes_FR.v0, ctes.pretransmissibility_internal_faces)
        Fk_face_RS, alpha_wv =  RS.LLF(M, fprop, Nk_faces, P_face, Ft_internal, Fk_faces, \
            np.ones_like(Ft_internal[0],dtype=bool))'''
        #ponteiro = np.zeros_like(Ft_internal[0], dtype=bool)
        #ponteiro[[0,1]] = True
        #Fk_face_RS[:,ponteiro] = MUSCL().update_flux_upwind(fprop.P[np.newaxis,:], Fk_faces[:,ponteiro], ponteiro)

        Fk_face_RS[:,0] = Fk_faces[:,0,0]
        Fk_face_RS[:,1] = Fk_faces[:,1,0]

        #Fk_face_RS[:,[0,1]] = Fk_face_contour_RS

        'Obtaining Flux at each CV side - by finding faces that compounds the CV \
        this only works for 1D problems'
        self.vols_vec = -np.ones((ctes.n_volumes,2),dtype=int)
        lines = np.arange(ctes.n_internal_faces)
        self.vols_vec[ctes_FR.v0[:,0],1] = lines
        self.vols_vec[ctes_FR.v0[:,1],0] = lines
        Fk_vols_RS_neig = Fk_face_RS[:,ctes_FR.vols_vec]
        Fk_vols_RS_neig[:,self.vols_vec<0] = 0 #Fk_face_contour_RS #FOR THE BURGERS

        #For the bastian
        #Fk_vols_RS_neig[-1,-1] = 0
        #Fk_vols_RS_neig[0,-1,-1] = 2 * Fk_vols_RS_neig[0,-1,-1]
        return Fk_faces, Fk_vols_RS_neig, alpha_wv

    def dFlux_Continuous(self, Fk_SP, Fk_vols_RS_neig):
        Fk_D_l = Fk_SP[...,0]#Fk_D @ x_left
        Fk_D_r = Fk_SP[...,-1]#Fk_D @ x_right
        dFk_C = (Fk_vols_RS_neig[:,:,0] - Fk_D_l)[:,:,np.newaxis] * ctes_FR.dgLB[np.newaxis,:] + \
                (Fk_vols_RS_neig[:,:,1] - Fk_D_r)[:,:,np.newaxis] * ctes_FR.dgRB[np.newaxis,:]
        return dFk_C


    def minmod(self, dNks):
        s = (np.sign(dNks[...,0]) + np.sign(dNks[...,1]) + np.sign(dNks[...,2]))/3
        mmod = np.zeros_like(dNks[...,0])
        ids2 = abs(s)==1
        mmod[ids2] = (s * np.min(abs(dNks),axis=-1))[ids2]
        return mmod

    def minmodB(self, us, M, h):
        #TVB based minmod function
        mfunc = np.copy(us[...,0])
        #par = (M*h**2)
        par = 1e-330
        ids = (np.abs(mfunc) > par)
        mfunc[ids] = self.minmod(us[ids,:])
        return mfunc

    def slopeLimN(self, M, Nk_SP_in):
        #Compute cell averages
        Nk = (np.linalg.inv(ctes_FR.V)[np.newaxis,] @ Nk_SP_in[:,:,:, np.newaxis])[:,:,:,0]
        Nk_SP = self.projections(Nk, ctes_FR.n_points-1)

        'Projected n=0'
        Nk_avg = self.projections(Nk, 0)
        Nk_avg0 = Nk_avg[...,0]

        eps0 = 1e-330
        Nk_L = Nk_SP[...,0]
        Nk_R = Nk_SP[...,1]

        #only works for 1D
        Nk_avg_neig_1 = np.concatenate((Nk_avg0[:,0][:,np.newaxis], Nk_avg0[:,0:-1]),axis=1)
        Nk_avg_neig_2 = np.concatenate((Nk_avg0[:,1:], Nk_avg0[:,-1][:,np.newaxis]),axis=1)

        dNkL1 = (Nk_avg0 - Nk_L)
        dNkR1 = (Nk_R - Nk_avg0)
        dNk2 = (Nk_avg0 - Nk_avg_neig_1)
        dNk3 = (Nk_avg_neig_2 - Nk_avg0)

        dNk = np.concatenate((dNk2[...,np.newaxis], dNk3[...,np.newaxis]),axis=-1)
        dNkL = np.concatenate((dNk, dNkL1[...,np.newaxis]),axis=-1)
        dNkR = np.concatenate((dNk, dNkR1[...,np.newaxis]),axis=-1)
        Nk_lim_L = Nk_avg0 - self.minmod(dNkL)
        Nk_lim_R = Nk_avg0 + self.minmod(dNkR)

        ponteiro = np.zeros((ctes.n_components,ctes.n_volumes),dtype=bool)
        dL = abs((Nk_L-Nk_lim_L))#/Nk_L)
        dR = abs((Nk_R-Nk_lim_R))#/Nk_R)
        ponteiro[(dL>eps0) + (dR>eps0)] = True

        #Limit function
        c_vols = M.data['centroid_volumes'][:,0]
        xSP1 = c_vols - c_vols[0]
        c_faces = c_vols[ctes.v0]
        h = c_vols[0]*2

        x0 = np.ones((ctes_FR.n_points,ctes.n_volumes)) * (xSP1 + h/2)[np.newaxis,:]
        x_SPs = x0 + (ctes_FR.points * h/ctes.n_points)[:,np.newaxis]

        Nk_P1 = self.projections(Nk, 1)

        ux = (2/h)*(ctes_FR.Dr@Nk_P1[ponteiro,:].T)
        ux = ux.T
        ux1 = (Nk_avg0-Nk_avg_neig_1)/h
        ux2 = (Nk_avg_neig_2-Nk_avg0)/h
        us = np.concatenate((ux[...,0][...,np.newaxis],ux2[ponteiro,np.newaxis]),axis=-1)
        us =  np.concatenate((us,ux1[ponteiro,np.newaxis]),axis=-1)
        Nk_lim = np.copy(Nk_SP)
        Nk_lim[ponteiro] = Nk_avg[ponteiro] + ((x_SPs - x0).T[np.newaxis,:] * \
            np.ones_like(Nk_SP))[ponteiro]* (self.minmodB(us, 0.2, h)[...,np.newaxis])
        if any((Nk_lim<0).flatten()): import pdb; pdb.set_trace()
        return Nk_lim

    def MLP_slope_limiter(self, M, fprop, Nk_SP_in, wells):

        inds = np.array([0,-1])


        Nk = (np.linalg.inv(ctes_FR.V)[np.newaxis,] @ Nk_SP_in[:,:,:, np.newaxis])[:,:,:,0]
        Nk_SP = self.projections(Nk, ctes_FR.n_points-1)

        #self.machine_error = np.max(abs(Nk_SP - Nk_SP_in))
        #if self.machine_error<np.finfo(np.float64).eps: machine_error = np.finfo(np.float64).eps
        self.machine_error = 0 #1e-150 #1e-323 #1e-700 #1e-150

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

        Phi_P1 = self.P1_limiter(M, Nk, Nk_SP, Nk_P0_vertex, Nk_P1_vertex, Nk_neig, Nk_faces)

        phi_P2 = np.zeros_like(Phi_P1)
        phi_Pn = np.zeros_like(Phi_P1)
        Nk_P2 = Nk_SP

        #phi_Pn[:,[0,-1]]=1

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

        #tornando nao hierarquico
        #phi_P2[:] = 1
        #phi_Pn[:] = 1

        #Phi_P1[:] = 1
        Nk_SPlim = (Nk_P0 + Phi_P1[:,:, np.newaxis] * (Nk_P1 - Nk_P0) + \
            phi_P2[:,:,np.newaxis] * ((Nk_P2 - Nk_P1) + phi_Pn[:,:,np.newaxis] * (Nk_SP - Nk_P2)))

        Nk = 1 / sum(ctes_FR.weights) * np.sum(ctes_FR.weights * Nk_SPlim,axis=2)

        z = Nk[0:ctes.Nc,:] / np.sum(Nk[0:ctes.Nc,:], axis = 0)
        #if any(abs(z[1,90:]-0.3)>1e-4): import pdb; pdb.set_trace()

        '''if any(abs(Nk_SPlim[Nk_SPlim<0])):
            import pdb; pdb.set_trace()'''

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

        Phi_r = self.MLP_u1_mod(Nk_neig, Nk_P0_vertex, Linear_term)
        #Phi_r = self.MLP_u1(Nk_neig, Nk_P0_vertex, Linear_term)
        #Phi_r = self.MLP_vk(M, Nk_neig, Nk_P0_vertex, Linear_term)

        Phi_P1 = np.ones_like(Phi_r)
        Phi_P1[abs(Nk_P1_vertex - Nk_P0_vertex) >= self.machine_error] = Phi_r[abs(Nk_P1_vertex - Nk_P0_vertex) >= self.machine_error]
        Phi_P1 = np.min(Phi_P1, axis = 2)

        #FOR THE BURGERS PROBLEM
        #Phi_P1[:,-1] = 1
        #Phi_P1[:,0] = 1
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
        K1 = 0 #1e-10 #BL use K1=2 or K1=6 for P3
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

        #Phi_r_u2[Phi_r_u2 > 1] = 1
        #Phi_r_u2[Phi_r_u2 < 0] = 0

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

    def smooth_extrema_cond(self, fprop, Nk_P0_vertex, Nk_P1_vertex, Nk_Pm_vertex, Nk_neig):
        'Smooth extrema detector'
        High_order_term = Nk_Pm_vertex - Nk_P1_vertex
        Linear_term = Nk_P1_vertex - Nk_P0_vertex

        Nk_less_max = (Nk_Pm_vertex < np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        Nk_big_min = (Nk_Pm_vertex > np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec])
        #Nk_less_max[abs(Nk_Pm_vertex - np.max(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]) <= self.machine_error] = True
        #Nk_big_min[abs(Nk_Pm_vertex - np.min(Nk_neig, axis = 2)[:,ctes_FR.vols_vec]) <= self.machine_error] = True
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

        C1, C2 = self.smooth_extrema_cond(fprop, Nk_P0_vertex, Nk_P1_vertex, Nk_Pmvertex, Nk_neig)

        ''' Dehativation inequality for nearly constant region'''
        e = 1e-3
        C3 = (abs(Nk_Pmvertex - Nk_P0_vertex)/Nk_P0_vertex <= e)
        #import pdb; pdb.set_trace()
        Nf = self.Gibbs_Wilbraham_oscilations(M, Nk_faces)

        '''A3 step'''
        #P1_proj_cond[C3] = True
        C3_cells = np.min(1*(C3),axis=-1)

        '''A4-1 step'''
        P1_proj_cond_cells = np.min(1*(P1_proj_cond),axis=-1)

        'Type I trouble'
        condI = (P1_proj_cond_cells==1)
        aux_I = P1_proj_cond_cells[condI]
        aux_Nf_I = Nf[condI]
        aux_I[aux_Nf_I>=1] = 0
        P1_proj_cond_cells[condI] = aux_I

        '''A4-2 step'''
        smooth_extrema = C1+C2
        smooth_extrema_cells = np.min(1*smooth_extrema,axis=-1)

        'Type II trouble'
        #shock_condition = abs(p_neig_l-p_neig_r)>=0.1*p_neig_l#
        aux_II = smooth_extrema_cells[smooth_extrema_cells==1]
        aux_Nf_II = Nf[smooth_extrema_cells==1]
        aux_II[aux_Nf_II>=1] = 0
        smooth_extrema_cells[smooth_extrema_cells==1] = aux_II

        phi_Pm = 1*((P1_proj_cond_cells + smooth_extrema_cells).astype(bool))
        #phi_Pm = 1*((P1_proj_cond_cells).astype(bool))

        if any((phi_Pm[:,:,np.newaxis]*Nk_Pmvertex).flatten()<0): import pdb; pdb.set_trace()

        if any(phi_Pm.flatten()>1): import pdb; pdb.set_trace()
        return phi_Pm
