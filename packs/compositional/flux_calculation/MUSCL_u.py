from .riemann_solvers import RiemannSolvers
from packs.utils import constants as ctes
from packs.utils import QR_decomp
from packs.directories import data_loaded
import scipy.sparse as sp
import numpy as np
from packs.compositional import prep_MUSCL as ctes_MUSCL
from .flux_volumes import Flux
from .MUSCL import MUSCL
import time

class MUSCL_u(MUSCL):

    """ Class created for the second order MUSCL implementation for the \
    calculation of the advective terms """

    def run(self, M, fprop, wells, P_old, Nk_old, ftot, Pot_hid):

        ''' Global function that calls others '''
        self.Nk = np.copy(Nk_old)
        self.P_face = np.sum(P_old[ctes.v0], axis=1) * 0.5
        self.P_face = np.concatenate((self.P_face[:,np.newaxis], self.P_face[:,np.newaxis]),axis=1)
        dNk_vols = self.volume_gradient_reconstruction(M, fprop, wells)
        #self.P_face[ctes_MUSCL.faces_contour] = P_old[ctes.v0[ctes_MUSCL.faces_contour,0]]
        dNk_face, dNk_face_neig = self.get_faces_gradient(M, fprop, dNk_vols)

        lim = getattr(self, ctes.MUSCL['lim'])
        #Phi = self.Van_Leer_slope_limiter(dNk_face, dNk_face_neig)
        #Phi = self.Van_Albada1_slope_limiter(dNk_face, dNk_face_neig)
        Phi = lim(dNk_face, dNk_face_neig)

        #BURGERS
        #Phi[:] = 1
        #Phi[:,:,1] = -1

        Nk_face, z_face = self.get_extrapolated_compositions(fprop, Phi, dNk_face_neig)
        
        #G = self.update_gravity_term() # for now, it has no gravity
        alpha, Fk_vols_total = self.update_flux(M, wells, fprop, Nk_face,
            ftot, Pot_hid)
        #alpha = fprop.Fk_vols_total/self.Nk
        #import pdb; pdb.set_trace()

        return alpha, Fk_vols_total
    
    
    def volume_gradient_reconstruction(self, M, fprop, wells):
        Nk_neig =  self.Nk[:,np.newaxis,:] * ctes_MUSCL.allneig_and_vol[np.newaxis,:,:]
        Nk = Nk_neig.transpose(0,2,1)

        dNk = Nk_neig - Nk #b vector
        #import pdb; pdb.set_trace()
        
        x_neig = M.volumes.center(M.volumes.all).T[:,np.newaxis,:] * ctes_MUSCL.allneig_and_vol[np.newaxis,:]
        x = x_neig.transpose(0,2,1)
        dx = x_neig - x
        mod_dx = np.sqrt(np.sum(dx*dx,axis=0))
        a = 1
        w_pond = 1/(mod_dx**a) #ponderation term forming Y matrix
        w_pond[np.isinf(w_pond)] = 0
        U1 = dx[0:2]
        #x_bar = x_face_center - x_neig
        #x_hat2 = x_bar**2 + 2*x_bar*(dx) + dx**2

        dNkdx = np.zeros((ctes.n_components,ctes.n_volumes, ctes.n_volumes))
        dNkdy = np.zeros((ctes.n_components,ctes.n_volumes, ctes.n_volumes))
        'think how to vectorize this... its going to take a while...'
        for i in range(ctes.n_volumes):
            Y_mtrx = np.diag(w_pond[i][w_pond[i]>0])
            U = U1[:,i]
            U_vec = U[:,w_pond[i]>0].T
            dNk_v0 = dNk[:,i]
            b_vec = dNk_v0[:,w_pond[i]>0]
            a = Y_mtrx @ b_vec.T
            
            'Orthogonal decomposition'
            A_x = Y_mtrx*U_vec[:,0]
            A_y = Y_mtrx*U_vec[:,1]

            #t0_code = time.time()
            A_x_M = QR_decomp.Matrix(n_row=len(A_x[:,0]), n_col=len(A_x[0,:]), two_d_array=A_x.tolist())
            A_y_M = QR_decomp.Matrix(n_row=len(A_y[:,0]), n_col=len(A_y[0,:]), two_d_array=A_y.tolist())
            
            Q_x, R_x = QR_decomp.QR_GS(A_x_M)
            Q_y, R_y = QR_decomp.QR_GS(A_y_M)
            R_x = np.array(R_x.values); Q_x = np.array(Q_x.values)
            R_y = np.array(R_y.values); Q_y = np.array(Q_y.values)
            #t1_code = time.time()
            #dt_code = t1_code - t0_code
            #print('dt_code: ', dt_code)
            #import pdb; pdb.set_trace()

            'Solution'
            dNkdx[:,i,w_pond[i]>0] = ((np.linalg.inv(R_x)@(Q_x.T))@a).T # para os vizinhos de 1 CV
            dNkdy[:,i,w_pond[i]>0] = ((np.linalg.inv(R_y)@(Q_y.T))@a).T # para os vizinhos de 1 CV

        
        #zero in the contour volumes
        dNkdx[:,(w_pond>0).sum(axis=1)<4] = 0 #dNkds_vols[:,ctes_MUSCL.all_neig_by_axes==1]
        dNkdy[:,(w_pond>0).sum(axis=1)<4] = 0
        #dNkds_vols[:,:,ctes.v0[ctes_MUSCL.faces_contour].flatten()] = 0 # zero in the contour volumes
        dNk_vols = dNkdx * dx[0,:] + dNkdy * dx[1,:]
        dNk_vols[np.isnan(dNk_vols)] = 0
        return dNk_vols

    def get_faces_gradient(self, M, fprop, dNk_vols):
        #ds_face = (M.data['centroid_volumes'][ctes.v0[:,1],:] -  M.data['centroid_volumes'][ctes.v0[:,0],:])
        #ds_face_abs = abs(ds_face)
        dNk_face = (self.Nk[:,ctes.v0[:,1]] - self.Nk[:,ctes.v0[:,0]]) #* ctes_MUSCL.versor_ds_face[np.newaxis,:]
        dNk_face_vols = 2. * (dNk_vols[:,:,ctes.v0]).sum(axis=1)
        dNk_face_neig = dNk_face_vols - dNk_face[:,:,np.newaxis]
        #dNk_face2 = dNk_face[...,np.newaxis]  * np.ones_like(dNk_face_neig)
        #dNk_face_vols[abs(dNk_face_neig)<1e-25] = dNk_face2[abs(dNk_face_neig)<1e-25]
        return dNk_face, dNk_face_neig

    def Van_Leer(self, dNk_face, dNk_face_neig):
        np.seterr(divide='ignore', invalid='ignore')
        r_face = dNk_face[:,:,np.newaxis] / dNk_face_neig
        r_face[dNk_face_neig==0] = 0
        phi = (r_face + abs(r_face)) / (r_face + 1)
        phi[r_face<0] = 0 #so botei pra caso r==-1
        Phi = phi
        Phi[:,:,1] = -Phi[:,:,1]
        return Phi

    def Van_Albada1(self, dNk_face, dNk_face_neig):
        np.seterr(divide='ignore', invalid='ignore')
        r_face = dNk_face[:,:,np.newaxis] / dNk_face_neig
        r_face[dNk_face_neig==0] = 0
        phi = (r_face**2 + (r_face)) / (r_face**2 + 1)
        phi[r_face<0]=0 #so botei pra caso r==-1
        Phi = phi
        Phi[:,:,1] = -Phi[:,:,1]
        return Phi

    def minmod(self, dNk_face, dNk_face_neig):
        np.seterr(divide='ignore', invalid='ignore')
        r_face = dNk_face[:,:,np.newaxis] / dNk_face_neig
        r_face[dNk_face_neig==0] = 0
        phi = np.copy(r_face)
        phi[r_face>1]=1
        phi[r_face<0]=0 #so botei pra caso r==-1
        Phi = phi
        Phi[:,:,1] = -Phi[:,:,1]
        return Phi

    def get_extrapolated_compositions(self, fprop, Phi, dNk_face_neig):
        Phi[abs(dNk_face_neig)<1e-30] = 0
        Nk_face = self.Nk[:,ctes.v0] + Phi / 2 * dNk_face_neig
        #Nk_face[(Nk_face<0)*(abs(Nk_face)<1e-30)] = 0
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

    def update_flux_upwind(self, Pot_hid, Fk_face_upwind_all, v0, ponteiro):
        Fk_face_upwind = np.empty_like(Fk_face_upwind_all[:,:,0])

        Pot_hidj = Pot_hid[0,v0[:,0]][ponteiro] #- G[0,:,:,0]
        Pot_hidj_up = Pot_hid[0,v0[:,1]][ponteiro] #- G[0,:,:,1]

        Fk_face_upwind[:,Pot_hidj_up <= Pot_hidj] = \
            Fk_face_upwind_all[:,Pot_hidj_up <= Pot_hidj, 0]
        Fk_face_upwind[:,Pot_hidj_up > Pot_hidj] = \
            Fk_face_upwind_all[:,Pot_hidj_up > Pot_hidj, 1]

        return Fk_face_upwind

    def update_flux(self, M, wells, fprop, Nk_face, ftotal, Pot_hid):
        Fk_internal_faces = np.zeros((ctes.n_components,ctes.n_internal_faces))
        RS = RiemannSolvers(ctes.v0, ctes.pretransmissibility_internal_faces)
        Vp_face = np.concatenate((fprop.Vp[ctes.v0[:,0]], fprop.Vp[ctes.v0[:,1]]), axis=0)
        Fk_face = RS.get_Fk_face(fprop, M, Nk_face, self.P_face, Vp_face, ftotal)
        #Fk_face = (Nk_face**2/2)/(1/ctes.n_volumes) #burgers

        ponteiro = np.zeros(ctes.n_internal_faces,dtype=bool)
        'Corrigir o cálculo do fluxo por upw - aplicar antese Fk_face vem um só'
        solver = getattr(RS, ctes.RS)
        #alpha_wv = np.zeros((ctes.n_internal_faces, 5))

        Fk_internal_faces[:,~ponteiro], alpha_wv = solver(M, \
                fprop, Nk_face, self.P_face, ftotal, Fk_face, ~ponteiro)


        #FOR THE BURGERS PROBLEM AND BASTIAN
        '''Nk_face_contour = np.empty((ctes.n_components,1,2))
        dNk_face_0 = (self.Nk[:,1] - self.Nk[:,0])
        dNk_face_end = (self.Nk[:,-1] - self.Nk[:,-2])
        Nk_face_contour[:,0,1] = self.Nk[:,0] - 1/2*dNk_face_0
        Nk_face_contour[:,0,0] = self.Nk[:,-1] + 1/2*dNk_face_end
        P_face = self.P_face
        #Nk_face_contour[-1,0,0] = 1*fprop.Csi_j[0,-1,0]*fprop.Vp[0] #BASTIAN
        Vp_face_contour = fprop.Vp[:2]
        RS_contour = RiemannSolvers(np.array([ctes.n_volumes-1,0])[np.newaxis], np.array([ctes.pretransmissibility_internal_faces[0]]))
        Fk_faces_contour = (Nk_face_contour**2)/2/(1/ctes.n_volumes)
        Fk_face_contour_RS, alpha_wv2 =  RS_contour.LLF(M, fprop, Nk_face_contour, P_face[np.newaxis,0],
            ftotal[:,0][:,np.newaxis], Fk_faces_contour, np.ones(1,dtype=bool))
        '''
        ponteiro[ctes_MUSCL.faces_contour] = True
        Fk_internal_faces[:,ponteiro] = self.update_flux_upwind(fprop.P[np.newaxis,:], \
            Fk_face[:,ponteiro], ctes.v0, ponteiro)

        #Fk_internal_faces[:,0] = Fk_face[:,0,0] #comment for burgers
        #Fk_internal_faces[:,1] = Fk_face[:,1,0] #comment for burgers
        '-------- Perform volume balance to obtain flux through volumes -------'
        Fk_vols_total = Flux().update_flux_volumes(Fk_internal_faces)

        #BURGERS
        #Fk_vols_total[:,0] += Fk_face_contour_RS[0]
        #Fk_vols_total[:,-1] -= Fk_face_contour_RS[0]

        if any(np.isnan(Fk_vols_total).flatten()): import pdb; pdb.set_trace()
        #Fk_vols_total[:ctes.Nc][fprop.z==0] = 0
        if any(Fk_vols_total[:ctes.Nc][fprop.z==0]<0): import pdb; pdb.set_trace()
        return alpha_wv, Fk_vols_total
