import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from scipy import linalg
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

class TPFASolver:
    def __init__(self, dVjdNk, dVjdP):
        self.dVt_derivatives(dVjdNk, dVjdP)

    def get_pressure(self, M, wells, fprop, Pold, delta_t):
        T = self.update_transmissibility(M, wells, fprop, delta_t)
        D = self.update_independent_terms(M, fprop, Pold, wells, delta_t)
        Pnew = self.update_pressure(T, D)
        Ft_internal_faces = self.update_total_flux_internal_faces(M, fprop, Pnew)
        self.update_flux_wells(fprop, Pnew, wells, delta_t)
        return Pnew, Ft_internal_faces, self.q

    def dVt_derivatives(self, dVjdNk, dVjdP):
        self.dVtP = dVjdP.sum(axis=1)[0]
        self.dVtk = dVjdNk.sum(axis=1)

    def update_transmissibility(self, M, wells, fprop, delta_t):
        self.t0_internal_faces_prod = fprop.xkj_internal_faces * \
                                      fprop.Csi_j_internal_faces * \
                                      fprop.mobilities_internal_faces

        ''' Transmissibility '''
        t0 = (self.t0_internal_faces_prod).sum(axis = 1)
        t0 = t0 * ctes.pretransmissibility_internal_faces
        T = sp.csr_matrix((ctes.n_volumes, ctes.n_volumes))
        # Look for a way of doing this not using a loop!!!
        #import pdb; pdb.set_trace()

        for i in range(ctes.n_components):
            lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
            cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
            data = np.array([-t0[i,:], -t0[i,:], +t0[i,:], +t0[i,:]]).flatten()

            Ta = (sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes)))#.toarray()
            T += Ta.multiply(self.dVtk[i, :, np.newaxis])

        T *= delta_t
        ''' Transmissibility diagonal term '''
        #diag = np.diag((ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP))
        #T += diag
        T.setdiag(T.diagonal().flatten() + (ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP))
        self.T_noCC = np.copy(T.toarray())#    .toarray() #Transmissibility without contour conditions

        ''' Includding contour conditions '''
        T[wells['ws_p'],:] = 0
        T[wells['ws_p'], wells['ws_p']] = 1
        return T

    def pressure_independent_term(self, fprop, Pold):
        vector = ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP
        pressure_term = vector * Pold
        return pressure_term

    def capillary_and_gravity_independent_term(self, fprop):

        t0_j = self.t0_internal_faces_prod * \
            ctes.pretransmissibility_internal_faces
        t0_k = ctes.g * np.sum(fprop.rho_j_internal_faces * t0_j, axis=1)

        # Look for a better way to do this !!!
        cap = np.zeros([ctes.n_volumes])
        grav = np.zeros([ctes.n_volumes,ctes.n_volumes])
        if any((ctes.z - ctes.z[0]) != 0):
            for i in range(ctes.n_components):
                lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                data = np.array([t0_k[i,:], t0_k[i,:], -t0_k[i,:], -t0_k[i,:]]).flatten()
                t0_rho = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes)).toarray()
                grav += t0_rho * self.dVtk[i,:]

                for j in range(ctes.n_phases):
                    lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                    cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
                    data = np.array([t0_j[i,j,:], t0_j[i,j,:], -t0_j[i,j,:], -t0_j[i,j,:]]).flatten()
                    t0 = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes))*self.dVtk[i,:]
                    cap += t0 @ fprop.Pcap[j,:]

        gravity_term = grav @ ctes.z

        # capillary_term = np.sum(self.dVtk * np.sum (fprop.xkj *
        #         fprop.Csi_j * fprop.mobilities * fprop.Pcap, axis = 1), axis = 0)
        return cap, gravity_term

    def volume_discrepancy_independent_term(self, fprop):
        volume_discrepancy_term = fprop.Vp - fprop.Vt
        if np.max(abs(volume_discrepancy_term)) > 5e-4:
            #import pdb; pdb.set_trace()
            print('hit: ', np.max(abs(volume_discrepancy_term)))
        return volume_discrepancy_term

    def well_term(self, fprop, wells):
        self.q = np.zeros([ctes.n_components, ctes.n_volumes])
        well_term = np.zeros(ctes.n_volumes)
        if len(wells['ws_q']) > 0:
            self.q[:,wells['ws_q']] =  wells['values_q']
            well_term[wells['ws_q']] = np.sum(self.dVtk[:,wells['ws_q']] *
                self.q[:,wells['ws_q']], axis = 0)
        return well_term

    def update_independent_terms(self, M, fprop, Pold, wells, delta_t):
        self.pressure_term = self.pressure_independent_term(fprop, Pold)
        self.capillary_term, self.gravity_term = self.capillary_and_gravity_independent_term(fprop)
        self.volume_term = self.volume_discrepancy_independent_term(fprop)
        well_term = self.well_term(fprop, wells)
        independent_terms = self.pressure_term - self.volume_term  + delta_t * \
            well_term - delta_t * (self.capillary_term + self.gravity_term)
        independent_terms[wells['ws_p']] = wells['values_p'] + ctes.g * \
            fprop.rho_j[0,0,wells['ws_p']] * (ctes.z[wells['ws_p']] - ctes.z[ctes.bhp_ind])
        return independent_terms

    def update_pressure(self, T, D):
        #P = linalg.solve(T,D)
        P = spsolve(T,D)
        return P

    def update_total_flux_internal_faces(self, M, fprop, Pnew):
        Pot_hid = Pnew + fprop.Pcap
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        Ft_internal_faces = - np.sum(fprop.mobilities_internal_faces
            * ctes.pretransmissibility_internal_faces * ((Pot_hidj_up - Pot_hidj) -
            ctes.g * fprop.rho_j_internal_faces * (z_up - z)), axis = 1)
        return Ft_internal_faces

    def update_flux_wells(self, fprop, Pnew, wells, delta_t):
        wp = wells['ws_p']

        if len(wp)>=1:
            if Pnew[2]>Pnew[1]: import pdb; pdb.set_trace()
            well_term =  (self.T_noCC[wp,:] @ Pnew - self.pressure_term[wp] +
                self.volume_term[wp]) / delta_t  + self.capillary_term[wp] + \
                self.gravity_term[wp]
            mob_ratio = fprop.mobilities[:,:,wp] / np.sum(fprop.mobilities[:,:,wp], axis = 1)
            self.q[:,wp] = np.sum(fprop.xkj[:,:,wp] * mob_ratio * fprop.Csi_j[:,:,wp] * well_term, axis = 1)
            fprop.q_phase = mob_ratio * well_term


# -------------------------------------------------------------------------------------

class NewtonSolver:
    def __init__(self, M, wells, fprop, delta_t, Pot_hid, Nk_old):
        self.M = M
        self.wells = wells
        self.fprop = fprop
        self.delta_t = delta_t
        self.Pot_hid = Pot_hid
        self.Nk_old = Nk_old


    def update_transmissibility_FI(self):
        self.t0_internal_faces_prod = self.fprop.xkj_internal_faces * \
                                      self.fprop.Csi_j_internal_faces * \
                                      self.fprop.mobilities_internal_faces

        ''' Transmissibility '''
        #t0 = (self.t0_internal_faces_prod).sum(axis = 1)
        t0 = self.t0_internal_faces_prod
        t0 = t0 * ctes.pretransmissibility_internal_faces
        return t0


    def well_term(self, fprop, wells): # Preciso para calcular os termos de poço
        self.q = np.zeros([ctes.n_components, ctes.n_volumes])
        well_term = np.zeros(ctes.n_volumes)
        if len(wells['ws_q']) > 0:
            self.q[:,wells['ws_q']] =  wells['values_q']
            well_term[wells['ws_q']] = np.sum(self.dVtk[:,wells['ws_q']] *
                self.q[:,wells['ws_q']], axis = 0)
        return well_term

    def residual_calculation(self):
        # Equacao de conservacao da massa
        transmissibility = self.update_transmissibility_FI()
        d_Pot_hid = np.zeros_like(self.fprop.mobilities_internal_faces)

        '''Ajustar para ficar automatico com as faces'''
        for i in range(3):
            for j in range(len(d_Pot_hid[0][0])):
                d_Pot_hid[0][i][j] = self.Pot_hid[i][j+1] - self.Pot_hid[i][j]

        '''
        # Tentativa de organizar calculo acima
        Pot_hidj = self.Pot_hid[0,ctes.v0[:,0]]
        Pot_hidj_up = self.Pot_hid[0,ctes.v0[:,1]]
        teste = Pot_hidj_up - Pot_hidj
        '''

        flux = transmissibility * d_Pot_hid
        fluxos_internos = flux.sum(axis = 1)

        ''' Errado, so funciona para 1D - AJUSTAR '''
        fluxos_contornos_dir = np.zeros([ctes.Nc+1,1])
        fluxos_contornos_esq = np.zeros([ctes.Nc+1,1])
        fluxos_totais = np.insert(fluxos_internos, [ctes.Nc+1], fluxos_contornos_dir, axis = 1)
        fluxos_totais = np.insert(fluxos_totais, [0], fluxos_contornos_esq, axis = 1)

        d_fluxos_totais = np.zeros([ctes.n_components, ctes.n_volumes])
        for i in range(len(fluxos_totais)):
            d_fluxos_totais[i][:] = fluxos_totais[i][1:] - fluxos_totais[i][:-1]

        ''' Adicionar termo de poço '''
        residuo_massa = (self.fprop.Nk - self.Nk_old) - self.delta_t * d_fluxos_totais

        # Restricao de volume poroso
        aux = self.fprop.Nj / (self.fprop.Csi_j + 1e-15)
        residuo_poroso_varavei = aux.sum(axis = 1)/self.fprop.Vp - 1.0
        residuo_poroso_bruno = aux.sum(axis = 1) - self.fprop.Vp
        ''' Preferencia: Varavei '''

        # Residuo total
        residuo = np.append(residuo_poroso_varavei, residuo_massa, axis = 0)
        import pdb; pdb.set_trace()
        return residuo

    def jacobian_calculation(self):
        # Calculo da matriz jacobiana

        '''
        Funcoes que vou precisar:
        - derivada da eq de restricao de volume poroso:
             - P
             - N
        - derivada da eq de conservacao da massa:
             - P
             - N
        '''

        pass

    def solver(self):
        #Nk_newton = np.copy(self.fprop.Nk)
        #P_Newton = np.copy(self.fprop.P)
        #P_Newton = P_Newton.reshape(1, len(P_Newton))
        #variaveis_primarias = np.append(P_Newton, Nk_newton, axis = 0)

        stop_criteria = 1e-4
        contador = 0
        residuo = self.residual_calculation()
        #import pdb; pdb.set_trace()

        while np.max(residuo) > stop_criteria:
            # calculo dos residuos, jacobiana, e solucao do sistema
            Nk_newton = np.copy(self.fprop.Nk)
            P_Newton = np.copy(self.fprop.P)
            P_Newton = P_Newton.reshape(1, ctes.n_volumes)
            variaveis_primarias = np.append(P_Newton, Nk_newton, axis = 0)

            residuo = self.residual_calculation()
            jacobiana = self.jacobian_calculation()
            import pdb; pdb.set_trace()


            # atualizar variaveis

            contador += 1
            if contador > 100:
                import pdb; pdb.set_trace()


        return residuo
