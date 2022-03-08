import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from scipy import linalg
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve


class NewtonSolver:
    def __init__(self, M, wells, fprop, delta_t, Pot_hid, Nk_old, dVjdNk, dVjdP):
        self.M = M
        self.wells = wells
        self.fprop = fprop
        self.delta_t = delta_t
        self.Pot_hid = Pot_hid
        self.Nk_old = Nk_old
        self.dVt_derivatives(dVjdNk, dVjdP)

    def dVt_derivatives(self, dVjdNk, dVjdP):
        self.dVtP = dVjdP.sum(axis=1)[0]
        self.dVtk = dVjdNk.sum(axis=1)

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

        '''
        d_Pot_hid = np.zeros_like(self.fprop.mobilities_internal_faces)
        for i in range(3):
            for j in range(len(d_Pot_hid[0][0])):
                d_Pot_hid[0][i][j] = self.Pot_hid[i][j+1] - self.Pot_hid[i][j]
        '''
        Pot_hidj_test = self.Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up_test = self.Pot_hid[:,ctes.v0[:,1]]
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        d_Pot_hid = Pot_hidj_up_test - Pot_hidj_test - ctes.g * self.fprop.rho_j_internal_faces * (z_up - z)

        flux = transmissibility * d_Pot_hid
        Fk_internal_faces = flux.sum(axis = 1)

        '''
        fluxos_contornos_dir = np.zeros([ctes.Nc+1,1])
        fluxos_contornos_esq = np.zeros([ctes.Nc+1,1])
        fluxos_totais = np.insert(fluxos_internos, [ctes.Nc+1], fluxos_contornos_dir, axis = 1)
        fluxos_totais = np.insert(fluxos_totais, [0], fluxos_contornos_esq, axis = 1)
        d_fluxos_totais = np.zeros([ctes.n_components, ctes.n_volumes])
        for i in range(len(fluxos_totais)):
            d_fluxos_totais[i][:] = fluxos_totais[i][1:] - fluxos_totais[i][:-1]
        '''

        Fk_vols_total = self.component_flux_volumes(Fk_internal_faces)
        #Eduarda considera fluxo saindo negativo e entrando positivo?

        ''' Termo de poço ta certo ? falta o termo do produtor '''
        well_term = self.well_term(self.fprop, self.wells) # Poço injetor apenas ?
        import pdb; pdb.set_trace()
        residuo_massa = (self.fprop.Nk - self.Nk_old) - self.delta_t * d_fluxos_totais - self.delta_t * self.q

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

    def component_flux_volumes(self, Fk_internal_faces):
        ''' Function to compute component molar flux balance through the control \
        volume interfaces '''
        cx = np.arange(ctes.n_components)
        lines = np.array([np.repeat(cx,len(ctes.v0[:,0])), np.repeat(cx,len(ctes.v0[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(ctes.v0[:,0],ctes.n_components), np.tile(ctes.v0[:,1], ctes.n_components)]).flatten()
        data = np.array([-Fk_internal_faces, Fk_internal_faces]).flatten()
        Fk_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_components, ctes.n_volumes)).toarray()

        return Fk_vols_total
