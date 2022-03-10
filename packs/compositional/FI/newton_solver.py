import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from scipy import linalg
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
#from ..equation_of_state import PengRobinson as PR


class NewtonSolver:
    def __init__(self, M, wells, fprop, delta_t, Pot_hid, Nk_old):
        self.M = M
        self.wells = wells
        self.fprop = fprop
        self.delta_t = delta_t
        self.Pot_hid = Pot_hid
        self.Nk_old = Nk_old
        self.EOS = ctes.EOS_class(fprop.T)

    def solver(self):
        #Nk_newton = np.copy(self.fprop.Nk)
        #P_Newton = np.copy(self.fprop.P)
        #P_Newton = P_Newton.reshape(1, len(P_Newton))
        #variaveis_primarias = np.append(P_Newton, Nk_newton, axis = 0)

        stop_criteria = 1e-4
        contador = 0
        residuo = self.residual_calculation()

        while np.max(residuo) > stop_criteria:
            # calculo dos residuos, jacobiana, e solucao do sistema
            Nk_newton = np.copy(self.fprop.Nk)
            P_Newton = np.copy(self.fprop.P)
            P_Newton = P_Newton.reshape(1, ctes.n_volumes)
            variaveis_primarias = np.append(P_Newton, Nk_newton, axis = 0)
            ''' Rever para usar matriz esparsa '''

            residuo = self.residual_calculation()
            jacobiana = self.jacobian_calculation()
            import pdb; pdb.set_trace()


            # atualizar variaveis

            contador += 1
            if contador > 100:
                import pdb; pdb.set_trace()


        return residuo

    def residual_calculation(self):
        # Equacao de conservacao da massa
        transmissibility = self.transmissibility_FI()

        Pot_hidj_test = self.Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up_test = self.Pot_hid[:,ctes.v0[:,1]]
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        d_Pot_hid = Pot_hidj_up_test - Pot_hidj_test - ctes.g * self.fprop.rho_j_internal_faces * (z_up - z)

        flux = transmissibility * d_Pot_hid
        Fk_internal_faces = flux.sum(axis = 1)
        Fk_vols_total = self.component_flux_volumes(Fk_internal_faces)
        '''Eduarda considera fluxo saindo negativo e entrando positivo?'''

        self.well_term(self.wells)
        residuo_massa = (self.fprop.Nk - self.Nk_old) - self.delta_t * Fk_vols_total - self.delta_t * self.q

        # Restricao de volume poroso
        aux = self.fprop.Nj / (self.fprop.Csi_j + 1e-15)
        residuo_poroso_varavei = aux.sum(axis = 1)/self.fprop.Vp - 1.0
        residuo_poroso_bruno = aux.sum(axis = 1) - self.fprop.Vp
        # Preferencia: Varavei

        # Residuo total
        residuo = np.append(residuo_poroso_varavei, residuo_massa, axis = 0)
        return residuo

    def transmissibility_FI(self):
        self.t0_internal_faces_prod = self.fprop.xkj_internal_faces * \
                                      self.fprop.Csi_j_internal_faces * \
                                      self.fprop.mobilities_internal_faces

        ''' Transmissibility '''
        #t0 = (self.t0_internal_faces_prod).sum(axis = 1)
        t0 = self.t0_internal_faces_prod
        t0 = t0 * ctes.pretransmissibility_internal_faces
        return t0

    def well_term(self, wells): # Preciso para calcular os termos de poço
        self.q = np.zeros([ctes.n_components, ctes.n_volumes])
        if len(wells['ws_q']) > 0:
            self.q[:,wells['ws_q']] =  wells['values_q']

    def component_flux_volumes(self, Fk_internal_faces):
        ''' Function to compute component molar flux balance through the control \
        volume interfaces '''
        cx = np.arange(ctes.n_components)
        lines = np.array([np.repeat(cx,len(ctes.v0[:,0])), np.repeat(cx,len(ctes.v0[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(ctes.v0[:,0],ctes.n_components), np.tile(ctes.v0[:,1], ctes.n_components)]).flatten()
        data = np.array([-Fk_internal_faces, Fk_internal_faces]).flatten()
        Fk_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_components, ctes.n_volumes)).toarray()

        return Fk_vols_total

    def jacobian_calculation(self):
        # Calculo da matriz jacobiana
        d_residuoporoso_dP = self.d_residuoporoso_dP()
        d_residuoporoso_dNk = self.d_residuoporoso_dNk()

        import pdb; pdb.set_trace()
        '''
        Funcoes que vou precisar:
        - derivada da eq de restricao de volume poroso:
             - P > OK
             - N
        - derivada da eq de conservacao da massa:
             - P
             - N
        '''

        pass

    def d_residuoporoso_dP (self):
        x = self.fprop.xkj[0:ctes.Nc,0,:]
        y = self.fprop.xkj[0:ctes.Nc,1,:]
        Nl = self.fprop.Nj[0,0,:]
        Nv = self.fprop.Nj[0,1,:]
        P = self.fprop.P
        T = self.fprop.T

        dlnfildP, dlnfildnij, dZodP_parcial, dZldnij_parcial, Zo = \
                self.EOS.get_phase_derivatives(P, T, x, Nl, np.ones(ctes.n_volumes))
        dlnfivdP, dlnfivdnij, dZvdP_parcial, dZvdnij_parcial, Zv = \
                self.EOS.get_phase_derivatives(P, T, y, Nv, np.zeros(ctes.n_volumes))

        Nj_Vp = self.fprop.Nj/self.fprop.Vp
        dVw_dP = -ctes.Cw * self.fprop.Csi_W0 / (self.fprop.Csi_W ** 2)
        dVo_dP = ctes.R * self.fprop.T * (self.fprop.P * dZodP_parcial - Zo) / (self.fprop.P * self.fprop.P)
        dVv_dP = ctes.R * self.fprop.T * (self.fprop.P * dZvdP_parcial - Zv) / (self.fprop.P * self.fprop.P)

        dRv_dP = Nj_Vp[0,2] * dVw_dP + Nj_Vp[0,0] * dVo_dP + Nj_Vp[0,1] * dVv_dP
        return dRv_dP

    def d_residuoporoso_dNk (self):
        # Essa funcao passará por mudanças para a inclusão do térmico. Ver item C.1 da Tese do Varavei
        
        pass
