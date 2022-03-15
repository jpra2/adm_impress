import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from scipy import linalg
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
#from ..equation_of_state import PengRobinson as PR


class NewtonSolver:
    def __init__(self, fprop):
        #self.M = M
        #self.wells = wells
        #self.fprop = fprop
        #self.delta_t = delta_t
        #self.Pot_hid = Pot_hid
        #self.Nk_old = Nk_old
        self.EOS = ctes.EOS_class(fprop.T)

    def solver(self, wells, fprop, delta_t, Nk_old, G):
        #Nk_newton = np.copy(fprop.Nk)
        #P_Newton = np.copy(fprop.P)
        #P_Newton_aux = P_Newton.reshape(1, len(P_Newton))
        #variaveis_primarias = np.append(P_Newton_aux, Nk_newton, axis = 0)

        Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]
        stop_criteria = 1e-4
        contador = 0
        residuo = self.residual_calculation(fprop, Pot_hid, wells, Nk_old, delta_t)

        while np.max(residuo) > stop_criteria:
            # calculo dos residuos, jacobiana, e solucao do sistema
            #Nk_newton = np.copy(fprop.Nk)
            #P_Newton = np.copy(fprop.P)
            #P_Newton = P_Newton.reshape(1, ctes.n_volumes)
            #variaveis_primarias = np.append(P_Newton, Nk_newton, axis = 0)
            ''' Rever para usar matriz esparsa '''


            residuo = self.residual_calculation(fprop, Pot_hid, wells, Nk_old, delta_t)
            import pdb; pdb.set_trace()
            jacobiana = self.jacobian_calculation(fprop)
            import pdb; pdb.set_trace()


            # atualizar variaveis do fprop
            # G = self.update_gravity_term(fprop)
            # Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]

            contador += 1
            if contador > 100:
                import pdb; pdb.set_trace()


        return residuo

    def residual_calculation(self, fprop, Pot_hid, wells, Nk_old, delta_t):
        # Equacao de conservacao da massa
        transmissibility = self.transmissibility_FI(fprop)

        Pot_hidj_test = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up_test = Pot_hid[:,ctes.v0[:,1]]
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        d_Pot_hid = Pot_hidj_up_test - Pot_hidj_test - ctes.g * fprop.rho_j_internal_faces * (z_up - z)

        flux = transmissibility * d_Pot_hid
        Fk_internal_faces = flux.sum(axis = 1)
        Fk_vols_total = self.component_flux_volumes(Fk_internal_faces)
        '''Eduarda considera fluxo saindo negativo e entrando positivo?'''

        self.well_term(wells)
        residuo_massa = (fprop.Nk - Nk_old) - delta_t * Fk_vols_total - delta_t * self.q

        # Restricao de volume poroso
        aux = fprop.Nj / fprop.Csi_j
        residuo_poroso_varavei = aux.sum(axis = 1) / fprop.Vp - 1.0
        residuo_poroso_bruno = aux.sum(axis = 1) - fprop.Vp
        # Preferencia: Varavei

        # Residuo total
        residuo = np.append(residuo_poroso_varavei, residuo_massa, axis = 0)
        return residuo

    def transmissibility_FI(self, fprop):
        self.t0_internal_faces_prod = fprop.xkj_internal_faces * \
                                      fprop.Csi_j_internal_faces * \
                                      fprop.mobilities_internal_faces

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

    def jacobian_calculation(self, fprop):
        # Calculo da matriz jacobiana
        d_residuoporoso_dP = self.d_residuoporoso_dP(fprop)
        d_residuoporoso_dNk = self.d_residuoporoso_dNk(fprop)

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

    def d_residuoporoso_dP(self, fprop):
        x = fprop.xkj[0:ctes.Nc,0,:]
        y = fprop.xkj[0:ctes.Nc,1,:]
        Nl = fprop.Nj[0,0,:]
        Nv = fprop.Nj[0,1,:]
        P = fprop.P
        T = fprop.T

        dlnfildP, dlnfildnij, dZodP_parcial, dZldnij_parcial, Zo = \
                self.EOS.get_phase_derivatives(P, T, x, Nl, np.ones(ctes.n_volumes))
        dlnfivdP, dlnfivdnij, dZvdP_parcial, dZvdnij_parcial, Zv = \
                self.EOS.get_phase_derivatives(P, T, y, Nv, np.zeros(ctes.n_volumes))

        Nj_Vp = fprop.Nj / fprop.Vp
        dVw_dP = -ctes.Cw * fprop.Csi_W0 / (fprop.Csi_W ** 2)
        dVo_dP = ctes.R * fprop.T * (fprop.P * dZodP_parcial - Zo) / (fprop.P * fprop.P)
        dVv_dP = ctes.R * fprop.T * (fprop.P * dZvdP_parcial - Zv) / (fprop.P * fprop.P)

        dRv_dP = Nj_Vp[0,2] * dVw_dP + Nj_Vp[0,0] * dVo_dP + Nj_Vp[0,1] * dVv_dP
        return dRv_dP


    def d_residuoporoso_dNk(self, fprop):
        # Essa funcao passará por mudanças para a inclusão do térmico. Ver item C.1 da Tese do Varavei
        dRv_dNk = np.zeros_like(fprop.Nk)
        dRv_dNk[-1] = (1 / fprop.Vp) * (1 / fprop.Csi_W) # agua (obs: duvida se tem o Vp)
        NkT = np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0)

        dZodNk = 1 # IMPLEMENTAR
        dZvdNk = 1 # IMPLEMENTAR
        teste = self.EOS.dA_dxij(fprop)
        dRR_dNk = self.d_RachfordRice_d_Nk(fprop, NkT)
        dRR_dV = self.d_RachfordRice_d_vapour_fraction(fprop)

        import pdb; pdb.set_trace()
        dRv_dNk[:-1] = (1 / fprop.Vp) * ((ctes.R * fprop.T * NkT / fprop.P)*\
                (fprop.L * dZodNk + fprop.V * dZvdNk) + (1 / fprop.Csi_j[0,0])*fprop.L + \
                (1 / fprop.Csi_j[0,1])*fprop.V + NkT*((1 / fprop.Csi_j[0,1]) - (1 / fprop.Csi_j[0,0])) *\
                (dRR_dNk / dRR_dV))

        import pdb; pdb.set_trace()
        pass

    def d_RachfordRice_d_vapour_fraction(self, fprop):
        K_aux = fprop.xkj[:,1]/fprop.xkj[:,0]
        K = K_aux[:-1]
        dRR_dV = - np.sum((K - 1) ** 2 * fprop.z / (1 + fprop.V * (K - 1)) ** 2, axis = 0)
        return dRR_dV

    def d_RachfordRice_d_Nk(self, fprop, NkT):
        K_aux = fprop.xkj[:,1]/fprop.xkj[:,0]
        K = K_aux[:-1]
        RR = np.sum((K - 1) * fprop.z / (1 + fprop.V * (K - 1)), axis = 0)
        dRR_dNk = (1/NkT) * ((K-1)/(1 + fprop.V*(K-1)) - RR)
        return dRR_dNk
