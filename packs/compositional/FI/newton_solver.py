import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from scipy import linalg
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from packs.utils import relative_permeability2, phase_viscosity, capillary_pressure
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
        self.phase_viscosity_class = getattr(phase_viscosity,
        data_loaded['compositional_data']['phase_viscosity'])

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
            jacobiana = self.jacobian_calculation(fprop, Pot_hid, delta_t)

            import pdb; pdb.set_trace()
            Here = True


            jacobiana_test = np.delete(jacobiana, (ctes.n_components + 1) * ctes.bhp_ind, axis=0)
            jacobiana_test = np.delete(jacobiana_test, (ctes.n_components + 1) * ctes.bhp_ind, axis=1)
            residuo_test = np.delete(residuo, (ctes.n_components + 1) * ctes.bhp_ind, axis=0)
            try:
                Jacobian_inv = np.linalg.inv(jacobiana_test)
            except: import pdb; pdb.set_trace()
            delta_x = -Jacobian_inv.dot(residuo_test)

            jacobiana_test2 = jacobiana.copy()
            jacobiana_test2[(ctes.n_components + 1) * ctes.bhp_ind , (ctes.n_components + 1) * ctes.bhp_ind] += 1e16
            try:
                Jacobian_inv_2 = np.linalg.inv(jacobiana_test2)
            except: import pdb; pdb.set_trace()
            delta_x_2 = -Jacobian_inv_2.dot(residuo)
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
        '''Eduarda considera fluxo saindo negativo e entrando positivo'''
        self.well_term(wells)
        residuo_massa = (fprop.Nk - Nk_old) - delta_t * Fk_vols_total - delta_t * self.q

        # Restricao de volume poroso
        aux = fprop.Nj / fprop.Csi_j
        residuo_poroso_varavei = aux.sum(axis = 1) / fprop.Vp - 1.0
        residuo_poroso_bruno = aux.sum(axis = 1) - fprop.Vp
        # Preferencia: Bruno

        # Residuo total
        '''residuo = np.append(residuo_poroso_bruno, residuo_massa, axis = 0)'''
        residuo = np.zeros([ctes.n_volumes*(ctes.n_components + 1),1])
        for i in range(ctes.n_volumes):
            residuo_local = np.zeros([ctes.n_components + 1, 1])
            residuo_local[0] = residuo_poroso_bruno[:,i]
            residuo_local[1:] = residuo_massa[...,i][:,np.newaxis]
            residuo[i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1)] = residuo_local


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

    def jacobian_calculation(self, fprop, Pot_hid, delta_t):
        # Calculo da matriz jacobiana
        d_residuoporoso_dP = self.d_residuoporoso_dP(fprop)
        d_residuoporoso_dNk = self.d_residuoporoso_dNk(fprop)
        d_consvmassa_dP_P, d_consvmassa_dP_E, d_consvmassa_dP_W, d_consvmassa_dNk_P, \
            d_consvmassa_dNk_E, d_consvmassa_dNk_W = self.d_consvmassa_dP_dNk(fprop, Pot_hid, delta_t)


        ''' Lógica só funciona para 1D - Fazer no Scipy e de maneira geral '''
        jacobian = np.zeros([ctes.n_volumes*(ctes.n_components + 1), ctes.n_volumes*(ctes.n_components + 1)])
        for i in range(ctes.n_volumes):
            J_local = np.zeros([ctes.n_components + 1, ctes.n_components + 1])
            J_local[0,0] = d_residuoporoso_dP[i]
            J_local[0,1:] = d_residuoporoso_dNk[...,i]
            J_local[1:,0] = d_consvmassa_dP_P[...,i]
            J_local[1:,1:] = d_consvmassa_dNk_P[...,i]

            jacobian[i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1), \
                     i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1)] = J_local

            if i != 0:
                J_local_W = np.zeros([ctes.n_components + 1, ctes.n_components + 1])
                J_local_W[1:,0] = d_consvmassa_dP_W[...,i]
                J_local_W[1:,1:] = d_consvmassa_dNk_W[...,i]

                jacobian[i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1), \
                         (i-1)*(ctes.n_components + 1) : (i)*(ctes.n_components + 1)] = J_local_W

            if i != ctes.n_volumes - 1:
                J_local_E = np.zeros([ctes.n_components + 1, ctes.n_components + 1])
                J_local_E[1:,0] = d_consvmassa_dP_E[...,i]
                J_local_E[1:,1:] = d_consvmassa_dNk_E[...,i]

                jacobian[i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1), \
                         (i+1)*(ctes.n_components + 1) : (i+2)*(ctes.n_components + 1)] = J_local_E

        return jacobian

    def d_residuoporoso_dP(self, fprop):
        x = fprop.xkj[0:ctes.Nc,0,:]
        y = fprop.xkj[0:ctes.Nc,1,:]
        Nl = fprop.Nj[0,0,:]
        Nv = fprop.Nj[0,1,:]
        Nw = fprop.Nj[0,2,:]
        P = fprop.P
        T = fprop.T

        dlnfildP, dlnfildnij, dZldP_parcial, dZldnij_parcial, Zl = \
                self.EOS.get_phase_derivatives(P, T, x, Nl, np.ones(ctes.n_volumes))
        dlnfivdP, dlnfivdnij, dZvdP_parcial, dZvdnij_parcial, Zv = \
                self.EOS.get_phase_derivatives(P, T, y, Nv, np.zeros(ctes.n_volumes))
        dnldP, dnvdP, dnldNk, dnvdNk, dnildP, dnivdP, dnildNk, dnivdNk = \
                self.EOS.dnij_dNk_dP(dlnfildP, dlnfivdP, dlnfildnij, dlnfivdnij, Nl, Nv)
        dnwdP = np.zeros_like(dnldP)

        dVw_dP = -ctes.Cw * fprop.Csi_W0 / (fprop.Csi_W ** 2)
        dVl_dP = ctes.R * T * (P * dZldP_parcial - Zl) / (P**2)
        dVv_dP = ctes.R * T * (P * dZvdP_parcial - Zv) / (P**2)
        dporosity_dP = ctes.porosity * ctes.Cf

        dRv_dP = dnldP*(1/fprop.Csi_j[0,0]) + Nl*dVl_dP + dnvdP*(1/fprop.Csi_j[0,1]) +\
            Nv*dVv_dP + dnwdP*(1/fprop.Csi_j[0,2]) + Nw*dVw_dP - ctes.Vbulk*dporosity_dP
        return dRv_dP


    def d_residuoporoso_dNk(self, fprop):
        # Essa funcao passará por mudanças para a inclusão do térmico. Ver item C.1 da Tese do Varavei
        dRv_dNk = np.zeros_like(fprop.Nk)
        dRv_dNk[-1] = (1 / fprop.Csi_W) # agua

        # Formulação do Varavei
        '''NkT = np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0) # quando incluir a agua no termico, vai mudar esse calculo
        dZl_dx, dZv_dy = self.EOS.dZ_dxij(fprop)
        dx_dP, dy_dP, dnildP, dnivdP, dnldP, dnvdP, dnildNk, dnivdNk, dnldNk, dnvdNk = \
                    self.EOS.dxkj_dnij_dP(fprop)
        dx_dNk, dy_dNk = self.dxij_dNk(fprop, NkT, dnildNk, dnivdNk, dnldNk, dnvdNk)
        #dZodNk, dZvdNk = self.dZ_dNk(dZl_dx, dZv_dy, dx_dNk, dy_dNk)
        self.dZ_dNk(dZl_dx, dZv_dy, dx_dNk, dy_dNk)
        dRR_dNk = self.d_RachfordRice_d_Nk(fprop, NkT)
        dRR_dV = self.d_RachfordRice_d_vapour_fraction(fprop)

        dRv_dNk[:-1] = ((ctes.R*fprop.T*NkT / fprop.P)*(fprop.L*self.dZldNk + \
                fprop.V*self.dZvdNk) + (1/fprop.Csi_j[0,0])*fprop.L + \
                (1/fprop.Csi_j[0,1])*fprop.V + NkT*((1/fprop.Csi_j[0,1]) - \
                (1/fprop.Csi_j[0,0])) * (dRR_dNk / dRR_dV))'''

        # Formulação do Bruno
        x = fprop.xkj[0:ctes.Nc,0,:]
        y = fprop.xkj[0:ctes.Nc,1,:]
        Nl = fprop.Nj[0,0,:]
        Nv = fprop.Nj[0,1,:]
        P = fprop.P
        T = fprop.T

        dlnfildP, dlnfildnij, dZldP_parcial, dZldnij_parcial, Zl = \
                self.EOS.get_phase_derivatives(P, T, x, Nl, np.ones(ctes.n_volumes))
        dlnfivdP, dlnfivdnij, dZvdP_parcial, dZvdnij_parcial, Zv = \
                self.EOS.get_phase_derivatives(P, T, y, Nv, np.zeros(ctes.n_volumes))
        dnldP, dnvdP, dnldNk, dnvdNk, dnildP, dnivdP, dnildNk, dnivdNk = \
                self.EOS.dnij_dNk_dP(dlnfildP, dlnfivdP, dlnfildnij, dlnfivdnij, Nl, Nv)
        dZldP, dZvdP, dZldNk, dZvdNk = self.EOS.dZ_dP_dNk(dZldP_parcial,
                dZvdP_parcial, dZldnij_parcial, dZvdnij_parcial, dnildP,
                dnivdP, dnildNk, dnivdNk)
        dVl_dNk = (ctes.R * T / P) * dZldNk
        dVv_dNk = (ctes.R * T / P) * dZvdNk

        dRv_dNk[:-1] = dnldNk*(1/fprop.Csi_j[0,0]) + Nl*dVl_dNk + \
                    dnvdNk*(1/fprop.Csi_j[0,1]) + Nv*dVv_dNk
        return dRv_dNk

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

    def dxij_dNk(self, fprop, NkT, dnildNk, dnivdNk, dnldNk, dnvdNk):
        x = fprop.xkj[0:ctes.Nc,0,:]
        y = fprop.xkj[0:ctes.Nc,1,:]
        K = y/x

        # Formulação do Varavei:
        '''
        dx_dNk = np.zeros([ctes.Nc, ctes.Nc, ctes.n_volumes])
        dy_dNk = np.zeros([ctes.Nc, ctes.Nc, ctes.n_volumes])
        for i in range(ctes.Nc):
            dx_dNk[i,] = (1 / (1 + fprop.V*(K[i] - 1))) * (( - fprop.z[i]) / NkT)
            dx_dNk[i,i] = (1 / (1 + fprop.V*(K[i] - 1))) * ((1 - fprop.z[i]) / NkT)
            dy_dNk[i,] = dx_dNk[i,] * K[i]

        # Ajustar isso - fazer pelo menos mais 2 testes
        T = (1 / (1 + fprop.V*(K - 1))) / NkT
        TT = T[:,np.newaxis,:] * np.ones([ctes.Nc, ctes.Nc, ctes.n_volumes])
        I = np.identity(ctes.Nc)[...,np.newaxis] * np.ones([ctes.Nc, ctes.Nc, ctes.n_volumes])
        z_aux = fprop.z[:,np.newaxis] * np.ones([ctes.Nc, ctes.Nc, ctes.n_volumes])
        dx_dNk_teste = (I - z_aux) * TT
        '''
        # Formulação do Bruno:
        dx_dNk = (1 / fprop.Nj[0,0]) * (dnildNk - x[:,np.newaxis]*dnldNk)
        dy_dNk = (1 / fprop.Nj[0,1]) * (dnivdNk - y[:,np.newaxis]*dnvdNk)
        return dx_dNk, dy_dNk

    def dZ_dNk(self, dZl_dx, dZv_dy, dx_dNk, dy_dNk):

        self.dZldNk = (dZl_dx[:,np.newaxis] * dx_dNk[np.newaxis,:,:]).sum(axis=1)
        self.dZvdNk = (dZv_dy[:,np.newaxis] * dy_dNk[np.newaxis,:,:]).sum(axis=1)
        #return dZodNk, dZvdNk

    def d_consvmassa_dP_dNk(self, fprop, Pot_hid, delta_t):
        dCsi_j_dP, dCsi_j_dnij = self.dCsi_j_dP_dnij(fprop)
        dx_dP, dy_dP, dnildP, dnivdP, dnldP, dnvdP, dnildNk, dnivdNk, dnldNk, dnvdNk = \
                    self.EOS.dxkj_dnij_dP(fprop)
        NkT = np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0) # quando incluir a agua no termico, vai mudar esse calculo
        dx_dNk, dy_dNk = self.dxij_dNk(fprop, NkT, dnildNk, dnivdNk, dnldNk, dnvdNk) #[[dx1/dN1, dx1/dN2, ..., dx1/dNnc], [dx2/dN1, dx2/dN2, ..., dx2/dNnc], ..., [dxnc/dN1, dxnc/dN2, ..., dxnc/dNnc]] para cada bloco da malha
        drho_dP, drho_dNk = self.drho_dP_dNk(fprop, dx_dP, dy_dP, dCsi_j_dP, dx_dNk, dy_dNk)
        dx_dnij, dy_dnij = self.dx_dy_dnij(fprop)

        dxij_dP = np.zeros([ctes.n_phases, ctes.Nc + ctes.load_w, ctes.n_volumes])
        dxij_dP[0] = dx_dP
        dxij_dP[1] = dy_dP

        dkrs_dSj, kro, krg, krw = self.dkrs_dSj(fprop)
        relative_permeability = np.array([kro, krg, krw])
        dSj_dP, dSj_dNk = self.dSj_dP_dNk(fprop, dCsi_j_dP, dnldP, dnvdP, dnldNk, dnvdNk)
        dkrj_dP, dkrj_dNk = self.dkrj_dP_dNk(dkrs_dSj, dSj_dP, dSj_dNk)

        Csi_j = fprop.Csi_j.copy()
        phase_viscosity = self.phase_viscosity_class(fprop, Csi_j)
        dmi_dP, dmi_dNk = phase_viscosity.derivative_phase_viscosity_dP_dNk(fprop, dx_dP, dy_dP, dCsi_j_dP, dx_dnij, dy_dnij, dCsi_j_dnij, dnildP, dnivdP, dnildNk, dnivdNk)

        dmobilities_dP = (fprop.mis*dkrj_dP - relative_permeability*dmi_dP)/(fprop.mis**2)
        aux = (fprop.mobilities * fprop.Csi_j)[:,:,np.newaxis] * dxij_dP
        aux = np.transpose(aux, (0,2,1,3))
        dtransmissibility_dP = aux + fprop.xkj * fprop.Csi_j * dmobilities_dP \
                + fprop.xkj * fprop.mobilities * dCsi_j_dP

        # Derivada em funcao de Nk
        dmobilities_dNk = (fprop.mis[:,:,np.newaxis] * dkrj_dNk - \
            relative_permeability[:,np.newaxis] * dmi_dNk)/(fprop.mis[:,:,np.newaxis]**2)

        dtransmissibility_dNk = np.zeros([ctes.Nc + ctes.load_w, ctes.Nc + ctes.load_w, ctes.n_phases, ctes.n_volumes])
        # Só referente a fase óleo:
        aux = np.zeros([ctes.Nc + ctes.load_w, ctes.Nc + ctes.load_w, ctes.n_volumes])
        aux[0:ctes.Nc, 0:ctes.Nc] = dx_dNk * fprop.Csi_j[0,0] * fprop.mobilities[0,0]
        dtransmissibility_dNk[:,:,0,:] = aux + \
            fprop.xkj[:,0][:,np.newaxis] * self.dCsi_j_dNk[0] * fprop.mobilities[0,0] +\
            fprop.xkj[:,0][:,np.newaxis] * fprop.Csi_j[0,0] * dmobilities_dNk[0,0]
        # Só referente a fase gas:
        aux[0:ctes.Nc, 0:ctes.Nc] = dy_dNk * fprop.Csi_j[0,1] * fprop.mobilities[0,1]
        dtransmissibility_dNk[:,:,1,:] = aux + \
            fprop.xkj[:,1][:,np.newaxis] * self.dCsi_j_dNk[1] * fprop.mobilities[0,1] +\
            fprop.xkj[:,1][:,np.newaxis] * fprop.Csi_j[0,1] * dmobilities_dNk[0,1]
        # Só referente a fase agua:
        aux[0:ctes.Nc, 0:ctes.Nc] = 0 * fprop.Csi_j[0,2] * fprop.mobilities[0,2]
        dtransmissibility_dNk[:,:,2,:] = aux + \
            fprop.xkj[:,2][:,np.newaxis] * self.dCsi_j_dNk[2] * fprop.mobilities[0,2] +\
            fprop.xkj[:,2][:,np.newaxis] * fprop.Csi_j[0,2] * dmobilities_dNk[0,2]
        # dtransmissibility_dNk.shape = (Nc1, Nc2, ..., Nci); (dNc1, dNc2, ..., dNci); fases; blocos

        self.dPot_hid_dP(drho_dP, drho_dNk, Pot_hid, fprop)
        transmissibility = self.transmissibility_FI(fprop)

        '''Logica testada para 1D'''
        dtransmissibility_faceup_dP_P = dtransmissibility_dP[:,:,:,ctes.v0[:,0]]
        dtransmissibility_facedown_dP_P = dtransmissibility_dP[:,:,:,ctes.v0[:,1]]
        bool = Pot_hid[:,ctes.v0[:,0]] < Pot_hid[:,ctes.v0[:,1]]
        dtransmissibility_faceup_dP_P[:,:,bool] = 0.0
        dtransmissibility_facedown_dP_P[:,:,~bool] = 0.0

        Pot_hidj_test = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up_test = Pot_hid[:,ctes.v0[:,1]]
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        d_Pot_hid = Pot_hidj_up_test - Pot_hidj_test - ctes.g * fprop.rho_j_internal_faces * (z_up - z)

        dfluxo_faceup_dP_P = dtransmissibility_faceup_dP_P * \
                ctes.pretransmissibility_internal_faces * d_Pot_hid + \
                transmissibility * self.dPot_hid_faceup_dP_P
        dfluxo_facedown_dP_P = dtransmissibility_facedown_dP_P * \
                ctes.pretransmissibility_internal_faces * d_Pot_hid + \
                transmissibility * self.dPot_hid_facedown_dP_P

        # Derivada do termo de fluxo em relação a pressão - diagonal principal
        cx = np.arange(ctes.n_phases)
        lines = np.array([np.repeat(cx,len(ctes.v0[:,0])), np.repeat(cx,len(ctes.v0[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(ctes.v0[:,0],ctes.n_phases), np.tile(ctes.v0[:,1], ctes.n_phases)]).flatten()
        dflux_dP_P = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        for i in range(ctes.n_components):
            data = np.array([dfluxo_faceup_dP_P[:,i], - dfluxo_facedown_dP_P[:,i]]).flatten()
            dflux_dP_P[i] = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()

        # Derivada da eq. de cons. da massa em relação a pressão do bloco (diagonal principal)
        dwell_dP_P = 0.0
        ''' FALTA O TERMO DE POÇO!! ATENÇÃO '''
        d_consvmassa_dP_P = - delta_t * (dflux_dP_P.sum(axis=1) + dwell_dP_P)

        # Blocos vizinhos
        dtransmissibility_faceup_dP_E = dtransmissibility_dP[:,:,:,ctes.v0[:,1]]
        dtransmissibility_facedown_dP_W = dtransmissibility_dP[:,:,:,ctes.v0[:,0]]
        bool = Pot_hid[:,ctes.v0[:,0]] < Pot_hid[:,ctes.v0[:,1]]
        dtransmissibility_faceup_dP_E[:,:,~bool] = 0.0
        dtransmissibility_facedown_dP_W[:,:,bool] = 0.0

        dfluxo_faceup_dP_E = dtransmissibility_faceup_dP_E * \
                ctes.pretransmissibility_internal_faces * d_Pot_hid + \
                transmissibility * self.dPot_hid_faceup_dP_E
        dfluxo_facedown_dP_W = dtransmissibility_facedown_dP_W * \
                ctes.pretransmissibility_internal_faces * d_Pot_hid + \
                transmissibility * self.dPot_hid_facedown_dP_W
        dfluxo_facedown_dP_E = np.zeros_like(dfluxo_faceup_dP_E)
        dfluxo_faceup_dP_W = np.zeros_like(dfluxo_facedown_dP_W)

        # Derivada do termo de fluxo em relação a pressão - blocos vizinhos (E e W)
        dflux_dP_E = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        dflux_dP_W = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        for i in range(ctes.n_components):
            data_E = np.array([dfluxo_faceup_dP_E[:,i], - dfluxo_facedown_dP_E[:,i]]).flatten()
            data_W = np.array([dfluxo_faceup_dP_W[:,i], - dfluxo_facedown_dP_W[:,i]]).flatten()
            dflux_dP_E[i] = sp.csc_matrix((data_E, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()
            dflux_dP_W[i] = sp.csc_matrix((data_W, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()

        # Derivada da eq. de cons. da massa em relação a pressão do bloco vizinho (fora da diagonal principal)
        dwell_dP_E = 0.0
        dwell_dP_W = 0.0
        d_consvmassa_dP_E = - delta_t * (dflux_dP_E.sum(axis=1) + dwell_dP_E)
        d_consvmassa_dP_W = - delta_t * (dflux_dP_W.sum(axis=1) + dwell_dP_P)

        '''Logica testada para 1D'''
        # Derivada em função da composição - diagonal principal
        dtransmissibility_faceup_dNk_P = dtransmissibility_dNk[:,:,:,ctes.v0[:,0]]
        dtransmissibility_facedown_dNk_P = dtransmissibility_dNk[:,:,:,ctes.v0[:,1]]
        bool = Pot_hid[:,ctes.v0[:,0]] < Pot_hid[:,ctes.v0[:,1]]
        dtransmissibility_faceup_dNk_P[:,:,bool] = 0.0
        dtransmissibility_facedown_dNk_P[:,:,~bool] = 0.0

        dfluxo_faceup_dNk_P = dtransmissibility_faceup_dNk_P * \
                ctes.pretransmissibility_internal_faces * d_Pot_hid + \
                transmissibility[:,np.newaxis] * self.dPot_hid_faceup_dNk_P
        dfluxo_facedown_dNk_P = dtransmissibility_facedown_dNk_P * \
                ctes.pretransmissibility_internal_faces * d_Pot_hid + \
                transmissibility[:,np.newaxis] * self.dPot_hid_facedown_dNk_P
        # dfluxo_faceup_dNk_P.shape = (Nc1, Nc2, ..., Nci); (dNc1, dNc2, ..., dNci); fases; faces

        # Derivada do termo de fluxo em relação a composição - diagonal principal
        cx = np.arange(ctes.n_phases)
        lines = np.array([np.repeat(cx,len(ctes.v0[:,0])), np.repeat(cx,len(ctes.v0[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(ctes.v0[:,0],ctes.n_phases), np.tile(ctes.v0[:,1], ctes.n_phases)]).flatten()
        dflux_dNk_P = np.zeros([ctes.n_components, ctes.n_components, ctes.n_phases, ctes.n_volumes])
        for i in range(ctes.n_components):
            for k in range(ctes.n_components):
                data = np.array([dfluxo_faceup_dNk_P[i,k], - dfluxo_facedown_dNk_P[i,k]]).flatten()
                dflux_dNk_P[i,k] = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()

        # Derivada da eq. de cons. da massa em relação a pressão do bloco (diagonal principal)
        dwell_dNk_P = 0.0
        dacumulation_dNk = np.identity(ctes.n_components)[...,np.newaxis] * np.ones([ctes.n_volumes])
        ''' FALTA O TERMO DE POÇO!! ATENÇÃO '''
        #d_consvmassa_dNk_P = - delta_t * (dflux_dNk_P.sum(axis=2) + dwell_dNk_P)
        d_consvmassa_dNk_P = dacumulation_dNk - delta_t * (dflux_dNk_P.sum(axis=2) + dwell_dNk_P)
        #test = d_consvmassa_dNk_P + dacumulation_dNk
        #[[dRc1/dN1, dRc1/dN2, ..., dRc1/dNnc], [dRc2/dN1, dRc2/dN2, ..., dRc2/dNnc], ..., [dRcnc/dN1, dRcnc/dN2, ..., dRcnc/dNnc]] para cada bloco da malha

        # Blocos vizinhos
        dtransmissibility_faceup_dNk_E = dtransmissibility_dNk[:,:,:,ctes.v0[:,1]]
        dtransmissibility_facedown_dNk_W = dtransmissibility_dNk[:,:,:,ctes.v0[:,0]]
        dtransmissibility_faceup_dNk_E[:,:,~bool] = 0.0
        dtransmissibility_facedown_dNk_W[:,:,bool] = 0.0

        dfluxo_faceup_dNk_E = dtransmissibility_faceup_dNk_E * \
                ctes.pretransmissibility_internal_faces * d_Pot_hid + \
                transmissibility[:,np.newaxis] * self.dPot_hid_faceup_dNk_E
        dfluxo_facedown_dNk_W = dtransmissibility_facedown_dNk_W * \
                ctes.pretransmissibility_internal_faces * d_Pot_hid + \
                transmissibility[:,np.newaxis] * self.dPot_hid_facedown_dNk_W
        dfluxo_facedown_dNk_E = np.zeros_like(dfluxo_faceup_dNk_E)
        dfluxo_faceup_dNk_W = np.zeros_like(dfluxo_facedown_dNk_W)

        # Derivada do termo de fluxo em relação a pressão - blocos vizinhos (E e W)
        dflux_dNk_E = np.zeros([ctes.n_components, ctes.n_components, ctes.n_phases, ctes.n_volumes])
        dflux_dNk_W = np.zeros([ctes.n_components, ctes.n_components, ctes.n_phases, ctes.n_volumes])

        for i in range(ctes.n_components):
            for k in range(ctes.n_components):
                data_E = np.array([dfluxo_faceup_dNk_E[i,k], - dfluxo_facedown_dNk_E[i,k]]).flatten()
                data_W = np.array([dfluxo_faceup_dNk_W[i,k], - dfluxo_facedown_dNk_W[i,k]]).flatten()
                dflux_dNk_E[i,k] = sp.csc_matrix((data_E, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()
                dflux_dNk_W[i,k] = sp.csc_matrix((data_W, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()

        # Derivada da eq. de cons. da massa em relação a pressão do bloco vizinho (fora da diagonal principal)
        dwell_dNk_E = 0.0
        dwell_dNk_W = 0.0
        d_consvmassa_dNk_E = - delta_t * (dflux_dNk_E.sum(axis=2) + dwell_dNk_E)
        d_consvmassa_dNk_W = - delta_t * (dflux_dNk_W.sum(axis=2) + dwell_dNk_W)


        return d_consvmassa_dP_P, d_consvmassa_dP_E, d_consvmassa_dP_W, d_consvmassa_dNk_P, \
                d_consvmassa_dNk_E, d_consvmassa_dNk_W

    def dCsi_j_dP_dnij(self, fprop):
        x = fprop.xkj[0:ctes.Nc,0,:]
        y = fprop.xkj[0:ctes.Nc,1,:]
        Nl = fprop.Nj[0,0,:]
        Nv = fprop.Nj[0,1,:]
        P = fprop.P
        T = fprop.T

        dlnfildP, dlnfildnij, dZldP_parcial, dZldnij_parcial, Zl = \
                self.EOS.get_phase_derivatives(P, T, x, Nl, np.ones(ctes.n_volumes))
        dlnfivdP, dlnfivdnij, dZvdP_parcial, dZvdnij_parcial, Zv = \
                self.EOS.get_phase_derivatives(P, T, y, Nv, np.zeros(ctes.n_volumes))
        dnldP, dnvdP, dnldNk, dnvdNk, dnildP, dnivdP, dnildNk, dnivdNk = \
                self.EOS.dnij_dNk_dP(dlnfildP, dlnfivdP, dlnfildnij, dlnfivdnij, Nl, Nv)
        dZldP, dZvdP, dZldNk, dZvdNk = self.EOS.dZ_dP_dNk(dZldP_parcial,
                dZvdP_parcial, dZldnij_parcial, dZvdnij_parcial, dnildP,
                dnivdP, dnildNk, dnivdNk)

        dCsi_j_dP = np.zeros_like(fprop.Csi_j)
        dCsi_j_dP[0,0] = (Zl - P*dZldP_parcial) / (ctes.R*T*(Zl**2))
        dCsi_j_dP[0,1] = (Zv - P*dZvdP_parcial) / (ctes.R*T*(Zv**2))
        dCsi_j_dP[0,2] = fprop.Csi_W * ctes.Cw

        dCsi_j_dnij = np.zeros([ctes.n_phases, ctes.Nc, ctes.n_volumes])
        dCsi_j_dnij[0] = P * dZldnij_parcial / (ctes.R*T*(Zl**2))
        dCsi_j_dnij[1] = P * dZvdnij_parcial / (ctes.R*T*(Zv**2))
        '''
        self.dCsi_j_dNk = np.zeros([ctes.n_phases, ctes.Nc + ctes.load_w, ctes.n_volumes])
        self.dCsi_j_dNk[0,0:ctes.Nc] = -P/(ctes.R*T*(Zl**2)) * self.dZldNk
        self.dCsi_j_dNk[1,0:ctes.Nc] = -P/(ctes.R*T*(Zv**2)) * self.dZvdNk
        '''
        self.dCsi_j_dNk = np.zeros([ctes.n_phases, ctes.Nc + ctes.load_w, ctes.n_volumes])
        self.dCsi_j_dNk[0,0:ctes.Nc] = -P/(ctes.R*T*(Zl**2)) * dZldNk
        self.dCsi_j_dNk[1,0:ctes.Nc] = -P/(ctes.R*T*(Zv**2)) * dZvdP
        # Posso testar a formulação do Bruno
        return dCsi_j_dP, dCsi_j_dnij

    def drho_dP_dNk(self, fprop, dx_dP, dy_dP, dCsi_j_dP, dx_dNk, dy_dNk):
        x = fprop.xkj[0:ctes.Nc,0,:]
        y = fprop.xkj[0:ctes.Nc,1,:]

        drho_dP = np.zeros_like(fprop.rho_j)
        drho_dP[0,0] = dCsi_j_dP[0,0] * (x * ctes.Mw[:, np.newaxis]).sum(axis=0) + \
                fprop.Csi_j[0,0] * (ctes.Mw[:, np.newaxis] * dx_dP[0:ctes.Nc]).sum(axis=0)
        drho_dP[0,1] = dCsi_j_dP[0,1] * (y * ctes.Mw[:, np.newaxis]).sum(axis=0) + \
                fprop.Csi_j[0,1] * (ctes.Mw[:, np.newaxis] * dy_dP[0:ctes.Nc]).sum(axis=0)
        drho_dP[0,2] = ctes.Mw_w * dCsi_j_dP[0,2]

        drho_dNk = np.zeros([ctes.n_phases, ctes.Nc + ctes.load_w, ctes.n_volumes])
        #dx_dNk_new = np.zeros([ctes.Nc + ctes.load_w, ctes.Nc + ctes.load_w, ctes.n_volumes])
        dxkj_dNk = np.zeros([ctes.n_phases, ctes.Nc + ctes.load_w, \
                             ctes.Nc + ctes.load_w, ctes.n_volumes])
        #dx_dNk_new[0:ctes.Nc, 0:ctes.Nc, :] = dx_dNk
        dxkj_dNk[0, 0:ctes.Nc, 0:ctes.Nc, :] = dx_dNk
        dxkj_dNk[1, 0:ctes.Nc, 0:ctes.Nc, :] = dy_dNk
        Mw = np.append(ctes.Mw, ctes.Mw_w)
        aux_oil = (Mw * dxkj_dNk[0].T).sum(axis = 2)
        aux_gas = (Mw * dxkj_dNk[1].T).sum(axis = 2)

        drho_dNk[0] = self.dCsi_j_dNk[0] * (x * ctes.Mw[:,np.newaxis]).sum(axis=0) + \
                        fprop.Csi_j[0,0] * aux_oil.T
        drho_dNk[1] = self.dCsi_j_dNk[1] * (y * ctes.Mw[:,np.newaxis]).sum(axis=0) + \
                        fprop.Csi_j[0,1] * aux_gas.T

        return drho_dP, drho_dNk

    def dkrs_dSj(self, fprop):
        self.relative_permeability_class = getattr(relative_permeability2,
        data_loaded['compositional_data']['relative_permeability'])
        self.relative_permeability = self.relative_permeability_class()

        So = fprop.So.copy()
        Sg = fprop.Sg.copy()
        Sw = fprop.Sw.copy()
        saturations = np.array([So, Sg, Sw])
        kro,krg,krw, Sor = self.relative_permeability(fprop, saturations)
        krj = np.zeros([1, ctes.n_phases, len(So)])
        if ctes.load_k:
            krj[0,0,:] = kro
            krj[0,1,:] = krg
        if ctes.load_w:
            krj[0, ctes.n_phases-1,:] = krw
        dkrsdSj = self.relative_permeability.dkrs_dSj(krj, saturations)
        return dkrsdSj, kro, krg, krw

    def dSj_dP_dNk(self, fprop, dCsi_j_dP, dnldP, dnvdP, dnldNk, dnvdNk):
        dSw_dP = (- 1 / ((fprop.Vp * fprop.Csi_W)**2)) * (fprop.Nj[0,-1]*\
            fprop.Vp * dCsi_j_dP[0,-1] + fprop.Nj[0,-1]*fprop.Csi_W*ctes.Vbulk*ctes.porosity*ctes.Cf)

        x = fprop.xkj[0:ctes.Nc,0,:]
        y = fprop.xkj[0:ctes.Nc,1,:]
        Nl = fprop.Nj[0,0,:]
        Nv = fprop.Nj[0,1,:]
        P = fprop.P
        T = fprop.T
        NkT = np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0) # quando incluir a agua no termico, vai mudar esse calculo

        dlnfildP, dlnfildnij, dZodP_parcial, dZldnij_parcial, Zo = \
                self.EOS.get_phase_derivatives(P, T, x, Nl, np.ones(ctes.n_volumes))
        dlnfivdP, dlnfivdnij, dZvdP_parcial, dZvdnij_parcial, Zv = \
                self.EOS.get_phase_derivatives(P, T, y, Nv, np.zeros(ctes.n_volumes))
        '''
        # Formulação do Varavei:
        dSo_dP = (ctes.R * fprop.T * NkT / fprop.P)*fprop.L * (dZodP_parcial - Zo/fprop.P)
        dSv_dP = (ctes.R * fprop.T * NkT / fprop.P)*fprop.V * (dZvdP_parcial - Zv/fprop.P)

        dNoCsio_dP = (dnldP * fprop.Csi_j[0,0] - fprop.Nj[0,0] * dCsi_j_dP[0,0]) / \
                        (fprop.Csi_j[0,0]**2)
        dNvCsio_dP = (dnvdP * fprop.Csi_j[0,1] - fprop.Nj[0,1] * dCsi_j_dP[0,1]) / \
                        (fprop.Csi_j[0,1]**2)
        beta = np.sum(fprop.Nj / fprop.Csi_j, axis = 1)

        dSo_dP = (dNoCsio_dP*beta - (fprop.Nj[0,0]/fprop.Csi_j[0,0]) * \
                  (dNoCsio_dP + dNvCsio_dP)) / (beta**2)
        dSv_dP = (dNvCsio_dP*beta - (fprop.Nj[0,1]/fprop.Csi_j[0,1]) * \
                  (dNoCsio_dP + dNvCsio_dP)) / (beta**2)
        '''
        # Formulação do Bruno:
        NjT = np.sum(fprop.Nj[:,0:ctes.n_phases-1,:], axis = 1)
        dL_dP = (1 / NjT) * (dnldP - fprop.L*(dnldP + dnvdP))
        dV_dP = (1 / NjT) * (dnvdP - fprop.V*(dnldP + dnvdP))
        dL_dNk = (1 / NjT) * (dnldNk - fprop.L*(dnldNk + dnvdNk))
        dV_dNk = (1 / NjT) * (dnvdNk - fprop.V*(dnldNk + dnvdNk))

        dLCsi_dP = (1/fprop.Csi_j[0,0])*(dL_dP - (fprop.L/fprop.Csi_j[0,0])*dCsi_j_dP[0,0])
        dVCsi_dP = (1/fprop.Csi_j[0,1])*(dV_dP - (fprop.V/fprop.Csi_j[0,1])*dCsi_j_dP[0,1])
        dLCsi_dNk = (1/fprop.Csi_j[0,0])*(dL_dNk - (fprop.L/fprop.Csi_j[0,0])*self.dCsi_j_dNk[0,0:ctes.Nc])
        dVCsi_dNk = (1/fprop.Csi_j[0,1])*(dV_dNk - (fprop.V/fprop.Csi_j[0,1])*self.dCsi_j_dNk[1,0:ctes.Nc])

        dSo_dP = (1 - fprop.Sw) * ((dLCsi_dP * (fprop.L / fprop.Csi_j[0,0] + \
             fprop.V / fprop.Csi_j[0,1]) - (fprop.L / fprop.Csi_j[0,0])*(dLCsi_dP + dVCsi_dP)) / \
             ((fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1])**2)) - \
             ((fprop.L / fprop.Csi_j[0,0]) / (fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1]))*dSw_dP
        dSv_dP = (1 - fprop.Sw) * ((dVCsi_dP * (fprop.L / fprop.Csi_j[0,0] + \
             fprop.V / fprop.Csi_j[0,1]) - (fprop.V / fprop.Csi_j[0,1])*(dLCsi_dP + dVCsi_dP)) / \
             ((fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1])**2)) - \
             ((fprop.V / fprop.Csi_j[0,1]) / (fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1]))*dSw_dP

        dSj_dP = np.zeros_like(fprop.Csi_j)
        dSj_dP[0,0] = dSo_dP
        dSj_dP[0,1] = dSv_dP
        dSj_dP[0,2] = dSw_dP

        dSj_dNk = np.zeros([ctes.n_phases, ctes.Nc + ctes.load_w, ctes.n_volumes])
        dSj_dNk[2,-1] = 1 / (fprop.Csi_j[:,-1] * fprop.Vp)
        dSj_dNk[0,0:ctes.Nc] = (1 - fprop.Sw) * ((dLCsi_dNk * (fprop.L / fprop.Csi_j[0,0] + \
             fprop.V / fprop.Csi_j[0,1]) - (fprop.L / fprop.Csi_j[0,0])*(dLCsi_dNk + dVCsi_dNk)) / \
             ((fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1])**2))
        dSj_dNk[1,0:ctes.Nc] = (1 - fprop.Sw) * ((dVCsi_dNk * (fprop.L / fprop.Csi_j[0,0] + \
             fprop.V / fprop.Csi_j[0,1]) - (fprop.V / fprop.Csi_j[0,1])*(dLCsi_dNk + dVCsi_dNk)) / \
             ((fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1])**2))

        return dSj_dP, dSj_dNk

    def dkrj_dP_dNk(self, dkrs_dSj, dSj_dP, dSj_dNk):
        dkrj_dP = np.zeros_like(dSj_dP)
        dkrj_dNk = np.zeros_like(dSj_dNk)
        # Formulação do Varavei:
        dkrj_dP[0,0] = dkrs_dSj[0,0,1] * dSj_dP[0,1] + dkrs_dSj[0,0,2] * dSj_dP[0,2]
        dkrj_dP[0,1] = dkrs_dSj[0,1,1] * dSj_dP[0,1]
        dkrj_dP[0,2] = dkrs_dSj[0,2,2] * dSj_dP[0,2]

        dkrj_dNk[0] = dkrs_dSj[:,2,2] * dSj_dNk[2] + dkrs_dSj[:,0,1] * dSj_dNk[1]
        dkrj_dNk[1] = dkrs_dSj[:,1,1] * dSj_dNk[1]
        dkrj_dNk[2] = dkrs_dSj[:,2,2] * dSj_dNk[2]
        '''
        # Formulação do Bruno:
        dkrj_dP[0,0] = dkrs_dSj[0,0,0]*dSj_dP[0,0] + dkrs_dSj[0,0,1]*dSj_dP[0,1] + dkrs_dSj[0,0,2]*dSj_dP[0,2]
        dkrj_dP[0,1] = dkrs_dSj[0,1,0]*dSj_dP[0,0] + dkrs_dSj[0,1,1]*dSj_dP[0,1] + dkrs_dSj[0,1,2]*dSj_dP[0,2]
        dkrj_dP[0,2] = dkrs_dSj[0,2,0]*dSj_dP[0,0] + dkrs_dSj[0,2,1]*dSj_dP[0,1] + dkrs_dSj[0,2,2]*dSj_dP[0,2]
        '''
        return dkrj_dP, dkrj_dNk

    def dx_dy_dnij(self, fprop):
        x = fprop.xkj[0:ctes.Nc,0,:]
        y = fprop.xkj[0:ctes.Nc,1,:]
        I = np.identity(ctes.Nc)[...,np.newaxis] * \
            np.ones([ctes.Nc, ctes.Nc, ctes.n_volumes])
        dx_dnij = (I - x[:,np.newaxis,:]*np.ones([ctes.Nc, ctes.Nc, ctes.n_volumes])) \
                    / fprop.Nj[0,0]
        dy_dnij = (I - y[:,np.newaxis,:]*np.ones([ctes.Nc, ctes.Nc, ctes.n_volumes])) \
                    / fprop.Nj[0,1]
        return dx_dnij, dy_dnij

    def dPot_hid_dP(self, drho_dP, drho_dNk, Pot_hid, fprop):
        '''Testado para 1D'''
        #dPot_hid_dP = np.ones_like(drho_dP) - drho_dP * ctes.g * ctes.z
        #dPot_hid_dNk = ctes.g * ctes.z * drho_dNk
        ones = np.ones_like(drho_dP[:,:,ctes.v0[:,0]])
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]

        #self.dPot_hid_faceup_dP_P = - dPot_hid_dP[:,:,ctes.v0[:,0]]
        self.dPot_hid_faceup_dP_P = - ones - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,0]] \
            / (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) * drho_dP[:,:,ctes.v0[:,0]]
        #self.dPot_hid_facedown_dP_P = dPot_hid_dP[:,:,ctes.v0[:,1]]
        self.dPot_hid_facedown_dP_P = ones - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,1]] \
            / (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) * drho_dP[:,:,ctes.v0[:,1]]
        #self.dPot_hid_faceup_dP_E = dPot_hid_dP[:,:,ctes.v0[:,1]]
        self.dPot_hid_faceup_dP_E = ones - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,1]] \
            / (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) * drho_dP[:,:,ctes.v0[:,1]]
        #self.dPot_hid_facedown_dP_W = - dPot_hid_dP[:,:,ctes.v0[:,0]]
        self.dPot_hid_facedown_dP_W = - ones - ctes.g *(z_up - z) *(fprop.Vp[ctes.v0[:,0]] \
            / (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) * drho_dP[:,:,ctes.v0[:,0]]

        #self.dPot_hid_faceup_dNk_P = dPot_hid_dNk[:,:,ctes.v0[:,0]]
        self.dPot_hid_faceup_dNk_P = - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,0]] \
            /(fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) *drho_dNk[:,:,ctes.v0[:,0]]
        self.dPot_hid_faceup_dNk_P = np.transpose(self.dPot_hid_faceup_dNk_P, (1,0,2))
        #self.dPot_hid_facedown_dNk_P = - dPot_hid_dNk[:,:,ctes.v0[:,1]]
        self.dPot_hid_facedown_dNk_P = - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,1]] \
            /(fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) *drho_dNk[:,:,ctes.v0[:,1]]
        self.dPot_hid_facedown_dNk_P = np.transpose(self.dPot_hid_facedown_dNk_P, (1,0,2))
        #self.dPot_hid_faceup_dNk_E = - dPot_hid_dNk[:,:,ctes.v0[:,1]]
        self.dPot_hid_faceup_dNk_E = - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,1]] \
            /(fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) *drho_dNk[:,:,ctes.v0[:,1]]
        self.dPot_hid_faceup_dNk_E = np.transpose(self.dPot_hid_faceup_dNk_E, (1,0,2))
        #self.dPot_hid_facedown_dNk_W = dPot_hid_dNk[:,:,ctes.v0[:,0]]
        self.dPot_hid_facedown_dNk_W = - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,0]] \
            /(fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) *drho_dNk[:,:,ctes.v0[:,0]]
        self.dPot_hid_facedown_dNk_W = np.transpose(self.dPot_hid_facedown_dNk_W, (1,0,2))
