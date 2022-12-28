import numpy as np
from packs.directories import data_loaded
from packs.utils import constants as ctes
from scipy import linalg
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from packs.utils import relative_permeability2, phase_viscosity, capillary_pressure
import time
from scipy import sparse
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import cg
from packs.directories import data_loaded


class NewtonSolver:
    def __init__(self, fprop):

        if ctes.load_k: self.EOS = ctes.EOS_class(fprop.T)
        self.phase_viscosity_class = getattr(phase_viscosity,
        data_loaded['compositional_data']['phase_viscosity'])
        if data_loaded['compositional_data']['relative_permeability'] == 'StoneII':
            self.krw0, self.krow0, self.krog0, self.krg0, self.n_og, self.n_ow, self.n_w, self.n_g = \
                            data_loaded['compositional_data']['relative_permeability_data'].values()
            self.Sorw, self.Sorg, self.Swr, self.Sgr = \
                            data_loaded['compositional_data']['residual_saturations'].values()

    def solver(self, wells, fprop, delta_t, Nk_old, G, StabilityCheck, p1, M, face_properties, phase_densities, dVjdNk, dVjdP, P_old):

        self.well_term(wells)
        Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]
        stop_criteria = 1e-6 # TESTAR
        contador = 0
        contador_global = 0
        wp = wells['ws_p']

        residuo = self.residual_calculation(fprop, Pot_hid, wells, Nk_old, delta_t)

        # Adimensionalização dos residuos
        residuo_adm = residuo.copy()
        """bool_inds = np.ones_like(residuo,dtype=bool)
        bool_inds[0:len(residuo):ctes.n_components+1] = False
        residuo_adm[bool_inds] = residuo[bool_inds] / Nk_old.T.flatten()
        aux = residuo_adm[bool_inds]
        aux[Nk_old.T.flatten() == 0] = 0.0
        residuo_adm[bool_inds] = aux
        residuo_adm[~bool_inds] = residuo[~bool_inds] / fprop.Vp"""

        while np.max(abs(residuo_adm)) > stop_criteria:

            ''' Rever para usar matriz esparsa '''
            jacobiana = self.jacobian_calculation(fprop, Pot_hid, delta_t, wells)

            '---------------------- Solução do sistema --------------------------'
            jacobiana_boundary = jacobiana.copy()
            jacobiana_boundary[(ctes.n_components + 1) * wp , :] = 0.0
            #jacobiana_boundary[(ctes.n_components + 1) * ctes.bhp_ind , :] = 0.0
            jacobiana_boundary[(ctes.n_components + 1) * wp , (ctes.n_components + 1) * wp] = 1.0
            #jacobiana_boundary[(ctes.n_components + 1) * ctes.bhp_ind , (ctes.n_components + 1) * ctes.bhp_ind] = 1.0

            residuo_boundary = residuo.copy()
            residuo_boundary[(ctes.n_components + 1) * wp] = 0.0
            #residuo_boundary[(ctes.n_components + 1) * ctes.bhp_ind] = 0.0

            Jb = sparse.csr_matrix(jacobiana_boundary)
            delta_x = spsolve(Jb, -residuo_boundary)


            '--------------------- Atualização das variáveis ---------------------'
            detla_x_new = np.reshape(delta_x, (ctes.n_volumes, ctes.n_components + 1)).T
            # Gambiarra
            #detla_x_new[0] = detla_x_new[0]/2

            """
            '----------------- Pre condicionador de Jacobi ----------------------'
            Jacobi_PC = np.zeros_like(jacobiana_boundary)
            for i in range(ctes.n_volumes):
                Jacobi_PC[(ctes.n_components+1)*i:(ctes.n_components+1)*(i+1), (ctes.n_components+1)*i:(ctes.n_components+1)*(i+1)] = \
                    jacobiana_boundary[(ctes.n_components+1)*i:(ctes.n_components+1)*(i+1), (ctes.n_components+1)*i:(ctes.n_components+1)*(i+1)]
            jacobiana_PC = np.linalg.inv(Jacobi_PC).dot(jacobiana_boundary)
            #Jb_PC = sparse.csr_matrix(jacobiana_PC)
            residuo_PC = np.linalg.inv(Jacobi_PC).dot(residuo)

            #delta_x_2 = spsolve(Jb_test, -residuo_test)
            delta_x_Jacobi = np.linalg.solve(jacobiana_PC, -residuo_PC) # O MAIS RÁPIDO
            detla_x_new = np.reshape(delta_x_Jacobi, (ctes.n_volumes, ctes.n_components + 1)).T
            #delta_x_Jacobi_GMRES, code = gmres(teste2, -residuo_test)
            #delta_x_Jacobi_CG, code = cg(Jb_test, -residuo_test)
            #detla_x_new = np.reshape(delta_x_Jacobi_CG, (ctes.n_volumes, ctes.n_components + 1)).T
            'Gambiarra'
            detla_x_new[0] = detla_x_new[0]/2
            """


            detla_x_new[1:, wp] = 0.0 # Considerar que:
            # Para a hipotese de considerar o bloco de P prescrita calculado separadamente
            # Para o caso 6K, usar essa lógica
            # Para o caso do BL = comentar

            fprop.P += detla_x_new[0]
            fprop.Nk += detla_x_new[1:]
            fprop.Nk[fprop.Nk < -1e-300] = 0.0
            fprop.Nk[fprop.Nk < 1e-20] = 0.0
            if any(fprop.Nk.flatten()<0): import pdb; pdb.set_trace()
            if ctes.load_k: fprop.z = fprop.Nk[0:ctes.Nc,:] / fprop.Nk[0:ctes.Nc,:].sum(axis=0)

            '----------------- Perform Phase stability test and flash -------------'
            #if fprop.Sg[0]<1: import pdb; pdb.set_trace()
            if any(fprop.Sg>1): import pdb; pdb.set_trace()
            if ctes.load_k and ctes.compressible_k:
                #if self.z
                self.p2 = StabilityCheck(fprop.P, fprop.T)
                fprop.L, fprop.V, fprop.xkj[0:ctes.Nc, 0, :], \
                fprop.xkj[0:ctes.Nc, 1, :], fprop.Csi_j[:,0,:], \
                fprop.Csi_j[:,1,:], fprop.rho_j[:,0,:], fprop.rho_j[:,1,:]  =  \
                self.p2.run_init(fprop.P, np.copy(fprop.z))#, wells)

                if len(wells['ws_inj'])>0 and not self.p2.constant_K:

                    if any(wells['inj_cond']=='reservoir'):
                        injP_bool = wells['ws_p'] == wells['ws_inj']
                        injP = wells['ws_p'][injP_bool]

                        z = (wells['z']).T

                        p_well = StabilityCheck(fprop.P[wells['ws_inj']], fprop.T)
                        L, V, x, y, Csi_L, Csi_V, rho_L, rho_V  =  \
                        p_well.run_init(fprop.P[wells['ws_inj']],z[0:ctes.Nc])

                        """
                        # Nao sei se precisa disso!
                        if any(injP_bool):
                            qk_molar = fprop.qk_molar[:,injP]
                            #wells['values_q_vol'] = np.zeros_like((ctes.n_components,len(injP)))
                            wells['values_q_vol'] = qk_molar / \
                                ((Csi_V * V + Csi_L * L))#[wells['inj_cond']=='reservoir']
                        else:
                            wells['values_q'][:,wells['inj_cond']=='reservoir'] = (Csi_V * V + Csi_L * L) * self.q_vol
                            #wells['values_q_vol'][:,wells['inj_cond']=='reservoir'] = self.q_vol
                        """
                    #else:
                        #wells['inj_p_term'] = []
                        #qk_molar = wells['values_q'][:,wells['inj_cond']=='surface']
                        #wells['values_q_vol'][:,wells['inj_cond']=='surface'] = qk_molar / \
                        #    ((Csi_V * V + Csi_L * L))[wells['inj_cond']=='surface']

            #if any(fprop.L!=0): import pdb; pdb.set_trace()


            '----------------------- Update fluid properties ----------------------'
            p1.run_inside_loop(M, fprop)
            G = self.update_gravity_term(fprop)
            Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]
            face_properties(fprop, G)
            phase_densities(fprop)


            '----- Calcular termo de poço com pressão prescrita caso nao adote o bloco separado----'
            wp = wells['ws_p']
            transmissibility = self.transmissibility_FI(fprop)
            Pot_hidj_new = fprop.P[ctes.v0[:,0]]
            Pot_hidj_up_new = fprop.P[ctes.v0[:,1]]
            z = ctes.z[ctes.v0[:,0]]
            z_up = ctes.z[ctes.v0[:,1]]
            d_Pot_hid_new = Pot_hidj_up_new - Pot_hidj_new - ctes.g * fprop.rho_j_internal_faces * (z_up - z)
            flux_new = transmissibility * d_Pot_hid_new
            Fk_internal_faces_new = flux_new.sum(axis = 1)
            Fk_vols_total_new = self.component_flux_volumes(Fk_internal_faces_new)

            mob_ratio = fprop.mobilities[:,:,wp] / np.sum(fprop.mobilities[:,:,wp], axis = 1)
            q_term = fprop.xkj[:,:,wp] * mob_ratio * fprop.Csi_j[:,:,wp]
            well_term = -(fprop.Vt[wp] - fprop.Vp[wp])/delta_t

            #self.q[:,wp] += np.sum(q_term * well_term, axis = 1)
            #self.q[:,wp] = - Fk_vols_total_new[:,wp]
            # Fica comentado caso o bloco de P prescrita calculado separadamente
            # Para o caso 6K, comentar
            # Para o caso do BL, calcular dessa forma


            contador += 1

            residuo = self.residual_calculation(fprop, Pot_hid, wells, Nk_old, delta_t)

            # Adimensionalização dos residuos
            residuo_adm = residuo.copy()
            """bool_inds = np.ones_like(residuo,dtype=bool)
            bool_inds[0:len(residuo):ctes.n_components+1] = False
            residuo_adm[bool_inds] = residuo[bool_inds] / Nk_old.T.flatten()
            aux = residuo_adm[bool_inds]
            aux[Nk_old.T.flatten() == 0] = 0.0
            residuo_adm[bool_inds] = aux
            residuo_adm[~bool_inds] = residuo[~bool_inds] / fprop.Vp"""


            residuo_adm[(ctes.n_components + 1) * wp] = 0.0 # Considerar que:
            #residuo_adm[(ctes.n_components + 1) * ctes.bhp_ind] = 0.0 # Considerar que:
            for j in wp:
                it = int((ctes.n_components + 1) * j) # Considerar que:
                #it = int((ctes.n_components + 1) * ctes.bhp_ind) # Considerar que:
                residuo_adm[it:it + (ctes.n_components + 1)] = 0.0 # Considerar que:
            # Para a hipotese de considerar o bloco de P prescrita calculado separadamente
            # Para o caso 6K, usar essa lógica
            # Para o caso do BL = comentar


            if contador > 100:
                import pdb; pdb.set_trace()

        print(f'{contador} iterações de Newton')


        '----- Calcular termo de poço com pressão prescrita caso o bloco seja separado----'
        dVjdNk, dVjdP = self.dVt_derivatives(fprop)
        self.dVtP = dVjdP.sum(axis=1)[0]
        self.dVtk = dVjdNk.sum(axis=1)

        #well_term = -(fprop.Vt[wp] - fprop.Vp[wp])/delta_t # terceira alternativa ao calculo do poço
        well_term = ((ctes.Vbulk[wp] * ctes.porosity[wp] * ctes.Cf - self.dVtP[wp]) * (fprop.P[wp] - P_old[wp]) + \
                        (fprop.Vp[wp] - fprop.Vt[wp]))/delta_t - (self.dVtk[:,wp] * Fk_vols_total_new[:,wp]).sum(axis=0)

        """if (len(wp)>1) and (well_term[0]<0):
            well_term = -well_term
            import pdb; pdb.set_trace()"""


        mob_ratio = fprop.mobilities[:,:,wp] / np.sum(fprop.mobilities[:,:,wp], axis = 1)
        q_term = fprop.xkj[:,:,wp] * mob_ratio * fprop.Csi_j[:,:,wp]

        self.q[:,wp] = np.sum(q_term * well_term, axis = 1)

        fprop.Nk[:,wp] = Nk_old[:,wp] + delta_t * (Fk_vols_total_new[:,wp] + self.q[:,wp])
        fprop.Nk[fprop.Nk < -1e-300] = 0.0
        fprop.Nk[fprop.Nk < 1e-10] = 0.0
        if any(fprop.Nk.flatten()<0): import pdb; pdb.set_trace()
        fprop.z = fprop.Nk[0:ctes.Nc,:] / fprop.Nk[0:ctes.Nc,:].sum(axis=0)

        '----------------- Perform Phase stability test and flash -------------'
        #if fprop.Sg[0]<1: import pdb; pdb.set_trace()
        if any(fprop.Sg>1): import pdb; pdb.set_trace()
        if ctes.load_k and ctes.compressible_k:

            self.p2 = StabilityCheck(fprop.P, fprop.T)
            fprop.L, fprop.V, fprop.xkj[0:ctes.Nc, 0, :], \
            fprop.xkj[0:ctes.Nc, 1, :], fprop.Csi_j[:,0,:], \
            fprop.Csi_j[:,1,:], fprop.rho_j[:,0,:], fprop.rho_j[:,1,:]  =  \
            self.p2.run_init(fprop.P, np.copy(fprop.z))#, wells)

            if len(wells['ws_inj'])>0 and not self.p2.constant_K:

                if any(wells['inj_cond']=='reservoir'):
                    injP_bool = wells['ws_p'] == wells['ws_inj']
                    injP = wells['ws_p'][injP_bool]

                    z = (wells['z']).T

                    p_well = StabilityCheck(fprop.P[wells['ws_inj']], fprop.T)
                    L, V, x, y, Csi_L, Csi_V, rho_L, rho_V  =  \
                    p_well.run_init(fprop.P[wells['ws_inj']],z[0:ctes.Nc])

                    """
                    if any(injP_bool):
                        qk_molar = self.q[:,injP]
                        #qk_molar = fprop.qk_molar[:,injP]
                        #wells['values_q_vol'] = np.zeros_like((ctes.n_components,len(injP)))
                        wells['values_q_vol'] = qk_molar / \
                            ((Csi_V * V + Csi_L * L))#[wells['inj_cond']=='reservoir']
                    else:
                        wells['values_q'][:,wells['inj_cond']=='reservoir'] = (Csi_V * V + Csi_L * L) * self.q_vol
                        #wells['values_q_vol'][:,wells['inj_cond']=='reservoir'] = self.q_vol
                    """
                #else:
                    #wells['inj_p_term'] = []
                    #qk_molar = wells['values_q'][:,wells['inj_cond']=='surface']
                    #wells['values_q_vol'][:,wells['inj_cond']=='surface'] = qk_molar / \
                    #    ((Csi_V * V + Csi_L * L))[wells['inj_cond']=='surface']

        #if any(fprop.L!=0): import pdb; pdb.set_trace()

        '----------------------- Update fluid properties ----------------------'
        p1.run_inside_loop(M, fprop)
        G = self.update_gravity_term(fprop)
        Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]
        face_properties(fprop, G)
        phase_densities(fprop)
        '--------------------------- FIM do termo de poço --------------------------'



        # Retornar termos que são necessários para o cálculo do CFL
        wave_velocity = Fk_vols_total_new/fprop.Nk
        wave_velocity[fprop.Nk==0] = 0.0
        if any(np.isnan(wave_velocity.flatten())) or any(np.isinf(wave_velocity.flatten())): import pdb; pdb.set_trace()

        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        Ft_internal_faces = - np.sum(fprop.mobilities_internal_faces
            * ctes.pretransmissibility_internal_faces * ((Pot_hidj_up_new - Pot_hidj_new) -
            ctes.g * fprop.rho_j_internal_faces * (z_up - z)), axis = 1)

        fprop.qk_molar = self.q
        #import pdb; pdb.set_trace()
        return Fk_vols_total_new, wave_velocity, Ft_internal_faces


    def residual_calculation(self, fprop, Pot_hid, wells, Nk_old, delta_t):
        # Equacao de conservacao da massa
        transmissibility = self.transmissibility_FI(fprop)
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        d_Pot_hid = Pot_hidj_up - Pot_hidj - ctes.g * fprop.rho_j_internal_faces * (z_up - z)

        flux = transmissibility * d_Pot_hid
        Fk_internal_faces = flux.sum(axis = 1)
        Fk_vols_total = self.component_flux_volumes(Fk_internal_faces)
        '''Eduarda considera fluxo saindo negativo e entrando positivo'''

        residuo_massa = (fprop.Nk - Nk_old) - delta_t * (Fk_vols_total + self.q)

        # Restricao de volume poroso
        aux = fprop.Nj / fprop.Csi_j
        residuo_poroso_varavei = aux.sum(axis = 1) / fprop.Vp - 1.0
        residuo_poroso_bruno = aux.sum(axis = 1) - fprop.Vp
        #residuo_poroso_bruno = fprop.Vt - fprop.Vp
        # Preferencia: Bruno

        # Residuo total
        residuo = np.zeros([ctes.n_volumes*(ctes.n_components + 1),1])
        """
        for i in range(ctes.n_volumes):
            residuo_local = np.zeros([ctes.n_components + 1, 1])
            residuo_local[0] = residuo_poroso_bruno[:,i]
            residuo_local[1:] = residuo_massa[...,i][:,np.newaxis]
            residuo[i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1)] = residuo_local
        """
        bool_inds = np.ones_like(residuo,dtype=bool)
        bool_inds[0:len(residuo):ctes.n_components+1] = False
        residuo[bool_inds] = residuo_massa.T.flatten()
        residuo[~bool_inds] = residuo_poroso_bruno[0]
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
        #data = np.array([-Fk_internal_faces, Fk_internal_faces]).flatten()
        data = np.array([Fk_internal_faces, -Fk_internal_faces]).flatten() # TESTE
        Fk_vols_total = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_components, ctes.n_volumes)).toarray()
        return Fk_vols_total

    def jacobian_calculation(self, fprop, Pot_hid, delta_t, wells):
        # Calculo da matriz jacobiana
        d_residuoporoso_dP = self.d_residuoporoso_dP(fprop)
        d_residuoporoso_dNk = self.d_residuoporoso_dNk(fprop)
        d_consvmassa_dP_P, d_consvmassa_dP_E, d_consvmassa_dP_W, d_consvmassa_dNk_P, \
            d_consvmassa_dNk_E, d_consvmassa_dNk_W = self.d_consvmassa_dP_dNk(fprop, Pot_hid, delta_t, wells)

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

            index_left = ctes.v0[:,0][ctes.v0[:,1] == i]
            index_right = ctes.v0[:,1][ctes.v0[:,0] == i]

            for k in index_right:
                J_local_E = np.zeros([ctes.n_components + 1, ctes.n_components + 1])
                bool_E = (ctes.v0[:,0] == i) & (ctes.v0[:,1] == k)
                face_index_E = np.argwhere(bool_E)

                J_local_E[1:,0] = d_consvmassa_dP_E[...,face_index_E[0,0]]
                J_local_E[1:,1:] = d_consvmassa_dNk_E[...,face_index_E[0,0]]

                jacobian[i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1), \
                         (k)*(ctes.n_components + 1) : (k+1)*(ctes.n_components + 1)] = J_local_E


            for j in index_left:
                J_local_W = np.zeros([ctes.n_components + 1, ctes.n_components + 1])
                bool_W = (ctes.v0[:,0] == j) & (ctes.v0[:,1] == i)
                face_index_W = np.argwhere(bool_W)

                J_local_W[1:,0] = d_consvmassa_dP_W[...,face_index_W[0,0]]
                J_local_W[1:,1:] = d_consvmassa_dNk_W[...,face_index_W[0,0]]

                jacobian[i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1), \
                         (j)*(ctes.n_components + 1) : (j+1)*(ctes.n_components + 1)] = J_local_W


            """
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
            """
        if any(np.isnan(jacobian).flatten()): import pdb; pdb.set_trace()

        """
        #import pdb; pdb.set_trace()
        '''------------------------ Teste ----------------------'''
        #Teste = sp.csr_matrix((ctes.n_volumes*(ctes.n_components + 1), ctes.n_volumes*(ctes.n_components + 1)))
        Teste = np.zeros([ctes.n_volumes*(ctes.n_components + 1), ctes.n_volumes*(ctes.n_components + 1)])
        J_local = np.zeros([ctes.n_components + 1, ctes.n_components + 1])
        J_local_W = np.zeros([ctes.n_components + 1, ctes.n_components + 1])
        J_local_E = np.zeros([ctes.n_components + 1, ctes.n_components + 1])
        for i in range(ctes.n_volumes):
            J_local[0,0] = d_residuoporoso_dP[i]
            J_local[0,1:] = d_residuoporoso_dNk[...,i]
            J_local[1:,0] = d_consvmassa_dP_P[...,i]
            J_local[1:,1:] = d_consvmassa_dNk_P[...,i]

            J_local_W[1:,0] = d_consvmassa_dP_W[...,i]
            J_local_W[1:,1:] = d_consvmassa_dNk_W[...,i]

            J_local_E[1:,0] = d_consvmassa_dP_E[...,i]
            J_local_E[1:,1:] = d_consvmassa_dNk_E[...,i]

            index_right = ctes.v0[:, 1][ctes.v0[:, 0] == i]
            index_left = ctes.v0[:, 0][ctes.v0[:, 1] == i]


            #import pdb; pdb.set_trace()
            Teste[i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1), \
                i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1)] = J_local

            if len(index_left) > 0:
                index_left = int(index_left)
                Teste[i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1), \
                     index_left*(ctes.n_components + 1) : (index_left+1)*(ctes.n_components + 1)] = J_local_W

            if len(index_right) > 0:
                index_right = int(index_right)
                Teste[i*(ctes.n_components + 1) : (i+1)*(ctes.n_components + 1), \
                     index_right*(ctes.n_components + 1) : (index_right+1)*(ctes.n_components + 1)] = J_local_E



            #lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
            #cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
            #data = np.array([-t0[i,:], -t0[i,:], +t0[i,:], +t0[i,:]]).flatten()

        import pdb; pdb.set_trace()
        """
        return jacobian

    def d_residuoporoso_dP(self, fprop):

        if ctes.load_k:
            x = fprop.xkj[0:ctes.Nc,0,:]
            y = fprop.xkj[0:ctes.Nc,1,:]
            Nl = fprop.Nj[0,0,:]
            Nv = fprop.Nj[0,1,:]
            if ctes.load_w:
                Nw = fprop.Nj[0,2,:]

        if not ctes.load_k and ctes.load_w:
            Nw = fprop.Nj[0,0,:]

        P = fprop.P
        T = fprop.T

        if ctes.load_k:
            dlnfildP, dlnfildnij, dZldP_parcial, dZldnij_parcial, Zl = \
                    self.EOS.get_phase_derivatives(P, T, x, Nl, np.ones(ctes.n_volumes))
            dlnfivdP, dlnfivdnij, dZvdP_parcial, dZvdnij_parcial, Zv = \
                    self.EOS.get_phase_derivatives(P, T, y, Nv, np.zeros(ctes.n_volumes))
            dnldP, dnvdP, dnldNk, dnvdNk, dnildP, dnivdP, dnildNk, dnivdNk = \
                    self.EOS.dnij_dNk_dP(dlnfildP, dlnfivdP, dlnfildnij, dlnfivdnij, Nl, Nv)
        dnwdP = np.zeros(ctes.n_volumes)

        if ctes.load_w: dVw_dP = -ctes.Cw * fprop.Csi_W0 / (fprop.Csi_W ** 2)
        #dVw_dP = -ctes.Cw * fprop.Csi_W0 / (fprop.Csi_W ** 2)
        if ctes.load_k:
            dVl_dP = ctes.R * T * (P * dZldP_parcial - Zl) / (P**2)
            dVv_dP = ctes.R * T * (P * dZvdP_parcial - Zv) / (P**2)
        dporosity_dP = ctes.porosity * ctes.Cf

        if ctes.load_k and ctes.load_w:
            dRv_dP = dnldP*(1/fprop.Csi_j[0,0]) + Nl*dVl_dP + dnvdP*(1/fprop.Csi_j[0,1]) +\
                Nv*dVv_dP + dnwdP*(1/fprop.Csi_j[0,2]) + Nw*dVw_dP - ctes.Vbulk*dporosity_dP
        if ctes.load_k and not ctes.load_w:
            dRv_dP = dnldP*(1/fprop.Csi_j[0,0]) + Nl*dVl_dP + dnvdP*(1/fprop.Csi_j[0,1]) +\
                Nv*dVv_dP - ctes.Vbulk*dporosity_dP
        if not ctes.load_k and ctes.load_w:
            dRv_dP = dnwdP*(1/fprop.Csi_j[0,0]) + Nw*dVw_dP - ctes.Vbulk*dporosity_dP

        return dRv_dP


    def d_residuoporoso_dNk(self, fprop):
        # Essa funcao passará por mudanças para a inclusão do térmico. Ver item C.1 da Tese do Varavei
        dRv_dNk = np.zeros_like(fprop.Nk)
        if ctes.load_w: dRv_dNk[-1] = (1 / fprop.Csi_W) # agua

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
        if ctes.load_k:
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

        if ctes.load_k and ctes.load_w:
            dRv_dNk[:-1] = dnldNk*(1/fprop.Csi_j[0,0]) + Nl*dVl_dNk + \
                    dnvdNk*(1/fprop.Csi_j[0,1]) + Nv*dVv_dNk
        if ctes.load_k and not ctes.load_w:
            dRv_dNk = dnldNk*(1/fprop.Csi_j[0,0]) + Nl*dVl_dNk + \
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
        dx_dNk[:,:,fprop.L==0] = 0.0
        dy_dNk = (1 / fprop.Nj[0,1]) * (dnivdNk - y[:,np.newaxis]*dnvdNk)
        dy_dNk[:,:,fprop.V==0] = 0.0
        return dx_dNk, dy_dNk

    def dZ_dNk(self, dZl_dx, dZv_dy, dx_dNk, dy_dNk):

        self.dZldNk = (dZl_dx[:,np.newaxis] * dx_dNk[np.newaxis,:,:]).sum(axis=1)
        self.dZvdNk = (dZv_dy[:,np.newaxis] * dy_dNk[np.newaxis,:,:]).sum(axis=1)
        #return dZodNk, dZvdNk

    def d_consvmassa_dP_dNk(self, fprop, Pot_hid, delta_t, wells):
        dCsi_j_dP, dCsi_j_dnij = self.dCsi_j_dP_dnij(fprop)
        if ctes.load_k:
            dx_dP, dy_dP, dnildP, dnivdP, dnldP, dnvdP, dnildNk, dnivdNk, dnldNk, dnvdNk = \
                    self.EOS.dxkj_dnij_dP(fprop)
            NkT = np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0) # quando incluir a agua no termico, vai mudar esse calculo
            dx_dNk, dy_dNk = self.dxij_dNk(fprop, NkT, dnildNk, dnivdNk, dnldNk, dnvdNk) #[[dx1/dN1, dx1/dN2, ..., dx1/dNnc], [dx2/dN1, dx2/dN2, ..., dx2/dNnc], ..., [dxnc/dN1, dxnc/dN2, ..., dxnc/dNnc]] para cada bloco da malha
            drho_dP, drho_dNk = self.drho_dP_dNk(fprop, dx_dP, dy_dP, dCsi_j_dP, dx_dNk, dy_dNk)
        if not ctes.load_k and ctes.load_w:
            dx_dP, dy_dP, dx_dNk, dy_dNk = [], [], [], []
            drho_dP, drho_dNk = self.drho_dP_dNk(fprop, dx_dP, dy_dP, dCsi_j_dP, dx_dNk, dy_dNk)

        if ctes.load_k: dx_dnij, dy_dnij = self.dx_dy_dnij(fprop)

        dxij_dP = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        if ctes.load_k:
            dxij_dP[:,0] = dx_dP
            dxij_dP[:,1] = dy_dP

        dkrs_dSj, kro, krg, krw = self.dkrs_dSj(fprop)
        if ctes.load_w:
            relative_permeability = np.array([kro, krg, krw])
        else:
            relative_permeability = np.array([kro, krg])

        if ctes.load_k:
            dSj_dP, dSj_dNk = self.dSj_dP_dNk(fprop, dCsi_j_dP, dnldP, dnvdP, dnldNk, dnvdNk)
        if not ctes.load_k and ctes.load_w:
            dnldP, dnvdP, dnldNk, dnvdNk = [], [], [], []
            dSj_dP, dSj_dNk = self.dSj_dP_dNk(fprop, dCsi_j_dP, dnldP, dnvdP, dnldNk, dnvdNk)

        if data_loaded['compositional_data']['relative_permeability'] == 'StoneII':
            dkrj_dP, dkrj_dNk = self.dkrj_dP_dNk_stone(dkrs_dSj, dSj_dP, dSj_dNk, fprop, relative_permeability)
        else:
            dkrj_dP, dkrj_dNk = self.dkrj_dP_dNk_corey(dkrs_dSj, dSj_dP, dSj_dNk)

        if ctes.load_k:
            phase_viscosity = self.phase_viscosity_class(fprop, fprop.Csi_j.copy())
            dmi_dP, dmi_dNk = phase_viscosity.derivative_phase_viscosity_dP_dNk(fprop, dx_dP, dy_dP, dCsi_j_dP, dx_dnij, dy_dnij, dCsi_j_dnij, dnildP, dnivdP, dnildNk, dnivdNk)
            'Problema do BL a viscosidade é constante. derivada fica igual a 0'
            #dmi_dP = np.zeros_like(dmi_dP)
            #dmi_dNk = np.zeros_like(dmi_dNk)
        if not ctes.load_k and ctes.load_w:
            dmi_dP = np.zeros([1, ctes.n_volumes])
            dmi_dNk = np.zeros([1, ctes.n_volumes])

        # Derivada em funcao de P
        if not ctes.load_k and ctes.load_w: relative_permeability = krw
        dmobilities_dP = (fprop.mis*dkrj_dP - relative_permeability*dmi_dP)/(fprop.mis**2)
        aux = fprop.mobilities * fprop.Csi_j * dxij_dP
        dtransmissibility_dP = aux + fprop.xkj * fprop.Csi_j * dmobilities_dP \
                + fprop.xkj * fprop.mobilities * dCsi_j_dP

        # Derivada em funcao de Nk
        dmobilities_dNk = (fprop.mis * dkrj_dNk - \
            relative_permeability[np.newaxis] * dmi_dNk)/(fprop.mis**2)

        dtransmissibility_dNk = np.zeros([ctes.n_components, ctes.n_components, ctes.n_phases, ctes.n_volumes])
        aux = np.zeros([ctes.n_components, ctes.n_components, ctes.n_volumes])
        if ctes.load_k:
            # Só referente a fase óleo:
            aux[0:ctes.Nc, 0:ctes.Nc] = dx_dNk * fprop.Csi_j[0,0] * fprop.mobilities[0,0]
            dtransmissibility_dNk[:,:,0,:] = aux + \
                fprop.xkj[:,0][:,np.newaxis] * self.dCsi_j_dNk[:,0] * fprop.mobilities[0,0] +\
                fprop.xkj[:,0][:,np.newaxis] * fprop.Csi_j[0,0] * dmobilities_dNk[:,0]
            # Só referente a fase gas:
            aux[0:ctes.Nc, 0:ctes.Nc] = dy_dNk * fprop.Csi_j[0,1] * fprop.mobilities[0,1]
            dtransmissibility_dNk[:,:,1,:] = aux + \
                fprop.xkj[:,1][:,np.newaxis] * self.dCsi_j_dNk[:,1] * fprop.mobilities[0,1] +\
                fprop.xkj[:,1][:,np.newaxis] * fprop.Csi_j[0,1] * dmobilities_dNk[:,1]
            # Só referente a fase agua:
            if ctes.load_w:
                aux[0:ctes.Nc, 0:ctes.Nc] = 0 * fprop.Csi_j[0,2] * fprop.mobilities[0,2]
                dtransmissibility_dNk[:,:,2,:] = aux + \
                    fprop.xkj[:,2][:,np.newaxis] * self.dCsi_j_dNk[:,2] * fprop.mobilities[0,2] +\
                    fprop.xkj[:,2][:,np.newaxis] * fprop.Csi_j[0,2] * dmobilities_dNk[:,2]
        # dtransmissibility_dNk.shape = (Nc1, Nc2, ..., Nci); (dNc1, dNc2, ..., dNci); fases; blocos
        if not ctes.load_k and ctes.load_w:
            aux[0,0] = 0 * fprop.Csi_j[0,0] * fprop.mobilities[0,0]
            dtransmissibility_dNk[:,:,0,:] = aux + \
                fprop.xkj[:,0][:,np.newaxis] * self.dCsi_j_dNk[:,0] * fprop.mobilities[0,0] +\
                fprop.xkj[:,0][:,np.newaxis] * fprop.Csi_j[0,0] * dmobilities_dNk[:,0]

        self.dPot_hid_dP(drho_dP, drho_dNk, Pot_hid, fprop)
        transmissibility = self.transmissibility_FI(fprop)

        '''Logica testada para 1D'''
        dtransmissibility_faceup_dP_P = dtransmissibility_dP[:,:,ctes.v0[:,0]]
        dtransmissibility_facedown_dP_P = dtransmissibility_dP[:,:,ctes.v0[:,1]]
        bool = Pot_hid[:,ctes.v0[:,0]] <= Pot_hid[:,ctes.v0[:,1]]
        dtransmissibility_faceup_dP_P[:,bool] = 0.0
        dtransmissibility_facedown_dP_P[:,~bool] = 0.0

        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        d_Pot_hid = Pot_hidj_up - Pot_hidj - ctes.g * fprop.rho_j_internal_faces * (z_up - z)

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
            data = np.array([dfluxo_faceup_dP_P[i], - dfluxo_facedown_dP_P[i]]).flatten()
            dflux_dP_P[i] = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()

        dwell_dP_P = 0.0 # Caso use modelo de poço, aqui muda
        # Derivada da eq. de cons. da massa em relação a pressão do bloco (diagonal principal)
        d_consvmassa_dP_P = - delta_t * (dflux_dP_P.sum(axis=1) + dwell_dP_P)

        # Blocos vizinhos
        dtransmissibility_faceup_dP_E = dtransmissibility_dP[:,:,ctes.v0[:,1]]
        dtransmissibility_facedown_dP_W = dtransmissibility_dP[:,:,ctes.v0[:,0]]
        bool = Pot_hid[:,ctes.v0[:,0]] <= Pot_hid[:,ctes.v0[:,1]]
        dtransmissibility_faceup_dP_E[:,~bool] = 0.0
        dtransmissibility_facedown_dP_W[:,bool] = 0.0

        dfluxo_faceup_dP_E = dtransmissibility_faceup_dP_E * \
                ctes.pretransmissibility_internal_faces * d_Pot_hid + \
                transmissibility * self.dPot_hid_faceup_dP_E
        dfluxo_facedown_dP_W = dtransmissibility_facedown_dP_W * \
                ctes.pretransmissibility_internal_faces * d_Pot_hid + \
                transmissibility * self.dPot_hid_facedown_dP_W
        dfluxo_facedown_dP_E = np.zeros_like(dfluxo_faceup_dP_E)
        dfluxo_faceup_dP_W = np.zeros_like(dfluxo_facedown_dP_W)

        # Derivada do termo de fluxo em relação a pressão - blocos vizinhos (E e W)
        'Fazer esse balanço é errado'
        """
        dflux_dP_E = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        dflux_dP_W = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])

        for i in range(ctes.n_components):
            data_E = np.array([dfluxo_faceup_dP_E[i], - dfluxo_facedown_dP_E[i]]).flatten()
            data_W = np.array([dfluxo_faceup_dP_W[i], - dfluxo_facedown_dP_W[i]]).flatten()
            dflux_dP_E[i] = sp.csc_matrix((data_E, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()
            dflux_dP_W[i] = sp.csc_matrix((data_W, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()
        """
        dflux_dP_E = dfluxo_faceup_dP_E - dfluxo_facedown_dP_E
        dflux_dP_W = dfluxo_faceup_dP_W - dfluxo_facedown_dP_W

        # Derivada da eq. de cons. da massa em relação a pressão do bloco vizinho (fora da diagonal principal)
        dwell_dP_E = 0.0 # Caso use modelo de poço, aqui muda
        dwell_dP_W = 0.0 # Caso use modelo de poço, aqui muda
        d_consvmassa_dP_E = - delta_t * (dflux_dP_E.sum(axis=1) + dwell_dP_E)
        d_consvmassa_dP_W = - delta_t * (dflux_dP_W.sum(axis=1) + dwell_dP_P)

        '''Logica testada para 1D'''
        # Derivada em função da composição - diagonal principal
        dtransmissibility_faceup_dNk_P = dtransmissibility_dNk[:,:,:,ctes.v0[:,0]]
        dtransmissibility_facedown_dNk_P = dtransmissibility_dNk[:,:,:,ctes.v0[:,1]]
        bool = Pot_hid[:,ctes.v0[:,0]] <= Pot_hid[:,ctes.v0[:,1]]
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
        dwell_dNk_P = 0.0 # Caso use modelo de poço, aqui muda
        dacumulation_dNk = np.identity(ctes.n_components)[...,np.newaxis] * np.ones([ctes.n_volumes])

        d_consvmassa_dNk_P = dacumulation_dNk - delta_t * (dflux_dNk_P.sum(axis=2) + dwell_dNk_P)
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

        # Derivada do termo de fluxo em relação a composição - blocos vizinhos (E e W)
        'Fazer esse balanço é errado'
        """
        dflux_dNk_E = np.zeros([ctes.n_components, ctes.n_components, ctes.n_phases, ctes.n_volumes])
        dflux_dNk_W = np.zeros([ctes.n_components, ctes.n_components, ctes.n_phases, ctes.n_volumes])

        for i in range(ctes.n_components):
            for k in range(ctes.n_components):
                data_E = np.array([dfluxo_faceup_dNk_E[i,k], - dfluxo_facedown_dNk_E[i,k]]).flatten()
                data_W = np.array([dfluxo_faceup_dNk_W[i,k], - dfluxo_facedown_dNk_W[i,k]]).flatten()
                dflux_dNk_E[i,k] = sp.csc_matrix((data_E, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()
                dflux_dNk_W[i,k] = sp.csc_matrix((data_W, (lines, cols)), shape = (ctes.n_phases, ctes.n_volumes)).toarray()
        """
        dflux_dNk_E = dfluxo_faceup_dNk_E - dfluxo_facedown_dNk_E
        dflux_dNk_W = dfluxo_faceup_dNk_W - dfluxo_facedown_dNk_W

        # Derivada da eq. de cons. da massa em relação a pressão do bloco vizinho (fora da diagonal principal)
        dwell_dNk_E = 0.0 # Caso use modelo de poço, aqui muda
        dwell_dNk_W = 0.0 # Caso use modelo de poço, aqui muda
        d_consvmassa_dNk_E = - delta_t * (dflux_dNk_E.sum(axis=2) + dwell_dNk_E)
        d_consvmassa_dNk_W = - delta_t * (dflux_dNk_W.sum(axis=2) + dwell_dNk_W)


        return d_consvmassa_dP_P, d_consvmassa_dP_E, d_consvmassa_dP_W, d_consvmassa_dNk_P, \
                d_consvmassa_dNk_E, d_consvmassa_dNk_W

    def dCsi_j_dP_dnij(self, fprop):
        if ctes.load_k:
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

        # Minhas deduções
        dCsi_j_dP = np.zeros_like(fprop.Csi_j)

        if ctes.load_k:
            dCsi_j_dP[0,0][fprop.L!=0] = (Zl[fprop.L!=0] - \
                P[fprop.L!=0]*dZldP_parcial[fprop.L!=0]) / (ctes.R*T*(Zl[fprop.L!=0]**2))
            dCsi_j_dP[0,1][fprop.V!=0] = (Zv[fprop.V!=0] - \
                P[fprop.V!=0]*dZvdP_parcial[fprop.V!=0]) / (ctes.R*T*(Zv[fprop.V!=0]**2))
        if ctes.load_k and ctes.load_w:
            dCsi_j_dP[0,2][ctes.load_w] = fprop.Csi_W0[ctes.load_w] * ctes.Cw
        if not ctes.load_k and ctes.load_w:
            dCsi_j_dP[0,0] = fprop.Csi_W0 * ctes.Cw

        """
        # Formulação do Bruno - resultados diferentes!!
        dCsi_j_dP = np.zeros_like(fprop.Csi_j)
        dCsi_j_dP[:,0,fprop.L!=0] = - (fprop.Csi_j[0,0,fprop.L!=0] ** 2) * ((ctes.R*T)/P[fprop.L!=0]) * \
                (np.sum(dZldnij_parcial[:,:,fprop.L!=0] * dnildP[:,fprop.L!=0], axis = 1) + \
                 (1/P[fprop.L!=0])*(P[fprop.L!=0]*dZldP_parcial[fprop.L!=0] - Zl[fprop.L!=0]))
        dCsi_j_dP[:,1,fprop.V!=0] = - (fprop.Csi_j[0,1,fprop.V!=0] ** 2) * ((ctes.R*T)/P[fprop.V!=0]) * \
                (np.sum(dZvdnij_parcial[:,:,fprop.V!=0] * dnivdP[:,fprop.V!=0], axis = 1) + \
                 (1/P[fprop.V!=0])*(P[fprop.V!=0]*dZvdP_parcial[fprop.V!=0] - Zv[fprop.V!=0]))
        if ctes.load_w: dCsi_j_dP[0,2][ctes.load_w] = fprop.Csi_W0[ctes.load_w] * ctes.Cw
        """

        if ctes.load_k:
            dCsi_j_dnij = np.zeros([ctes.Nc, ctes.n_phases, ctes.n_volumes])
            dCsi_j_dnij[:,0][:,fprop.L!=0] = P[fprop.L!=0] * dZldnij_parcial[:,:,fprop.L!=0] / (ctes.R*T*(Zl[fprop.L!=0]**2))
            dCsi_j_dnij[:,1][:,fprop.V!=0] = P[fprop.V!=0] * dZvdnij_parcial[:,:,fprop.V!=0] / (ctes.R*T*(Zv[fprop.V!=0]**2))
        if not ctes.load_k and ctes.load_w:
            dCsi_j_dnij = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])

        self.dCsi_j_dNk = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        if ctes.load_k:
            self.dCsi_j_dNk[0:ctes.Nc,0][:,fprop.L!=0] = \
                -P[fprop.L!=0]/(ctes.R*T*(Zl[fprop.L!=0]**2)) * dZldNk[:,fprop.L!=0]
            self.dCsi_j_dNk[0:ctes.Nc,1][:,fprop.V!=0] = \
                -P[fprop.V!=0]/(ctes.R*T*(Zv[fprop.V!=0]**2)) * dZvdNk[:,fprop.V!=0]


        """
        # Formulação do Bruno - mesmo resultado
        self.dCsi_j_dNk = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        d_molar_volume_dnij_oil = ((ctes.R*T)/P) * dZldnij_parcial
        d_molar_volume_dnij_gas = ((ctes.R*T)/P) * dZvdnij_parcial
        self.dCsi_j_dNk[0:ctes.Nc,0] = -(fprop.Csi_j[0,0] ** 2) * \
                            np.sum(d_molar_volume_dnij_oil[:,:,np.newaxis] * dnildNk, axis=1)
        self.dCsi_j_dNk[0:ctes.Nc,1] = -(fprop.Csi_j[0,1] ** 2) * \
                            np.sum(d_molar_volume_dnij_gas[:,:,np.newaxis] * dnivdNk, axis=1)
        """
        return dCsi_j_dP, dCsi_j_dnij

    def drho_dP_dNk(self, fprop, dx_dP, dy_dP, dCsi_j_dP, dx_dNk, dy_dNk):
        if ctes.load_k:
            x = fprop.xkj[0:ctes.Nc,0,:]
            y = fprop.xkj[0:ctes.Nc,1,:]

        drho_dP = np.zeros_like(fprop.rho_j)
        if ctes.load_k:
            drho_dP[0,0] = dCsi_j_dP[0,0] * (x * ctes.Mw[:, np.newaxis]).sum(axis=0) + \
                fprop.Csi_j[0,0] * (ctes.Mw[:, np.newaxis] * dx_dP[0:ctes.Nc]).sum(axis=0)
            drho_dP[0,1] = dCsi_j_dP[0,1] * (y * ctes.Mw[:, np.newaxis]).sum(axis=0) + \
                fprop.Csi_j[0,1] * (ctes.Mw[:, np.newaxis] * dy_dP[0:ctes.Nc]).sum(axis=0)
            if ctes.load_w: drho_dP[0,2] = ctes.Mw_w * dCsi_j_dP[0,2]
        if not ctes.load_k and ctes.load_w: drho_dP[0,0] = ctes.Mw_w * dCsi_j_dP[0,0]


        drho_dNk = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        if ctes.load_k:
            dxkj_dNk = np.zeros([ctes.n_components, ctes.n_components, ctes.n_phases, ctes.n_volumes])
            dxkj_dNk[0:ctes.Nc, 0:ctes.Nc, 0, :] = dx_dNk
            dxkj_dNk[0:ctes.Nc, 0:ctes.Nc, 1, :] = dy_dNk
            if ctes.load_w:
                Mw = np.append(ctes.Mw, ctes.Mw_w)
            else:
                Mw = ctes.Mw

            aux_oil = (Mw * dxkj_dNk[:,:,0].T).sum(axis = 2)
            aux_gas = (Mw * dxkj_dNk[:,:,1].T).sum(axis = 2)

            drho_dNk[:,0] = self.dCsi_j_dNk[:,0] * (x * ctes.Mw[:,np.newaxis]).sum(axis=0) + \
                        fprop.Csi_j[0,0] * aux_oil.T
            drho_dNk[:,1] = self.dCsi_j_dNk[:,1] * (y * ctes.Mw[:,np.newaxis]).sum(axis=0) + \
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
        if not ctes.load_w:
            krw = np.zeros_like(kro)
            krj = np.append(krj, krw[np.newaxis, np.newaxis, :], axis = 1)

        dkrsdSj = self.relative_permeability.dkrs_dSj(krj, saturations)
        if not ctes.load_w: dkrsdSj[:,:,2,:] = 0.0

        return dkrsdSj, kro, krg, krw

    def dSj_dP_dNk(self, fprop, dCsi_j_dP, dnldP, dnvdP, dnldNk, dnvdNk):

        if ctes.load_w:
            dSw_dP = (- 1 / ((fprop.Vp * fprop.Csi_W)**2)) * (fprop.Nj[0,-1]*fprop.Vp * \
                dCsi_j_dP[0,-1] + fprop.Nj[0,-1]*fprop.Csi_W*ctes.Vbulk*ctes.porosity*ctes.Cf)
        else:
            dSw_dP = np.zeros(ctes.n_volumes)

        if ctes.load_k:
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
        if ctes.load_k and ctes.load_w:
            NjT = np.sum(fprop.Nj[:,0:ctes.n_phases-1,:], axis = 1)
        else:
            NjT = np.sum(fprop.Nj, axis = 1)

        if ctes.load_k:
            dL_dP = (1 / NjT) * (dnldP - fprop.L*(dnldP + dnvdP))
            dV_dP = (1 / NjT) * (dnvdP - fprop.V*(dnldP + dnvdP))
            dL_dNk = (1 / NjT) * (dnldNk - fprop.L*(dnldNk + dnvdNk))
            dV_dNk = (1 / NjT) * (dnvdNk - fprop.V*(dnldNk + dnvdNk))

            dLCsi_dP = (1/fprop.Csi_j[0,0])*(dL_dP - (fprop.L/fprop.Csi_j[0,0])*dCsi_j_dP[0,0])
            dVCsi_dP = (1/fprop.Csi_j[0,1])*(dV_dP - (fprop.V/fprop.Csi_j[0,1])*dCsi_j_dP[0,1])
            dLCsi_dNk = (1/fprop.Csi_j[0,0])*(dL_dNk - (fprop.L/fprop.Csi_j[0,0])*self.dCsi_j_dNk[0:ctes.Nc,0])
            dVCsi_dNk = (1/fprop.Csi_j[0,1])*(dV_dNk - (fprop.V/fprop.Csi_j[0,1])*self.dCsi_j_dNk[0:ctes.Nc,1])

            dSo_dP = (1 - fprop.Sw) * ((dLCsi_dP * (fprop.L / fprop.Csi_j[0,0] + \
                fprop.V / fprop.Csi_j[0,1]) - (fprop.L / fprop.Csi_j[0,0])*(dLCsi_dP + dVCsi_dP)) / \
                ((fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1])**2)) - \
                ((fprop.L / fprop.Csi_j[0,0]) / (fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1]))*dSw_dP
            dSv_dP = (1 - fprop.Sw) * ((dVCsi_dP * (fprop.L / fprop.Csi_j[0,0] + \
                fprop.V / fprop.Csi_j[0,1]) - (fprop.V / fprop.Csi_j[0,1])*(dLCsi_dP + dVCsi_dP)) / \
                ((fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1])**2)) - \
                ((fprop.V / fprop.Csi_j[0,1]) / (fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1]))*dSw_dP

        dSj_dP = np.zeros_like(fprop.Csi_j)
        if ctes.load_k:
            dSj_dP[0,0] = dSo_dP
            dSj_dP[0,1] = dSv_dP
            if ctes.load_w: dSj_dP[0,2] = dSw_dP
        if not ctes.load_k and ctes.load_w: dSj_dP = dSw_dP


        dSj_dNk = np.zeros([ctes.n_components, ctes.n_phases, ctes.n_volumes])
        if ctes.load_w:
            dSj_dNk[-1,-1] = (1 / (fprop.Csi_j[:,-1] * fprop.Vp)) - \
                            (fprop.Nj[0,-1]/(fprop.Vp*(fprop.Csi_j[:,-1]**2)))*self.dCsi_j_dNk[-1,-1]
        #else:
            #dSw_dNk = np.zeros([ctes.n_components, 1, ctes.n_volumes])
            #dSj_dNk = np.append(dSj_dNk, dSw_dNk, axis = 1)
        if ctes.load_k:
            dSj_dNk[0:ctes.Nc,0] = (1 - fprop.Sw) * ((dLCsi_dNk * (fprop.L / fprop.Csi_j[0,0] + \
                fprop.V / fprop.Csi_j[0,1]) - (fprop.L / fprop.Csi_j[0,0])*(dLCsi_dNk + dVCsi_dNk)) / \
                ((fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1])**2))

            dSj_dNk[0:ctes.Nc,1] = (1 - fprop.Sw) * ((dVCsi_dNk * (fprop.L / fprop.Csi_j[0,0] + \
                fprop.V / fprop.Csi_j[0,1]) - (fprop.V / fprop.Csi_j[0,1])*(dLCsi_dNk + dVCsi_dNk)) / \
                ((fprop.L / fprop.Csi_j[0,0] + fprop.V / fprop.Csi_j[0,1])**2))

        return dSj_dP, dSj_dNk

    def dkrj_dP_dNk_corey(self, dkrs_dSj, dSj_dP, dSj_dNk):
        dkrj_dP = np.zeros_like(dSj_dP)
        dkrj_dNk = np.zeros_like(dSj_dNk)

        # Formulação do Varavei:
        'PARA O MODELO DE BROOKS E COREY'
        if ctes.load_k and ctes.load_w:
            dkrj_dP[0,0] = dkrs_dSj[:,0,1,:] * dSj_dP[0,1,:] + dkrs_dSj[:,0,2,:] * dSj_dP[0,2,:]
            dkrj_dP[0,1] = dkrs_dSj[:,1,1,:] * dSj_dP[0,1,:]
            dkrj_dP[0,2] = dkrs_dSj[:,2,2,:] * dSj_dP[0,2,:]

            dkrj_dNk[:,0] = dkrs_dSj[:,0,2,:] * dSj_dNk[:,2,:] + dkrs_dSj[:,0,1,:] * dSj_dNk[:,1,:]
            dkrj_dNk[:,1] = dkrs_dSj[:,1,1,:] * dSj_dNk[:,1,:]
            dkrj_dNk[:,2] = dkrs_dSj[:,2,2,:] * dSj_dNk[:,2,:]
        if ctes.load_k and not ctes.load_w:
            dkrj_dP[0,0] = dkrs_dSj[:,0,1,:] * dSj_dP[0,1,:]
            dkrj_dP[0,1] = dkrs_dSj[:,1,1,:] * dSj_dP[0,1,:]

            dkrj_dNk[:,0] = dkrs_dSj[:,0,1,:] * dSj_dNk[:,1,:]
            dkrj_dNk[:,1] = dkrs_dSj[:,1,1,:] * dSj_dNk[:,1,:]
        if not ctes.load_k and ctes.load_w:
            dkrj_dP = dkrs_dSj[:,2,2,:] * dSj_dP
            dkrj_dNk = dkrs_dSj[:,2,2,:] * dSj_dNk

        '''
        # Formulação do Bruno:
        dkrj_dP[0,0] = dkrs_dSj[0,0,:]*dSj_dP[0,0,:] + dkrs_dSj[0,1,:]*dSj_dP[0,1,:] + dkrs_dSj[0,2,:]*dSj_dP[0,2,:]
        dkrj_dP[0,1] = dkrs_dSj[1,0,:]*dSj_dP[0,0,:] + dkrs_dSj[1,1,:]*dSj_dP[0,1,:] + dkrs_dSj[1,2,:]*dSj_dP[0,2,:]
        dkrj_dP[0,2] = dkrs_dSj[2,0,:]*dSj_dP[0,0,:] + dkrs_dSj[2,1,:]*dSj_dP[0,1,:] + dkrs_dSj[2,2,:]*dSj_dP[0,2,:]
        '''
        return dkrj_dP, dkrj_dNk

    def dkrj_dP_dNk_stone(self, dkrs_dSj, dSj_dP, dSj_dNk, fprop, relative_permeability):
        dkrj_dP = np.zeros_like(dSj_dP)
        dkrj_dNk = np.zeros_like(dSj_dNk)

        # Formulação do Varavei:
        'PARA O MODELO DE STONE'
        dkrj_dP[0,1] = dkrs_dSj[np.newaxis][:,1,1,:] * dSj_dP[0,1,:]
        dkrj_dP[0,2] = dkrs_dSj[np.newaxis][:,2,2,:] * dSj_dP[0,2,:]

        dkrow_dSw = self.krow0 * self.n_ow * (((1 - fprop.Sw - self.Sorw)**(self.n_ow - 1)) / \
                    ((1 - self.Swr - self.Sorw)**self.n_ow)) * (-1)
        dkrog_dSg = self.krog0 * self.n_og * (((1 - fprop.Sg - self.Swr - self.Sorg)**(self.n_og - 1)) / \
                    ((1 - self.Swr - self.Sgr - self.Sorg)**self.n_og)) * (-1)
        dkrow_dP = dkrow_dSw * dSj_dP[0,2]
        dkrog_dP = dkrog_dSg * dSj_dP[0,1]
        krog = self.krog0 * ((1. - fprop.Sg - self.Sorg - self.Swr) / (1 - self.Swr - self.Sgr - self.Sorg)) ** self.n_og
        krow = self.krow0 * ((1 - fprop.Sw - self.Sorw) / (1 - self.Swr - self.Sorw)) ** self.n_ow

        dkrj_dP[0,0] = self.krow0 * (((1/self.krow0)*dkrow_dP + dkrj_dP[0,2])*(krog/self.krow0 + relative_permeability[1]) + \
                (krow/self.krow0 + relative_permeability[2])*((1/self.krow0)*dkrog_dP + dkrj_dP[0,1]) - \
                (dkrj_dP[0,2] + dkrj_dP[0,1]))


        dkrj_dNk[:,1] = dkrs_dSj[1,1] * dSj_dNk[:,1]
        dkrj_dNk[:,2] = dkrs_dSj[2,2] * dSj_dNk[:,2]

        dkrow_dNk = dkrow_dSw * dSj_dNk[:,2]
        dkrog_dNk = dkrog_dSg * dSj_dNk[:,1]
        dkrj_dNk[:,0] = self.krow0 * (((1/self.krow0)*dkrow_dNk + dkrj_dNk[:,2])*(krog/self.krow0 + relative_permeability[1]) + \
                (krow/self.krow0 + relative_permeability[2])*((1/self.krow0)*dkrog_dNk + dkrj_dNk[:,1]) - \
                (dkrj_dNk[:,2] + dkrj_dNk[:,1]))

        '''
        # Formulação do Bruno:
        dkrj_dP[0,0] = dkrs_dSj[0,0,:]*dSj_dP[0,0,:] + dkrs_dSj[0,1,:]*dSj_dP[0,1,:] + dkrs_dSj[0,2,:]*dSj_dP[0,2,:]
        dkrj_dP[0,1] = dkrs_dSj[1,0,:]*dSj_dP[0,0,:] + dkrs_dSj[1,1,:]*dSj_dP[0,1,:] + dkrs_dSj[1,2,:]*dSj_dP[0,2,:]
        dkrj_dP[0,2] = dkrs_dSj[2,0,:]*dSj_dP[0,0,:] + dkrs_dSj[2,1,:]*dSj_dP[0,1,:] + dkrs_dSj[2,2,:]*dSj_dP[0,2,:]
        '''
        return dkrj_dP, dkrj_dNk

    def dx_dy_dnij(self, fprop):
        x = fprop.xkj[0:ctes.Nc,0,:]
        y = fprop.xkj[0:ctes.Nc,1,:]
        I = np.identity(ctes.Nc)[...,np.newaxis] * \
            np.ones([ctes.Nc, ctes.Nc, ctes.n_volumes])

        dx_dnij = (I - x[:,np.newaxis,:]*np.ones([ctes.Nc, ctes.Nc, ctes.n_volumes])) \
                    / fprop.Nj[0,0]
        dx_dnij[:,:,fprop.L==0] = 0.0

        dy_dnij = (I - y[:,np.newaxis,:]*np.ones([ctes.Nc, ctes.Nc, ctes.n_volumes])) \
                    / fprop.Nj[0,1]
        dy_dnij[:,:,fprop.V==0] = 0.0
        return dx_dnij, dy_dnij

    def dPot_hid_dP(self, drho_dP, drho_dNk, Pot_hid, fprop):
        '''Testado para 1D'''
        ones = np.ones_like(drho_dP[:,:,ctes.v0[:,0]])
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        dVp_dP_P = ctes.Vbulk * ctes.porosity * ctes.Cf

        drho_faceup_dP_P = ((dVp_dP_P[ctes.v0[:,0]] * fprop.rho_j[:,:,ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,0]] * \
                drho_dP[:,:,ctes.v0[:,0]]) * (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]]) - \
                (fprop.Vp[ctes.v0[:,0]] * fprop.rho_j[:,:,ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]] * \
                fprop.rho_j[:,:,ctes.v0[:,1]])*dVp_dP_P[ctes.v0[:,0]]) / ((fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]]) ** 2)

        self.dPot_hid_faceup_dP_P = - ones - ctes.g * (z_up - z) * drho_faceup_dP_P
        """self.dPot_hid_faceup_dP_P = - ones - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,0]] \
            / (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) * drho_dP[:,:,ctes.v0[:,0]]"""

        drho_facedown_dP_P = ((dVp_dP_P[ctes.v0[:,1]] * fprop.rho_j[:,:,ctes.v0[:,1]] + fprop.Vp[ctes.v0[:,1]] * \
                drho_dP[:,:,ctes.v0[:,1]]) * (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]]) - \
                (fprop.Vp[ctes.v0[:,0]] * fprop.rho_j[:,:,ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]] * \
                fprop.rho_j[:,:,ctes.v0[:,1]])*dVp_dP_P[ctes.v0[:,1]]) / ((fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]]) ** 2)

        self.dPot_hid_facedown_dP_P = ones - ctes.g * (z_up - z) * drho_facedown_dP_P
        """self.dPot_hid_facedown_dP_P = ones - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,1]] \
            / (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) * drho_dP[:,:,ctes.v0[:,1]]"""

        drho_faceup_dP_E = ((dVp_dP_P[ctes.v0[:,1]] * fprop.rho_j[:,:,ctes.v0[:,1]] + fprop.Vp[ctes.v0[:,1]] * \
                drho_dP[:,:,ctes.v0[:,1]]) * (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]]) - \
                (fprop.Vp[ctes.v0[:,0]] * fprop.rho_j[:,:,ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]] * \
                fprop.rho_j[:,:,ctes.v0[:,1]])*dVp_dP_P[ctes.v0[:,1]]) / ((fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]]) ** 2)

        self.dPot_hid_faceup_dP_E = ones - ctes.g * (z_up - z) * drho_faceup_dP_E
        """self.dPot_hid_faceup_dP_E = ones - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,1]] \
            / (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) * drho_dP[:,:,ctes.v0[:,1]]"""

        drho_facedown_dP_W = ((dVp_dP_P[ctes.v0[:,0]] * fprop.rho_j[:,:,ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,0]] * \
                drho_dP[:,:,ctes.v0[:,0]]) * (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]]) - \
                (fprop.Vp[ctes.v0[:,0]] * fprop.rho_j[:,:,ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]] * \
                fprop.rho_j[:,:,ctes.v0[:,1]])*dVp_dP_P[ctes.v0[:,0]]) / ((fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]]) ** 2)

        self.dPot_hid_facedown_dP_W = - ones - ctes.g *(z_up - z) * drho_facedown_dP_W
        """self.dPot_hid_facedown_dP_W = - ones - ctes.g *(z_up - z) *(fprop.Vp[ctes.v0[:,0]] \
            / (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) * drho_dP[:,:,ctes.v0[:,0]]"""


        self.dPot_hid_faceup_dNk_P = - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,0]] \
            /(fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) * drho_dNk[:,:,ctes.v0[:,0]]

        self.dPot_hid_facedown_dNk_P = - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,1]] \
            /(fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) *drho_dNk[:,:,ctes.v0[:,1]]

        self.dPot_hid_faceup_dNk_E = - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,1]] \
            /(fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) *drho_dNk[:,:,ctes.v0[:,1]]

        self.dPot_hid_facedown_dNk_W = - ctes.g * (z_up - z) * (fprop.Vp[ctes.v0[:,0]] \
            /(fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])) *drho_dNk[:,:,ctes.v0[:,0]]


    def update_gravity_term(self, fprop):
        if any((ctes.z - ctes.z[0]) != 0):
            G = ctes.g * fprop.rho_j * ctes.z
        else:
            G = np.zeros_like(fprop.rho_j)
        return G

    def dVt_derivatives(self, fprop):
        dVjdNk = np.zeros((ctes.n_components, ctes.n_phases, ctes.n_volumes))
        dVjdP = np.empty((1, ctes.n_phases, ctes.n_volumes))
        if ctes.load_k:
            self.EOS = ctes.EOS_class(fprop.T)
            if not ctes.compressible_k:
                dVjdNk[0:ctes.Nc,0,:] = 1 / fprop.Csi_j[0,0,:]
                dVjdNk[0:ctes.Nc,1,:] =  np.zeros_like(dVjdNk[0:ctes.Nc,0])#1 / fprop.Csi_j[0,1,:]
                dVjdP[0,0,:] = np.zeros(ctes.n_volumes)
                dVjdP[0,1,:] = np.zeros_like(dVjdP[0,0,:])
            else:
                dVjdP[0,0,:], dVjdP[0,1,:], dVjdNk[0:ctes.Nc,0,:], dVjdNk[0:ctes.Nc,1,:] = \
                self.EOS.get_all_derivatives(fprop)

        if ctes.load_k and ctes.load_w:
            dVjdNk[ctes.n_components-1,2,:] = 1 / fprop.Csi_j[0,ctes.n_phases-1,:]
            dVjdP[0,2,:] = - fprop.Nk[ctes.n_components-1,:] * fprop.Csi_W0 * ctes.Cw / (fprop.Csi_W)**2

        if not ctes.load_k and ctes.load_w:
            dVjdNk = 1 / fprop.Csi_j
            dVjdP[0] = - fprop.Nk * fprop.Csi_W0 * ctes.Cw / (fprop.Csi_W)**2

        return dVjdNk, dVjdP
