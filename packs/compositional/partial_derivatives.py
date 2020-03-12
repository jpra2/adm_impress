import numpy as np
import sympy
from sympy.utilities import lambdify
from .stability_check import StabilityCheck
from .properties_calculation import PropertiesCalc

class PartialDerivatives:

    def __init__(self):
        self.n_phases = 2

    def d_dn_all_derivatives(self, fprop_block, Nphase, nkphase, kprop, l):
        delta = 0.0001
        dlnf_dn = np.zeros([fprop_block.Nc, fprop_block.Nc, self.n_phases])
        dZ_dn = np.zeros([fprop_block.Nc, self.n_phases])
        for ph in range(self.n_phases):
            if Nphase[ph] != 0:
                for i in range(0,fprop_block.Nc):
                    l_plus = np.copy(l[:,ph]); l_minus = np.copy(l[:,ph])
                    l_plus[i] = (nkphase[i,ph] + delta / 2) / Nphase[ph]
                    l_minus[i] = (nkphase[i,ph] - delta / 2) / Nphase[ph]
                    dlnf_dn[:,i,ph] = (PartialDerivatives.lnf(fprop_block, kprop, l_plus, ph)
                     - PartialDerivatives.lnf(fprop_block, kprop, l_minus, ph)) / delta
                    dZ_dn[i,ph] = (self.Z(fprop_block, kprop, l_plus, ph)
                     - self.Z(fprop_block, kprop, l_minus, ph)) / delta

        return dlnf_dn, dZ_dn

    def d_dP_all_derivatives(self, fprop_block, Nphase, kprop, l, b):
        delta = 0.0001
        dZ_dP = np.zeros(self.n_phases)
        dlnf_dP = np.zeros([fprop_block.Nc, self.n_phases])
        for ph in range(0, self.n_phases):
            if Nphase[ph] != 0:
                fprop_block.P = self.Pvolume + delta/2
                Z_plus = self.Z(fprop_block, kprop, l[:,ph], ph)
                lnf_plus = PartialDerivatives.lnf(fprop_block, kprop, l[:,ph], ph)
                fprop_block.P = self.Pvolume - delta / 2
                dZ_dP[ph] = (Z_plus - self.Z(fprop_block, kprop, l[:,ph], ph)) / delta
                dlnf_dP[:,ph] = (lnf_plus - PartialDerivatives.lnf(fprop_block, kprop, l[:,ph], ph)) / delta
        return dlnf_dP, dZ_dP

    def lnf(fprop_block, kprop, l, ph):
        lnf = fprop_block.lnphi(kprop, l, ph) + np.log(fprop_block.P * l)
        return lnf

    def Z(self, fprop_block, kprop, l, ph):
        A, B = fprop_block.coefficientsPR(kprop, l)
        return StabilityCheck.Z_PR(B, A, ph)

    def dVt_derivatives(self, fprop, fprop_block, kprop):
        Nphase_allvolumes = fprop.mole_numbers_o_and_g
        nkphase_allvolumes = fprop.component_phase_mole_numbers
        l_allvolumes = fprop.component_molar_fractions
        ksi_j = fprop.ksi_o_and_g

        n_vols = len(Nphase_allvolumes[0,0,:])

        """ Initializing dN derivatives """
        dnij_dNk = np.zeros([fprop_block.Nc, fprop_block.Nc, self.n_phases, n_vols])
        dZj_dnij = np.zeros([fprop_block.Nc, self.n_phases, n_vols])

        """ Initializing dP derivatives """
        dnij_dP = np.zeros([fprop.Nc, self.n_phases, n_vols])
        dZj_dP = np.zeros([1, self.n_phases, n_vols])
        P = np.copy(fprop.P)
        Zj = np.zeros([1,self.n_phases, n_vols])
        for b in range(n_vols):
            self.Pvolume = P[b]
            fprop_block.P = self.Pvolume
            self.y = l_allvolumes[:,0,b]
            self.x = l_allvolumes[:,1,b]
            l = np.zeros([fprop_block.Nc, self.n_phases])
            l[:,0] = self.y[0:fprop_block.Nc]; l[:,1] = self.x[0:fprop_block.Nc]
            Nphase = Nphase_allvolumes[0,:,b]
            nkphase = nkphase_allvolumes[:,:,b]
            dlnfij_dnij, dZj_dnij[:,:,b] = self.d_dn_all_derivatives(fprop_block, Nphase, nkphase, kprop, l)
            dlnfj_dP, dZj_dP[0,:,b] = self.d_dP_all_derivatives(fprop_block, Nphase, kprop, l, b)
            Zj[0,:,b] = np.array([self.Z(fprop_block, kprop, self.y, 0), self.Z(fprop_block, kprop, self.x, 1)])
            matrix = np.sum(dlnfij_dnij, axis = 2)

            dlnf_dP_vector = dlnfj_dP[:,0] - dlnfj_dP[:,1]
            dnij_dP[:,1,b] =  np.linalg.inv(matrix)@dlnf_dP_vector
            dnij_dP[:,0,b] = - dnij_dP[:,1,b]

            for k in range(0, fprop_block.Nc):
                dlnfl_dnk = dlnfij_dnij[:,k,1]
                dnij_dNk[:,k,1,b] = np.linalg.inv(matrix)@dlnfl_dnk
                dnij_dNk[:,k,0,b] = - dnij_dNk[:,k,1,b]
                dnij_dNk[k,k,0,b] = 1 - dnij_dNk[k,k,1,b]

        dvj_dnij = fprop.R * fprop.T / fprop.P * dZj_dnij
        dVj_dNk = np.sum(dnij_dNk * (1 / ksi_j + Nphase[:,np.newaxis] * dvj_dnij), axis = 0)
        dVt_dNk = np.sum(dVj_dNk, axis = 1)

        dvj_dP = fprop.R * fprop.T / fprop.P * \
                (dZj_dP - Zj / fprop.P)
        dVj_dP = np.sum(Nphase_allvolumes * dvj_dP + np.sum(dnij_dP * (1 / ksi_j +
                Nphase[:,np.newaxis] * dvj_dnij), axis = 0),axis = 0)
        dVt_dP = np.sum(dVj_dP,axis = 0)

        # dVt_dNk_ = self.dVt_deriv_all(data_impress, wells, fprop, kprop)
        # import pdb; pdb.set_trace()
        #fprop.P = P #comming back
        return dVt_dNk, dVt_dP

'''def dVt_deriv_all(self, data_impress, wells, fprop, kprop): #tentar fazer analítico
        nk_allvolumes = fprop.component_mole_numbers
        Nphase_allvolumes = fprop.mole_numbers_o_and_g
        nkphase_allvolumes = fprop.component_phase_mole_numbers
        #Vt = np.sum(fprop.phase_mole_numbers / fprop.phase_molar_densities, axis = 1).ravel()
        Vt = fprop.Vt
        delta = 0.001
        n_vols = len(nk_allvolumes[0,:])
        dVt_dNk = np.zeros([fprop.Nc,n_vols])
        for k in range(fprop.Nc):
            Nk_old = nk_allvolumes
            Nk_new = nk_allvolumes
            Nk_old[k,:] = nk_allvolumes[k,:]
            Nk_new[k,:] = nk_allvolumes[k,:] + delta

            nk_old = np.copy(nkphase_allvolumes)
            nk_new = np.copy(nkphase_allvolumes)
            nk_old[k,:,:] = nkphase_allvolumes[k,:,:]
            nk_new[k,:,:] = nkphase_allvolumes[k,:,:] + delta/2

            x_save = np.copy(fprop.x)
            y_save = np.copy(fprop.y)
            fprop.x = nk_old[0:fprop.Nc,1,:]/Nphase_allvolumes[0,1,:]
            fprop.y = nk_old[0:fprop.Nc,0,:]/Nphase_allvolumes[0,0,:]
            fprop.y[nk_old[0:fprop.Nc,0,:]==0] = 0

            z_all = Nk_old/np.sum(Nk_old,axis=0)
            for i in range(n_vols):
                z = z_all[0:fprop.Nc,i]
                P = fprop.P[i]
                x = fprop.x[:,i]
                y = fprop.y[:,i]

                fprop_block = StabilityCheck(z, P, fprop.T, fprop.R, fprop.Nc, kprop)
                fprop_block.L = fprop.L[i]
                fprop_block.V = fprop.V[i]
                fprop_block.z = z
                fprop_block.x = x
                fprop_block.y = y
                fprop_block.Mw_L, fprop_block.ksi_L, fprop_block.rho_L = fprop_block.other_properties(kprop, x)
                fprop_block.Mw_V, fprop_block.ksi_V, fprop_block.rho_V = fprop_block.other_properties(kprop, y)
                fprop.update_all_volumes(fprop_block, i)
            prop_old = PropertiesCalc(data_impress, wells, fprop)
            prop_old.run_inside_loop(data_impress, wells, fprop)
            Vt_old = fprop.Vt

            fprop.x = nk_new[0:fprop.Nc,1,:]/Nphase_allvolumes[0,1,:]
            fprop.y = nk_new[0:fprop.Nc,0,:]/Nphase_allvolumes[0,0,:]
            fprop.y[nk_old[0:fprop.Nc,0,:]==0] = 0
            z_all = Nk_new/np.sum(Nk_old,axis=0)
            for i in range(n_vols):
                z = z_all[0:fprop.Nc,i]
                P = fprop.P[i]
                x = fprop.x[:,i]
                y = fprop.y[:,i]

                fprop_block = StabilityCheck(z, P, fprop.T, fprop.R, fprop.Nc, kprop)
                fprop_block.L = fprop.L[i]
                fprop_block.V = fprop.V[i]
                fprop_block.z = z
                fprop_block.x = x
                fprop_block.y = y
                fprop_block.Mw_L, fprop_block.ksi_L, fprop_block.rho_L = fprop_block.other_properties(kprop, x)
                fprop_block.Mw_V, fprop_block.ksi_V, fprop_block.rho_V = fprop_block.other_properties(kprop, y)
                fprop.update_all_volumes(fprop_block, i)
            prop_new = PropertiesCalc(data_impress, wells, fprop)
            dVt_dNk = np.zeros([fprop.Nc,n_vols])

            prop_new.run_inside_loop(data_impress, wells, fprop)
            Vt_new = fprop.Vt
            dVt_dNk[k,:] = (Vt_new - Vt_old)/delta
            import pdb; pdb.set_trace()
            fprop.x = x; fprop.y = y
        return dVt_dNk'''


'''class PartialDerivativesSym:
    def __init__(self):
        self.n_phases = 2

    def Z_PR_sym(B, A, ph):
        Z = symbols('Z')
        a = Z**3 - (1-B)*Z**2 + (A-2*B-3*B**2)*Z-(A*B-B**2-B**3)
        coef = [1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)]
        Za = solve(a,Z)
        return Za

    def lnphi_sym(self, nkphase, Nphase, ph):
        if Nphase == 0:
            lnphi = np.zeros(self.Nc)
        else:
            l = nkphase/Nphase
            A, B = self.coefficientsPR(l)
            Z = StabilityCheck.Z_PR_sym(B, A, ph)

            lnphi1 =  self.b / self.bm * (Z[0] - 1) - sympy.log(Z[0] - B) - A / (2 * (2 ** (1/2))
                        * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[0] + (1 +
                        2 ** (1/2)) * B) / (Z[0] + (1 - 2 ** (1/2)) * B))
            lnphi2 = self.b / self.bm * (Z[1] - 1) - sympy.log(Z[1] - B) - A / (2 * (2 ** (1/2))
                        * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[1] + (1 +
                        2 ** (1/2)) * B) / (Z[1] + (1 - 2 ** (1/2)) * B))
            lnphi3 =  self.b / self.bm * (Z[2] - 1) - sympy.log(Z[2] - B) - A / (2 * (2 ** (1/2))
                        * B) * (2 * self.psi / self.aalpha - self.b / self.bm) * sympy.log((Z[2] + (1 +
                        2 ** (1/2)) * B) / (Z[2] + (1 - 2 ** (1/2)) * B))
            lnphi = np.array([lnphi1,lnphi2,lnphi3])

        return lnphi

    def run_sym(self, z, ph):
        nkphase = sympy.symarray('nkphase', len(z))
        Nphaseaa = self.Nphase[0,:,0]
        Nphase = Nphaseaa[ph]
        if Nphase!=0:
            lnphi_all = self.lnphi_sym(nkphase, Nphase, ph)
            import pdb; pdb.set_trace()
            A, B = self.coefficientsPR(z)
            Z = np.array(StabilityCheck.Z_PR_sym(B, A, ph))
            reais = np.full(len(Z), False, dtype=bool)
            for r in range(0,len(Z)):
                reais[r] = sympify(Z[r]).is_real

            Z_reais = Z[reais]
            Z_ans = min(Z_reais) * ph + max(Z_reais) * (1 - ph)
            a = np.argwhere(Z == Z_ans)

            lnphi = lnphi_all[a,:]
            # lnphi = sympy.simplify(lnphi)
            nkphase_value = z * Nphase

            dlnphi_dn = sympy.diff(lnphi,nkphase[i])
            func1 = lambdify(nkphase, dlnphi_dn,'numpy')
            dlnphi_dnk = np.array(func1(*nkphase_value))
            dlnphi_dnk_resh = dlnphi_dnk.T[:,0,0,:]

        else: dlnphi_dnk_resh = np.zeros([self.Nc,self.Nc])

        return dlnphi_dnk_resh '''
