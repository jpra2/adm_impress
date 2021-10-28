from .. import directories as direc
import numpy as np
from scipy.interpolate import interp1d

class BrooksAndCorey:

    def __init__(self):
        self.Sorw = float(direc.data_loaded['compositional_data']['residual_saturations']['Sorw'])
        self.Sorg = float(direc.data_loaded['compositional_data']['residual_saturations']['Sorg'])
        self.Sgr = float(direc.data_loaded['compositional_data']['residual_saturations']['Sgr'])
        self.Swr = float(direc.data_loaded['compositional_data']['residual_saturations']['Swr'])

        self.n_w = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_w'])
        self.n_o = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_o'])
        self.n_g = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_g'])

        # End-point relative permeability data
        self.krw0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krw0'])
        self.kro0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['kro0'])
        self.krg0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krg0'])

    def relative_permeabilities(self, fprop, saturations):
        #saturations = [So,Sg,Sw]

        Sorg = np.ones_like(saturations[2]) * self.Sorg
        Sorw = np.ones_like(saturations[2]) * self.Sorw
        Swr = np.ones(saturations[2].shape) * self.Swr

        Sorw[saturations[0] < self.Sorw] = saturations[0][saturations[0] < self.Sorw]
        Sorw[saturations[0] < self.Sorw] = saturations[0][saturations[0] < self.Sorw]
        Sorg[saturations[0] < self.Sorg] = saturations[0][saturations[0] < self.Sorg]
        Swr[saturations[2] < Swr] = saturations[2][saturations[2] < Swr]

        Sor = Sorw * (1 - saturations[1] / (1 - Swr - Sorg)) + \
                    Sorg * saturations[1] / (1 - Swr - Sorg)
        
        Sor[saturations[0] < Sor] = saturations[0][saturations[0] < Sor]

        krw = self.krw0 * ((saturations[2] - Swr) / (1 - Swr - Sorw - self.Sgr)) ** self.n_w
        kro = self.kro0 * ((saturations[0] - Sor) / (1 - Swr - Sorw - self.Sgr)) ** self.n_o
        krg = self.krg0 * ((saturations[1] - self.Sgr) / (1 - Swr - Sorw - self.Sgr)) ** self.n_g

        #self.Sor = Sor; self.Swr = Swr
        '''parachor_number = np.array([71, 191, 431]) #entry parameter in grams.mole
        self.phase_molar_densities = fprop.phase_molar_densities * 10**(-6) # mole/m³ to mole/cm³
        sigma = (np.sum(parachor_number[:,np.newaxis] * (self.phase_molar_densities[:,0,:] * fprop.component_molar_fractions[0:3,0,:]
                - self.phase_molar_densities[:,1,:] * fprop.component_molar_fractions[0:3,1,:]), axis = 0)) ** (4)

        f = (sigma / 8.9)**(1/7)
        Sgr = 0.15*f
        Sorg = 0.2*f
        Sg = (saturations[1] - Sgr) / (0.8 - Sgr)
        So = (0.8 - saturations[1] - Sorg) / (0.8 - Sorg)
        krw = np.zeros(So.shape)
        kro = f * So**2 + (1 - f) * So
        krg = f * Sg**2 + (1 - f) * Sg
        kro[kro<0] = 0
        krg[krg<0] = 0'''
        return kro, krg, krw, Sor

    def __call__(self, fprop, saturations):
        return self.relative_permeabilities(fprop, saturations)

    def dkrs_dSj(self, krs, saturations):
        old_settings = np.seterr(all = 'ignore', divide = 'ignore')

        Sorg = np.ones_like(saturations[2]) * self.Sorg
        Sorw = np.ones_like(saturations[2]) * self.Sorw
        Swr = np.ones(saturations[2].shape) * self.Swr

        Sorw[saturations[0] < self.Sorw] = saturations[0][saturations[0] < self.Sorw]
        Sorw[saturations[0] < self.Sorw] = saturations[0][saturations[0] < self.Sorw]
        Sorg[saturations[0] < self.Sorg] = saturations[0][saturations[0] < self.Sorg]
        Swr[saturations[2] < Swr] = saturations[2][saturations[2] < Swr]

        Sor = Sorw * (1 - saturations[1] / (1 - Swr - Sorg)) + \
                    Sorg * saturations[1] / (1 - Swr - Sorg)

        Sor[saturations[0] < Sor] = saturations[0][saturations[0] < Sor]

        den_w = (saturations[2] - Swr)
        den_o = (saturations[0] - Sor)
        den_g = (saturations[1] - self.Sgr)
        den = np.array([den_o, den_g, den_w])

        ns = np.empty((3))
        ns = np.array([self.n_o, self.n_g, self.n_w])
        dkrj_dSj = ns[np.newaxis,:,np.newaxis] * krs / den[np.newaxis,:]

        dkrj_dSj[:,den==0] = 0
        #dkrsdSj = np.empty((3,3,len(saturations[0])))
        dkrsdSj = -dkrj_dSj[np.newaxis,...] * np.ones((3,3,len(saturations[0])))

        np.einsum('iij->ij',dkrsdSj[0,:])[...] = dkrj_dSj[0,:]
        np.seterr(**old_settings)

        return dkrsdSj.transpose(0,2,1,3)

class StoneII:

    def __init__(self):
        self.Sorw = float(direc.data_loaded['compositional_data']['residual_saturations']['Sorw'])
        self.Sorg = float(direc.data_loaded['compositional_data']['residual_saturations']['Sorg'])
        self.Sgr = float(direc.data_loaded['compositional_data']['residual_saturations']['Sgr'])
        self.Swr = float(direc.data_loaded['compositional_data']['residual_saturations']['Swr'])

        self.n_w = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_w'])
        self.n_ow = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_ow'])
        self.n_og = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_og'])
        self.n_g = float(direc.data_loaded['compositional_data']['relative_permeability_data']['n_g'])

        # End-point relative permeability data
        self.krw0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krw0'])
        self.krow0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krow0'])
        self.krog0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krog0'])
        self.krg0 = float(direc.data_loaded['compositional_data']['relative_permeability_data']['krg0'])

    def relative_permeabilities(self, saturations):
        #saturations = [So,Sg,Sw]

        Sor = self.Sorw * (1 - saturations[1] / (1 - self.Swr - self.Sorg)) + \
                    self.Sorg * saturations[1] / (1 - self.Swr - self.Sorg)

        krw = self.krw0 * ((saturations[2] - self.Swr) / (1 - self.Swr - self.Sorw)) ** self.n_w
        krg = self.krg0 * ((saturations[1] - self.Sgr) / (1 - self.Swr - self.Sorg - self.Sgr)) ** self.n_g
        krow = self.krow0 * ((1 - saturations[2] - self.Sorw) / (1 - self.Swr - self.Sorw)) ** self.n_ow
        krog = self.krog0 * ((1. - saturations[1] - self.Sorg - self.Swr) / (1 - self.Swr - self.Sgr - self.Sorg)) ** self.n_og

        krw[saturations[2] <= self.Swr] = 0
        #krw[saturations[0]<= self.Swr] = self.krw0
        krow[saturations[2]<= self.Swr] = self.krow0
        krow[saturations[0]<= self.Sorw] = 0
        krog[saturations[0]<= self.Sorg] = 0

        kro = self.krow0 * ((krow/self.krow0 + krw) * (krog/self.krow0 + krg) - (krw + krg))
        kro[kro<0] = 0
        #self.krow = krow; self.krog = krog
        #kro[saturations[2]<Swr] = self.kro0 * ((saturations[0][saturations[2]<Swr] - self.Sor) / (1 - self.Sor - self.Sgr)) ** self.n_o

        '''Sw = np.array([Swr, 0.2899, 0.3778, 0.4667, 0.5556, 0.6444, 0.7, 0.7333, 0.8222, 0.9111, 1.])
        krw_ = np.array([0., 0.0020, 0.018, 0.0607, 0.1138, 0.2809, 0.4009, 0.4855, 0.7709, 1., 1.])
        krow_ = np.array([1., 0.6769, 0.4153, 0.2178, 0.0835, 0.0123, 0., 0., 0., 0., 0.])
        f_krw = interp1d(Sw, krw_)
        f_krow = interp1d(Sw, krow_)

        Sl = np.array([Swr, 0.2899, 0.3778, 0.4667, 0.5556, 0.6444, 0.7, 0.7333, 0.8222, 0.9111, 1.])
        krl = np.array([0., 0., 0., 0.011, 0.037, 0.0878, 0.1715, 0.2963, 0.4705, 0.88, 1.])
        krg_ = np.array([1., 0.56, 0.39, 0.35, 0.2, 0.1, 0.05, 0.03, 0.01, 0.001, 0.])
        f_krl = interp1d(Sl, krl)
        f_krg = interp1d(Sl, krg_)

        Sw = saturations[2]
        Sl = saturations[0] + saturations[2]
        krg = f_krg(Sl)
        krw = f_krw(Sw)
        kro = (f_krow(Sw) + krw)*(f_krl(Sl) + krg) - (krw + krg)'''

        return kro, krg, krw, Sor

    def __call__(self, fprop, saturations):
        return self.relative_permeabilities(saturations)

    def dkro_dSj(self, dkrowdSj, dkrwdSj, dkrogdSj, dkrgdSj, krw, krg, krog, krow):
        dkrodSj = self.krow0 * ((1/self.krow0 * dkrowdSj + dkrwdSj) *
                (krog/self.krow0 + krg) + (1/self.krow0 * dkrogdSj + dkrgdSj) *
                (krow/self.krow0 + krw) - (dkrwdSj + dkrgdSj))
        return dkrodSj

    def dkrs_dSj(self, krs, saturations):
        krw = krs[0,-1,:]
        krg = krs[0,1,:]

        krow = self.krow0 * ((1 - saturations[2] - self.Sorw) / (1 - self.Swr - self.Sorw)) ** self.n_ow
        krog = self.krog0 * ((1. - saturations[1] - self.Sorg - self.Swr) / (1 - self.Swr - self.Sgr - self.Sorg)) ** self.n_og

        krow[saturations[2]<= self.Swr] = self.krow0
        krow[saturations[0]<= self.Sorw] = 0
        krog[saturations[0]<= self.Sorg] = 0

        dkrwdSw = krw * self.n_w / (saturations[2] - self.Swr)
        dkrwdSw[saturations[2] <= self.Swr] = 0
        dkrgdSg = krg * self.n_g / (saturations[1] - self.Sgr)
        dkrgdSg[(saturations[1] <= self.Sgr)] = 0
        dkrowdSw = - krow * self.n_ow / (1 - saturations[2] - self.Sorw)
        dkrwdSo = -dkrwdSw
        dkrwdSg = -dkrwdSw
        dkrgdSo = -dkrgdSg
        dkrgdSw = -dkrgdSg
        dkrowdSo = - dkrowdSw
        dkrowdSg = - dkrowdSw
        dkrogdSg = - krog * self.n_og / (1. - saturations[1] - self.Sorg - self.Swr)
        dkrogdSw = - dkrogdSg
        dkrogdSo = - dkrogdSg

        dkrodSo = self.dkro_dSj(dkrowdSo, dkrwdSo, dkrogdSo, dkrgdSo, krw, krg, krog, krow)
        dkrodSw = self.dkro_dSj(dkrowdSw, dkrwdSw, dkrogdSw, dkrgdSw, krw, krg, krog, krow)
        dkrodSg = self.dkro_dSj(dkrowdSg, dkrwdSg, dkrogdSg, dkrgdSg, krw, krg, krog, krow)

        dkrsdSj = np.empty((3,3,len(saturations[0])))

        dkrsdSj[0,:] = np.array([dkrodSo,dkrodSg,dkrodSw])
        dkrsdSj[1,:] = np.array([dkrgdSo,dkrgdSg,dkrgdSw])
        dkrsdSj[2,:] = np.array([dkrwdSo,dkrwdSg,dkrwdSw])

        return dkrsdSj
