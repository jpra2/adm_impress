from packs.directories import data_loaded
from packs.utils import relative_permeability2, phase_viscosity, capillary_pressure
from packs.utils import constants as ctes
from .. import equation_of_state
import numpy as np
import sympy.utilities.lambdify as lambdify
import math


class PropertiesCalc:
    def __init__(self):
        self.relative_permeability_class = getattr(relative_permeability2,
        data_loaded['compositional_data']['relative_permeability'])
        self.relative_permeability = self.relative_permeability_class()
        self.phase_viscosity_class = getattr(phase_viscosity,
        data_loaded['compositional_data']['phase_viscosity'])
        self.EOS_class = getattr(equation_of_state,
        data_loaded['compositional_data']['equation_of_state'])

    def run_outside_loop(self, M, fprop, wells):
        ''' This function was created to calculate the fluid properties at t=0'''

        if ctes.load_w:
            fprop.Sw = self.update_water_properties(M, fprop)
            #fprop.Sw = self.Sw

        if ctes.load_k:
            fprop.So, fprop.Sg = self.update_saturations(M.data['saturation'],
                            fprop.Csi_j, fprop.L, fprop.V)
        else:
            fprop.So = np.zeros(ctes.n_volumes)
            fprop.Sg = np.zeros(ctes.n_volumes)
        #self.So = fprop.So; self.Sg = fprop.Sg

        self.set_initial_volume(fprop)
        self.set_initial_mole_numbers(fprop)

        fprop.mobilities = self.update_mobilities(fprop, fprop.So, fprop.Sg, fprop.Sw,
                          fprop.Csi_j, fprop.xkj)

        self.update_capillary_pressure(fprop)
        if ctes.FR: fprop.Nk_SP = self.set_Nk_Pspace(fprop)

    def run_inside_loop(self, M, fprop):
        ''' Function to update fluid and reservoir properties along the \
        simulation '''

        #fprop.Vp = self.update_porous_volume(fprop.P)

        if ctes.load_w:
            self.update_water_properties(M, fprop)

        self.update_mole_numbers(fprop)
        self.update_total_volume(fprop)

        #self.So, self.Sg = self.update_saturations(M.data['saturation'],
        #                    fprop.Csi_j, fprop.L, fprop.V)
        #fprop.So = self.So; fprop.Sg = self.Sg; fprop.Sw = self.Sw
        fprop.mobilities = self.update_mobilities(fprop, fprop.So, fprop.Sg, fprop.Sw,
                          fprop.Csi_j, fprop.xkj)

        self.update_capillary_pressure(fprop)


    def update_water_properties(self, M, fprop):
        if data_loaded['compositional_data']['water_data']['mobility']:
            M.data['saturation'], fprop.Csi_W, fprop.rho_W = \
            self.update_water_saturation(fprop, fprop.Nk[-1,:], fprop.P, fprop.Vp, fprop.Csi_W0)
        Sw = M.data['saturation']
        fprop.Csi_j[0, ctes.n_phases-1,:] = fprop.Csi_W
        fprop.rho_j[0,ctes.n_phases-1,:] = fprop.rho_W
        return Sw

    def set_initial_volume(self, fprop):
        self.Vo = fprop.Vp * fprop.So
        self.Vg = fprop.Vp * fprop.Sg
        self.Vw = fprop.Vp * fprop.Sw
        fprop.Vt = self.Vo + self.Vg + self.Vw


    def set_initial_mole_numbers(self, fprop):
        if ctes.load_k:
            fprop.Nj[0,0,:] = fprop.Csi_j[0,0,:] * self.Vo
            fprop.Nj[0,1,:] = fprop.Csi_j[0,1,:] * self.Vg
        if ctes.load_w:
            fprop.Nj[0,ctes.n_phases-1,:] = fprop.Csi_W * self.Vw

        Nkj = fprop.xkj * fprop.Nj
        fprop.Nk = np.sum(Nkj, axis = 1)


    def set_Nk_Pspace(self, fprop):
        from packs.compositional import prep_FR as ctes_FR
        'Function to get Nk at the solution points at t=0 for the FR approach. \
        For now, this only works for constant Nk distribution'
        #x = Symbol('x')
        #Nk_vols_no_wells = fprop.Nk[:,ctes.vols_no_wells]
        Nk_SP = np.ones((ctes.n_components,ctes.n_volumes,ctes_FR.n_points))
        Nk_SP = Nk_SP * fprop.Nk[:,:,np.newaxis]
        return Nk_SP

    def update_porous_volume(self, P):
        #fprop.porosity = ctes.porosity * (1 + ctes.Cf * (fprop.P - ctes.Pf))
        Vp = ctes.porosity * ctes.Vbulk * (1 + ctes.Cf*(P - ctes.Pf))
        return Vp

    def update_saturations(self, Sw, Csi_j, L, V):
        if ctes.load_k:
            Sg = (1. - Sw) * (V / Csi_j[0,1,:]) / (V / Csi_j[0,1,:] + L /
                Csi_j[0,0,:])
            So = 1 - Sw - Sg
        else: So = np.zeros_like(Sw); Sg = np.zeros_like(Sw)
        return So, Sg

    '''def update_saturations(self, fprop, Nj, Csi_j, Sw):
        if ctes.load_k:
            So = Nj[0,0,:]/Csi_j[0,0,:]/fprop.Vt
            Sg = Nj[0,1,:]/Csi_j[0,1,:]/fprop.Vt
        else: So = np.zeros(ctes.n_volumes) ; Sg = np.zeros(ctes.n_volumes)
        return So, Sg'''

    def update_mole_numbers(self, fprop):
        # Esta função foi criada separada pois, quando compressivel, o volume poroso pode
        #diferir do volume total, servindo como um termo de correção de erro na equação da pressão,
        #como se sabe. A equação do outside loop ela relaciona o volume total com o poroso, de modo
        #que eles nunca vão ser diferentes e um erro vai ser propagado por toda a simulação. Todavia,
        #ele funciona para o primeiro passo de tempo uma vez que a pressão não mudou e Vp = Vt ainda.

        if ctes.load_k:
            fprop.Nj[0,0,:] = np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0) * fprop.L
            fprop.Nj[0,1,:] = np.sum(fprop.Nk[0:ctes.Nc,:], axis = 0) * fprop.V

        if ctes.load_w:
            fprop.Nj[0,ctes.n_phases-1,:] = fprop.Nk[ctes.n_components-1,:]

    def update_total_volume(self, fprop):
        fprop.Vt = np.sum(fprop.Nj / fprop.Csi_j, axis = 1).ravel()

    def update_relative_permeabilities(self, fprop, So, Sg, Sw):
        saturations = np.array([So, Sg, Sw])
        kro,krg,krw, fprop.Sor = self.relative_permeability(fprop, saturations)
        krj = np.zeros([1, ctes.n_phases, len(So)])
        if ctes.load_k:
            krj[0,0,:] = kro
            krj[0,1,:] = krg
        if ctes.load_w:
            krj[0, ctes.n_phases-1,:] = krw
        return krj

    def relative_permeability_derivative_call(self, krs, saturations):
        dkrsdSj = self.relative_permeability.dkrs_dSj(krs, saturations)
        return dkrsdSj

    def update_phase_viscosities(self, fprop, Csi_j, xkj):
        phase_viscosities = np.empty_like(Csi_j)
        if ctes.load_k:
            #phase_viscosity = self.phase_viscosity_class(fprop, Csi_j)
            phase_viscosities[0,0:2,:] = 0.02*np.ones([2,len(Csi_j[0,0,:])]) #0.02 only for BL test. for BL_Darlan use 1e-3
            #phase_viscosities[0,0:2,:] = 0.001*np.ones([2,len(Csi_j[0,0,:])]) #only for Dietz test
            #phase_viscosities[0,0:2,:] = phase_viscosity(fprop, xkj)

        if ctes.load_w:
            phase_viscosities[0,ctes.n_phases-1,:] = data_loaded['compositional_data']['water_data']['mi_W']

        return phase_viscosities

    def update_mobilities(self, fprop, So, Sg, Sw, Csi_j, xkj):
        krs = self.update_relative_permeabilities(fprop, So, Sg, Sw)
        mis = self.update_phase_viscosities(fprop, Csi_j, xkj)
        mobilities = krs / mis
        return mobilities

    def update_capillary_pressure(self, fprop):
        """ not working yet """
        #get_capillary_pressure = getattr(capillary_pressure, data_loaded['compositional_data']['capillary_pressure'])
        #get_capillary_pressure = get_capillary_pressure(data_loaded, data_impress, fprop.Csi_j, fprop.xkj)
        #Pcow, Pcog = get_capillary_pressure(data_loaded, fprop.Sw, fprop.So, fprop.Sg)

        fprop.Pcap = np.zeros([ctes.n_phases,ctes.n_volumes])
        # Pcap[0,0,:] = Pcog
        # Pcap[0,1,:] = Pcow

    def update_water_saturation(self, fprop, Nw, P, Vp, Csi_W0):
        Csi_W = Csi_W0 * (1 + ctes.Cw * (P - ctes.Pw))
        rho_W = Csi_W * ctes.Mw_w
        Sw = Nw * (1 / Csi_W) / Vp
        return Sw, Csi_W, rho_W
