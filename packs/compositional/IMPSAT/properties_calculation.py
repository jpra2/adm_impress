from packs.directories import data_loaded
from packs.compositional.IMPEC.properties_calculation import PropertiesCalc as PropertiesCalcIMPEC

from packs.utils import relative_permeability2, phase_viscosity, capillary_pressure
from packs.utils import constants as ctes
from .. import equation_of_state
import numpy as np
import sympy.utilities.lambdify as lambdify
import math


class PropertiesCalc(PropertiesCalcIMPEC):

    def run_outside_loop(self, M, fprop, wells):
        ''' This function was created to calculate the fluid properties at t=0'''

        if ctes.load_w:
            self.update_water_properties(M, fprop)
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

        fprop.mobilities = self.update_mobilities(fprop, fprop.So, fprop.Sg, fprop.Sw,
                          fprop.Csi_j, fprop.xkj)

        self.update_capillary_pressure(fprop)
