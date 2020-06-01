from .pressure_solver import TPFASolver
from .flux_calculation import Flux
from ..solvers.solvers_scipy.solver_sp import SolverSp
from scipy import linalg
from .update_time import delta_time
import numpy as np
from ..utils import constants as ctes

class CompositionalFVM:

    def runIMPEC(self, M, wells, fprop, kprop, delta_t):
        self.update_gravity_term(fprop)
        self.get_faces_properties_upwind(fprop, kprop)
        self.get_phase_densities_internal_faces(fprop)
        r = 0.8# enter the while loop
        while (r!=1.):
            fprop.P, total_flux_internal_faces, self.q = TPFASolver().get_pressure(M, wells, fprop, kprop, delta_t, r)
            Flux().update_flux(fprop, kprop, total_flux_internal_faces)
            # For the composition calculation the time step might be different because it treats
            #composition explicitly and this explicit models are conditionally stable - wich can
            #be based on the CFL parameter.
            delta_t_new = delta_time.update_CFL(delta_t, wells, fprop)
            r = delta_t_new/delta_t
            delta_t = delta_t_new
        self.update_composition(fprop, kprop, delta_t)

        return delta_t

    def update_gravity_term(self, fprop):
        self.G = ctes.g * fprop.phase_densities * ctes.z

    def get_faces_properties_upwind(self, fprop, kprop):
        ''' Using one-point upwind approximation '''
        Pot_hid = fprop.P + fprop.Pcap - self.G[0,:,:]
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]

        fprop.mobilities_internal_faces = np.zeros([1, kprop.n_phases, ctes.n_internal_faces])
        mobilities_vols = fprop.mobilities[:,:,ctes.v0[:,0]]
        mobilities_vols_up = fprop.mobilities[:,:,ctes.v0[:,1]]
        fprop.mobilities_internal_faces[0,Pot_hidj_up <= Pot_hidj] = mobilities_vols[0,Pot_hidj_up <= Pot_hidj]
        fprop.mobilities_internal_faces[0,Pot_hidj_up > Pot_hidj] = mobilities_vols_up[0,Pot_hidj_up > Pot_hidj]

        fprop.phase_molar_densities_internal_faces = np.zeros([1, kprop.n_phases, ctes.n_internal_faces])
        phase_molar_densities_vols = fprop.phase_molar_densities[:,:,ctes.v0[:,0]]
        phase_molar_densities_vols_up = fprop.phase_molar_densities[:,:,ctes.v0[:,1]]
        fprop.phase_molar_densities_internal_faces[0,Pot_hidj_up <= Pot_hidj] = phase_molar_densities_vols[0,Pot_hidj_up <= Pot_hidj]
        fprop.phase_molar_densities_internal_faces[0,Pot_hidj_up > Pot_hidj] = phase_molar_densities_vols_up[0,Pot_hidj_up > Pot_hidj]

        fprop.component_molar_fractions_internal_faces = np.zeros([kprop.n_components, kprop.n_phases, ctes.n_internal_faces])
        component_molar_fractions_vols = fprop.component_molar_fractions[:,:,ctes.v0[:,0]]
        component_molar_fractions_vols_up = fprop.component_molar_fractions[:,:,ctes.v0[:,1]]
        fprop.component_molar_fractions_internal_faces[:,Pot_hidj_up <= Pot_hidj] = component_molar_fractions_vols[:,Pot_hidj_up <= Pot_hidj]
        fprop.component_molar_fractions_internal_faces[:,Pot_hidj_up > Pot_hidj] = component_molar_fractions_vols_up[:,Pot_hidj_up > Pot_hidj]

    def get_phase_densities_internal_faces(self, fprop):
        fprop.phase_densities_internal_faces = (fprop.Vp[ctes.v0[:,0]] * fprop.phase_densities[:,:,ctes.v0[:,0]] +
                                                fprop.Vp[ctes.v0[:,1]] * fprop.phase_densities[:,:,ctes.v0[:,1]]) /  \
                                                (fprop.Vp[ctes.v0[:,0]] + fprop.Vp[ctes.v0[:,1]])

    def update_composition(self, fprop, kprop, delta_t):
        fprop.component_mole_numbers = fprop.component_mole_numbers + delta_t * (self.q + fprop.component_flux_vols_total)
        fprop.z = fprop.component_mole_numbers[0:kprop.Nc,:] / np.sum(fprop.component_mole_numbers[0:kprop.Nc,:], axis = 0)
