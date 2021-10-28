from ..IMPEC.pressure_solver import TPFASolver
from ..IMPEC.flux_calculation import Flux, MUSCL, FR
from ..update_time import delta_time
import numpy as np
from packs.utils import constants as ctes
from packs.directories import data_loaded
from ..IMPEC.composition_solver import Euler, RK3
from .saturation_solver import saturation as Sat
import math
from ..IMPEC.compositionalIMPEC import CompositionalFVM as CompositionalIMPEC

class CompositionalFVM(CompositionalIMPEC):

    def __call__(self, M, wells, fprop, delta_t, t):
        G = self.update_gravity_term(fprop)
        Pot_hid = fprop.P + fprop.Pcap - G[0,:,:]
        # self.get_faces_properties_average(fprop)
        self.get_faces_properties_upwind(fprop, G)
        #self.get_faces_properties_weighted_average(fprop, G)
        self.get_phase_densities_internal_faces(fprop)

        r = 0.8 # enter the while loop
        dVjdNk, dVjdP = self.dVt_derivatives(fprop)
        psolve = TPFASolver(dVjdNk, dVjdP)

        P_old = np.copy(fprop.P)
        Nk_old = np.copy(fprop.Nk)

        while (r!=1.):
            fprop.P, total_flux_internal_faces, q = psolve.get_pressure(M, wells, fprop, P_old, delta_t)
            fprop.qk_molar = q

            So, Sg, Sw, Fk_vols_total, wave_velocity, qk, mobilities = Sat(M).implicit_solver(M, fprop,
            wells, Pot_hid, total_flux_internal_faces, dVjdNk, dVjdP, P_old, q, delta_t)

            delta_t_new = delta_time.update_CFL(delta_t, fprop, wells, Fk_vols_total, fprop.Nk, wave_velocity)
            r = delta_t_new/delta_t
            delta_t = delta_t_new

        fprop.Nk, fprop.z = Euler().update_composition(Nk_old, qk, Fk_vols_total, delta_t)

        fprop.wave_velocity = wave_velocity
        fprop.mobilities = mobilities
        fprop.total_flux_internal_faces = total_flux_internal_faces
        if any(fprop.xkj.sum(axis=0).flatten()>1+1e-10): import pdb; pdb.set_trace()
        if len(fprop.Nk[fprop.Nk<0])>0: import pdb; pdb.set_trace()
        #if fprop.P[0]<fprop.P[1]: import pdb; pdb.set_trace()
        #if fprop.z[-1,0] > 0: import pdb; pdb.set_trace()

        fprop.So = So
        fprop.Sg = Sg
        fprop.Sw = Sw
        return delta_t
