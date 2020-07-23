from packs.tpfa.biphasic.get_relative import get_relative_permeability
from packs.tpfa.biphasic.biphasic_fluid_properties import BiphasicFluidProperties
from packs.tpfa.biphasic.biphasic_utils import BiphasicUtils
from packs.tpfa.common_files import SimulationVariables
from packs.tpfa.monophasic import TpfaMonophasic
import scipy.sparse as sp
import numpy as np
import pdb

class TpfaBiphasicCons:

    def __init__(self):
        self.relative_permeability = get_relative_permeability()
        self.properties = BiphasicFluidProperties()
        self.biphasic_utils = BiphasicUtils()
        self.simulation_variables = SimulationVariables()

    def get_empty_current_biphasic_results(self):

        dty = [('loop', np.int), ('delta_t [s]', np.float), ('simulation_time [s]', np.float),
               ('oil_production [m3/s]', np.float), ('water_production [m3/s]', np.float),
               ('t [s]', np.float), ('wor', np.float), ('vpi', np.float), ('contador_vtk', np.int)]

        # return [np.array(['loop', 'delta_t [s]', 'simulation_time [s]',
        #     'oil_production [m3/s]', 'water_production [m3/s]', 't [s]', 'wor', 'vpi', 'contador_vtk'])]

        return 0

    def get_krw_and_kro(self, saturation):

        return self.relative_permeability(saturation)

    def get_mobilities_w_o(self, krw, kro):
        return self.biphasic_utils.mobility_w_volumes(krw, self.properties.mi_w), self.biphasic_utils.mobility_o_volumes(kro, self.properties.mi_o)

    def set_initial_upwind_internal_faces(self, saturation, volumes_adj_internal_faces, injector_volumes, internal_faces, faces_adj_volumes, boundary_faces, map_internal_faces):

        upwind_w = np.full((len(internal_faces), 2), False, dtype=bool)
        delta_sat = saturation[volumes_adj_internal_faces[:,0]] - saturation[volumes_adj_internal_faces[:,1]]
        delta_sat = delta_sat.flatten()
        pos = delta_sat >= 0
        upwind_w[pos, 0] = np.full(pos.sum(), True, dtype=bool)
        pos = ~pos
        upwind_w[pos, 1] = np.full(pos.sum(), True, dtype=bool)
        upwind_o = ~upwind_w.copy()

        upwind_w, upwind_o = self.upwind_wells(injector_volumes, faces_adj_volumes, boundary_faces, map_internal_faces, volumes_adj_internal_faces, upwind_w, upwind_o)
        self.test_upwind_dup(upwind_w, upwind_o)

        return upwind_w, upwind_o

    def upwind_wells(self, injector_volumes, faces_adj_volumes, boundary_faces, map_internal_faces, volumes_adj_internal_faces, upwind_w, upwind_o):

        wells_inj = injector_volumes
        if len(wells_inj) > 0:
            set_wells_inj = set(wells_inj)
            faces = np.unique(np.concatenate(faces_adj_volumes[wells_inj]))
            faces = np.setdiff1d(faces, boundary_faces)

            ids_faces_internal = map_internal_faces[faces]
            upwind_w[ids_faces_internal] = False
            upwind_o[ids_faces_internal] = False

            v0 = volumes_adj_internal_faces[ids_faces_internal]

            for volumes, i in zip(v0, ids_faces_internal):
                if set_wells_inj & set([volumes[0]]):
                    upwind_w[i, 0] = True
                    upwind_o[i, 0] = True
                elif set_wells_inj & set([volumes[1]]):
                    upwind_w[i, 1] = True
                    upwind_o[i, 1] = True

        return upwind_w, upwind_o

    def test_upwind_dup(self, upwind_w, upwind_o):

        verif1 = upwind_w[:,0] ^ upwind_w[:,1]
        verif2 = upwind_o[:,0] ^ upwind_o[:,1]
        verif1 = ~verif1
        verif2 = ~verif2
        if verif1.sum() > 0 or verif2.sum() > 0:
            raise ValueError('Duplicidade no upwind')

    def get_transmissibility_faces(self, areas_faces, internal_faces, boundary_faces, volumes_adj_internal_faces, volumes_adj_boundary_faces, upwind_w, upwind_o, mob_w, mob_o, keq_faces):

        transmissibility = np.zeros(len(areas_faces))
        transmissibility[internal_faces] = (mob_w[volumes_adj_internal_faces[upwind_w]] + mob_o[volumes_adj_internal_faces[upwind_o]])*keq_faces[internal_faces]
        transmissibility[boundary_faces] = (mob_w[volumes_adj_boundary_faces] + mob_o[volumes_adj_boundary_faces])*keq_faces[boundary_faces]

        return transmissibility

    def get_total_velocity_internal_faces(self, pressure_volumes, internal_faces, gravity_vector, volumes_adj_internal_faces, keq_internal_faces, mobility_w_internal_faces, mobility_o_internal_faces, centroid_volumes, rho_w, rho_o):
        gradient_pressure = self.simulation_variables.gradient_pressure(pressure_volumes, centroid_volumes, volumes_adj_internal_faces)
        mob_w_int_f = mobility_w_internal_faces.reshape(len(internal_faces), 1)
        mob_o_int_f = mobility_o_internal_faces.reshape(len(internal_faces), 1)
        keq = keq_internal_faces.reshape(len(internal_faces), 1)
        mob_t = mob_w_int_f + mob_o_int_f
        k1 = gradient_pressure*(mob_t)*keq
        k2 = (rho_w*mob_w_int_f + rho_o*mob_o_int_f)*keq
        k2 = gravity_vector * k2
        resp = -(k1 + k2)
        return resp

    def get_total_flux_internal_faces(self, total_velocity_internal_faces, u_direction_internal_faces, areas_internal_faces):

        f2 = np.empty(len(areas_internal_faces))

        for i, area in enumerate(areas_internal_faces):
            f2[i] = np.dot(total_velocity_internal_faces[i],  u_direction_internal_faces[i])*area

        return f2

    def get_flux_w_and_o_internal_faces(self, total_flux_internal_faces, rho_w, rho_o, mobility_w_internal_faces, mobility_o_internal_faces, nkga):

        flux_w = np.empty(len(total_flux_internal_faces))
        flux_o = flux_w.copy()

        fw = mobility_w_internal_faces/(mobility_w_internal_faces + mobility_o_internal_faces)
        fo = 1-fw

        totf = total_flux_internal_faces

        for i in range(len(total_flux_internal_faces)):
            flux_w[i] = fw[i]*(totf[i] - fo[i]*nkga[i]*(rho_w - rho_o))
            flux_o[i] = fo[i]*(totf[i] + fw[i]*nkga[i]*(rho_w - rho_o))

        return flux_w, flux_o

    def get_velocity_w_and_o_internal_faces(self, total_velocity_internal_faces, gravity_vector, keq_internal_faces, rho_w, rho_o, mobility_w_internal_faces, mobility_o_internal_faces):

        totf = total_velocity_internal_faces

        v_w = np.empty((len(totf), 3))
        v_o = v_w.copy()

        fw = mobility_w_internal_faces/(mobility_w_internal_faces + mobility_o_internal_faces)
        fo = 1-fw

        for i, keq in enumerate(keq_internal_faces):
            v_w[i] = fw[i]*(totf[i] - fo[i]*keq*gravity_vector*(rho_w - rho_o))
            v_o[i] = fo[i]*(totf[i] + fw[i]*keq*gravity_vector*(rho_w - rho_o))

        return v_w, v_o

    def test_flux(self, velocity, u_normal_direction, areas):

        flux = np.empty(len(velocity))

        for i in range(len(velocity)):
            flux[i] = np.dot(u_normal_direction[i], velocity[i])*areas[i]

        return flux

    def get_flux_phase_volumes(self, flux_internal_faces, all_wells, mobility_phase_volumes, volumes_adj_internal_faces, volumes, total_flux_volumes):

        flux_volumes = TpfaMonophasic.get_total_flux_volumes(flux_internal_faces, volumes, volumes_adj_internal_faces)
        flux_volumes[all_wells] -= total_flux_volumes[all_wells]*mobility_phase_volumes[all_wells]

        return flux_volumes
