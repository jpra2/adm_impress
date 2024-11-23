from packs.tpfa.biphasic.get_relative import get_relative_permeability
from packs.tpfa.biphasic.biphasic_fluid_properties import BiphasicFluidProperties
from packs.tpfa.biphasic.biphasic_utils import BiphasicUtils
from packs.tpfa.biphasic.simulation_infos import SimulationInfos
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
        self.simulation_infos = SimulationInfos()

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

    def get_transmissibility_faces(self, areas_faces, internal_faces, boundary_faces, volumes_adj_internal_faces, volumes_adj_boundary_faces, mobility_w_internal_faces, mobility_o_internal_faces, mob_w, mob_o, keq_faces):

        transmissibility = np.empty(len(areas_faces))
        transmissibility[internal_faces] = (mobility_w_internal_faces + mobility_o_internal_faces)*keq_faces[internal_faces]*areas_faces[internal_faces]
        transmissibility[boundary_faces] = (mob_w[volumes_adj_boundary_faces] + mob_o[volumes_adj_boundary_faces])*keq_faces[boundary_faces]*areas_faces[boundary_faces]

        return transmissibility

    def get_transmissibility_faces_mrst(self, areas_faces, internal_faces, boundary_faces, volumes_adj_internal_faces, volumes_adj_boundary_faces, mob_w, mob_o, keq_faces, k_faces, one_sided_total_mobility, h_faces):


        kmob = (k_faces / h_faces) * one_sided_total_mobility
        kmob = 1 / ( 1 / kmob[:,0] + 1 / kmob[:,1])


        # transmissibility = np.empty(len(areas_faces))
        # transmissibility[internal_faces] = kmob[internal_faces]*areas_faces[internal_faces]
        # transmissibility[boundary_faces] = kmob[boundary_faces]*areas_faces[boundary_faces]

        transmissibility = kmob * areas_faces

        return transmissibility

    def get_total_velocity_internal_faces_dep(self, pressure_volumes, internal_faces, gravity_vector, volumes_adj_internal_faces, keq_internal_faces, mobility_w_internal_faces, mobility_o_internal_faces, rho_w, rho_o, abs_u_normal_internal_faces, hi):
        # gradient_pressure = self.simulation_variables.gradient_pressure(pressure_volumes, centroid_volumes, volumes_adj_internal_faces)
        # pressure_direction = self.simulation_variables.pressure_direction(pressure_volumes, centroid_volumes, volumes_adj_internal_faces)
        delta_p = self.simulation_variables.delta_p(pressure_volumes, volumes_adj_internal_faces)
        delta_p = delta_p.reshape(len(delta_p), 1)
        ni = len(keq_internal_faces)
        # abs_u_normal_internal_faces2 = abs_u_normal_internal_faces.reshape(ni, 3, 1)
        # abs_u_normal_internal_faces2 = abs_u_normal_internal_faces.copy()
        pressure_direction = -delta_p*abs_u_normal_internal_faces
        mob_w_int_f = mobility_w_internal_faces.reshape(ni, 1)
        mob_o_int_f = mobility_o_internal_faces.reshape(ni, 1)
        keq = keq_internal_faces.reshape(ni, 1)
        mob_t = mob_w_int_f + mob_o_int_f
        # pdb.set_trace()
        hi2 = hi.sum(axis=1).reshape(ni, 1)

        k1 = pressure_direction*(mob_t)*keq

        k2 = (rho_w*mob_w_int_f + rho_o*mob_o_int_f)*gravity_vector*keq*hi2
        resp = k1 + k2
        return resp

    def get_total_velocity_internal_faces(self, pressure_volumes, volumes_adj_internal_faces, keq_internal_faces, mobility_w_internal_faces, mobility_o_internal_faces, abs_u_normal_internal_faces, g_total_velocity_internal_faces):
        # gradient_pressure = self.simulation_variables.gradient_pressure(pressure_volumes, centroid_volumes, volumes_adj_internal_faces)
        # pressure_direction = self.simulation_variables.pressure_direction(pressure_volumes, centroid_volumes, volumes_adj_internal_faces)
        delta_p = self.simulation_variables.delta_p(pressure_volumes, volumes_adj_internal_faces)
        delta_p = delta_p.reshape(len(delta_p), 1)
        ni = len(keq_internal_faces)
        # abs_u_normal_internal_faces2 = abs_u_normal_internal_faces.reshape(ni, 3, 1)
        # abs_u_normal_internal_faces2 = abs_u_normal_internal_faces.copy()
        pressure_direction = -delta_p*abs_u_normal_internal_faces
        mob_w_int_f = mobility_w_internal_faces.reshape(ni, 1)
        mob_o_int_f = mobility_o_internal_faces.reshape(ni, 1)
        keq = keq_internal_faces.reshape(ni, 1)
        mob_t = mob_w_int_f + mob_o_int_f
        # pdb.set_trace()

        k1 = pressure_direction*(mob_t)*keq
        resp = k1 + g_total_velocity_internal_faces
        return resp

    def get_velocity(self, total_flux_internal_faces, areas_internal_faces, abs_u_normal_internal_faces):

        r1 = total_flux_internal_faces/areas_internal_faces
        r1 = r1.reshape(len(r1), 1)
        resp = r1*abs_u_normal_internal_faces
        return resp

    def get_total_flux_internal_faces1(self, pressure_volumes, internal_faces, gravity_vector, volumes_adj_internal_faces, keq_internal_faces, mobility_w_internal_faces, mobility_o_internal_faces, rho_w, rho_o, abs_u_normal_internal_faces, hi, areas_internal_faces):
        # gradient_pressure = self.simulation_variables.gradient_pressure(pressure_volumes, centroid_volumes, volumes_adj_internal_faces)
        # pressure_direction = self.simulation_variables.pressure_direction(pressure_volumes, centroid_volumes, volumes_adj_internal_faces)
        delta_p = self.simulation_variables.delta_p(pressure_volumes, volumes_adj_internal_faces)
        delta_p = delta_p.reshape(len(delta_p), 1)
        ni = len(keq_internal_faces)
        # abs_u_normal_internal_faces2 = abs_u_normal_internal_faces.reshape(ni, 3, 1)
        # abs_u_normal_internal_faces2 = abs_u_normal_internal_faces.copy()
        pressure_direction = -delta_p*abs_u_normal_internal_faces
        mob_w_int_f = mobility_w_internal_faces.reshape(ni, 1)
        mob_o_int_f = mobility_o_internal_faces.reshape(ni, 1)
        keq = keq_internal_faces.reshape(ni, 1)
        mob_t = mob_w_int_f + mob_o_int_f
        # pdb.set_trace()
        hi2 = hi.sum(axis=1).reshape(ni, 1)

        k1 = pressure_direction*(mob_t)*keq

        k2 = (rho_w*mob_w_int_f + rho_o*mob_o_int_f)*gravity_vector*keq*hi2
        resp = ((k1 + k2)*abs_u_normal_internal_faces).sum(axis=1)
        resp = resp*areas_internal_faces
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

    def get_velocity_w_and_o_internal_faces(self, total_velocity_internal_faces, gravity_vector, keq_internal_faces, rho_w, rho_o, mobility_w_internal_faces, mobility_o_internal_faces, hi, abs_normal_internal_faces):

        totf = total_velocity_internal_faces
        ni = len(keq_internal_faces)

        mobw2 = mobility_w_internal_faces.reshape(ni, 1)
        mobo2 = mobility_o_internal_faces.reshape(ni, 1)
        keq2 = keq_internal_faces.reshape(ni, 1)
        hi2 = hi.sum(axis=1).reshape(ni, 1)
        fw2 = mobw2/(mobw2 + mobo2)
        fo2 = 1 - fw2
        k1 = keq2*gravity_vector*(rho_w - rho_o)*hi2
        v_w = fw2*(totf + mobo2*k1)
        v_o = fo2*(totf - mobw2*k1)

        return v_w, v_o

    def test_flux(self, velocity, u_normal_direction, areas):

        flux = np.empty(len(velocity))

        for i in range(len(velocity)):
            flux[i] = np.dot(u_normal_direction[i], velocity[i])*areas[i]

        return flux

    def get_flux_phase_volumes(self, flux_phase_internal_faces, all_wells, flux_frac_phase_volumes, volumes_adj_internal_faces, volumes, total_flux_volumes):

        fw_phase = flux_frac_phase_volumes
        flux_volumes = TpfaMonophasic.get_total_flux_volumes(flux_phase_internal_faces, volumes, volumes_adj_internal_faces)
        flux_volumes[all_wells] -= total_flux_volumes[all_wells]*fw_phase[all_wells]
        return flux_volumes

    def get_flux_phase_volumes_with_wells(self, flux_phase_internal_faces, all_wells, flux_frac_phase_volumes, volumes_adj_internal_faces, volumes, total_flux_volumes):

        # fw_phase = flux_frac_phase_volumes
        flux_volumes = TpfaMonophasic.get_total_flux_volumes(flux_phase_internal_faces, volumes, volumes_adj_internal_faces)
        # flux_volumes[all_wells] -= total_flux_volumes[all_wells]*fw_phase[all_wells]
        return flux_volumes

    def update_delta_t(self, flux_w_volumes, porosity, vol_volumes, total_flux_volumes, total_velocity_internal_faces, flux_frac_w_volumes, saturation, volumes_adj_internal_faces, hi, abs_u_normal_internal_faces):

        deltas_t = []
        deltas_t.append(self.update_delta_t_for_delta_sat_max(flux_w_volumes, porosity, vol_volumes))
        deltas_t.append(self.update_delta_t_dep0(total_flux_volumes, porosity, vol_volumes))
        deltas_t.append(self.update_delta_t_new(porosity, vol_volumes, total_velocity_internal_faces, flux_frac_w_volumes, saturation, volumes_adj_internal_faces, hi, abs_u_normal_internal_faces))
        deltas_t.append(self.update_delta_t2(np.linalg.norm(total_velocity_internal_faces*abs_u_normal_internal_faces, axis=1), hi)*0.8)

        delta_t = min(deltas_t)

        return delta_t

    def update_delta_t_for_delta_sat_max(self, flux_w_volumes, porosity, vol_volumes):

        delta_sat_max = self.simulation_infos.data['delta_sat_max']
        phis = porosity
        volume = vol_volumes

        test = flux_w_volumes/(volume*phis)
        try:
            test = np.absolute(test[np.absolute(test)>0]).max()
        except Exception as e:
            import pdb; pdb.set_trace()
        delta_t = delta_sat_max/test
        return delta_t

    def update_delta_t_new(self, porosity, vol_volumes, total_velocity_internal_faces, flux_frac_w_volumes, saturation, volumes_adj_internal_faces, hi, abs_u_normal_internal_faces):

        ####
        # de acordo com as faces
        ####



        lim_ds = 1e-10
        cfl = self.simulation_infos.data['cfl']

        phis = porosity
        volume = vol_volumes
        # velocity_faces = np.absolute(self.data_impress['velocity_faces'])
        # internal_faces = self.elements_lv0['internal_faces']
        fw_vol = flux_frac_w_volumes
        # saturation = self.data_impress['saturation']
        viz_int = volumes_adj_internal_faces

        dists_int = hi.sum(axis=1)
        vel_internal_faces = np.linalg.norm(total_velocity_internal_faces*abs_u_normal_internal_faces, axis=1)
        v0 = viz_int[:, 0]
        v1 = viz_int[:, 1]
        df = np.absolute(fw_vol[v1] - fw_vol[v0])
        ds = np.absolute(saturation[v1] - saturation[v0])
        ids = np.arange(len(ds))
        ids_2 = ids[ds > lim_ds]

        dists_int = dists_int[ids_2]
        vel_internal_faces = vel_internal_faces[ids_2]
        df = df[ids_2]
        ds = ds[ids_2]
        dfds = df/ds
        try:
            delta_t = (cfl*(dists_int/(vel_internal_faces*dfds))).min()
        except ValueError:
            delta_t = np.inf
        return delta_t

    def update_delta_t_dep0(self, total_flux_volumes, porosity, vol_volumes):
        ###
        ## de acordo com o fluxo nos volumes
        ###

        cfl = self.simulation_infos.data['cfl']

        flux_volumes = np.absolute(total_flux_volumes)
        phis = porosity
        volume = vol_volumes

        with np.errstate(divide='ignore'):
            delta_t = (cfl*(volume*phis)/flux_volumes).min()

        return delta_t

    def update_delta_t2(self, total_velocity_internal_faces, hi):
        cfl = self.simulation_infos.data['cfl']
        abs_velocity = np.absolute(total_velocity_internal_faces)
        dx = hi.sum(axis=1)
        with np.errstate(divide='ignore'):
            deltas_t = cfl*(dx/abs_velocity)
        delta_t = deltas_t.min()
        return delta_t

    def update_saturation(self, flux_w_volumes, porosity, vol_volumes, total_flux_volumes, total_velocity_internal_faces, flux_frac_w_volumes, saturation0, volumes_adj_internal_faces, hi, abs_u_normal_internal_faces, dt=None, loop=None):
        cont = 0
        max_loops = 100

        delta_t = self.update_delta_t(
            flux_w_volumes,
            porosity,
            vol_volumes,
            total_flux_volumes,
            total_velocity_internal_faces,
            flux_frac_w_volumes,
            saturation0,
            volumes_adj_internal_faces,
            hi,
            abs_u_normal_internal_faces
        )

        if dt is None:
            pass
        else:
            delta_t = min(dt, delta_t)

        verif = -1

        while verif != 0:
            verif, saturation = self.update_sat(
                saturation0,
                flux_w_volumes,
                vol_volumes,
                porosity,
                delta_t
            )
            if verif == 1:
                print(loop)
                import pdb; pdb.set_trace()
                delta_t = self.reduce_delta_t(delta_t)

            if cont == max_loops:
                raise ValueError('Loop maximo atingido')

            cont += 1

        return delta_t, saturation

    def update_sat(self, saturation0, flux_w_volumes, vol_volumes, porosity, delta_t):

        lim_flux_w = 0
        # lim_flux_w = 1e-8
        delta_sat_max = self.simulation_infos.data['delta_sat_max']
        Swc = self.relative_permeability.Swc
        Sor = self.relative_permeability.Sor
        deltt = 0.0001

        # import pdb; pdb.set_trace()

        saturations = saturation0.copy()
        ids = np.arange(len(saturations))

        fw_volumes = -flux_w_volumes
        # fw_volumes = flux_w_volumes

        volumes = vol_volumes
        phis = porosity

        # import pdb; pdb.set_trace()

        # ids_2 = ids[fw_volumes < 0]
        # if len(ids_2) > 0:
        #     self.data_impress['test_fw'][ids_2] = np.repeat(1.0, len(ids_2))
        #     self.data_impress.update_variables_to_mesh()
        #     self.mesh.core.print(file='test', extension='.vtk', config_input="input_cards/print_settings0.yml")
        #     import pdb; pdb.set_trace()
        #     self.data_impress['test_fw'] = np.repeat(0.0, len(self.data_impress['test_fw']))
        #     self.data_impress.update_variables_to_mesh(['test_fw'])

        ###########################
        ## teste
        test = ids[(saturations < 0) | (saturations > 1)]
        if len(test) > 0:
            raise ValueError(f'valor errado da saturacao {saturations[test]}')
        del test

        ids_var = ids[np.absolute(fw_volumes) > lim_flux_w]

        ###################
        ## teste variacao do fluxo de agua
        if len(ids_var) == 0:
            import pdb; pdb.set_trace()
        ###################

        fw_volumes = fw_volumes[ids_var]
        volumes = volumes[ids_var]
        phis = phis[ids_var]
        sats = saturations[ids_var]

        sats += (fw_volumes*delta_t)/(phis*volumes)

        delta_sat = sats - saturations[ids_var]
        ids2 = np.arange(len(delta_sat))

        #############
        ## teste variacao maxima de saturacao
        test = ids2[np.absolute(delta_sat) > delta_sat_max+0.000001]
        if len(test) > 0:
            import pdb; pdb.set_trace()
            return 1, True
        del test
        ##############

        saturations[ids_var] = sats

        # if np.allclose(saturations, saturations0):
        #     import pdb; pdb.set_trace()


        min_sat = saturations.min()
        max_sat = saturations.max()

        ## teste dos limites
        if min_sat < Swc - deltt or max_sat > 1-Sor + deltt:
            # import pdb; pdb.set_trace()
            return 1, True

        return 0, saturations

    def reduce_delta_t(self, dt):
        delta_t = dt*1/2
        print(f'\nreducing delta_t: {dt} -> {delta_t} \n')
        return delta_t

    def get_production_w_o(self, flux_total_volumes, mob_w, mob_o, wells_injector, wells_producer, delta_t, vol_volumes, porosity):

        fw_volumes = mob_w/(mob_w + mob_o)
        fo_volumes = 1 - fw_volumes
        prod_o = -flux_total_volumes[wells_producer]*fo_volumes[wells_producer]
        prod_o = prod_o.sum()*delta_t

        prod_w = -flux_total_volumes[wells_producer]*fw_volumes[wells_producer]
        prod_w = prod_w.sum()*delta_t

        wor = prod_w/prod_o

        vol_all = vol_volumes*porosity
        vol_all = vol_all.sum()

        dvpi = flux_total_volumes[wells_injector].sum()
        dvpi = abs((dvpi*delta_t)/vol_all)

        return prod_w, prod_o, wor, dvpi

    def update_upwind_phases(self, centroid_volumes, internal_faces, volumes_adj_internal_faces, saturation, flux_w_internal_faces, flux_o_internal_faces, total_flux_internal_faces, gravity=True):
        '''
            paper Starnoni
            with gravity
        '''


        k0 = 9e-3
        min_abs = np.min(np.absolute(total_flux_internal_faces))
        v0 = volumes_adj_internal_faces
        # saturation = self.data_impress['saturation']
        flux_w_faces = flux_w_internal_faces.copy()
        flux_o_faces = flux_o_internal_faces.copy()
        flux_total = total_flux_internal_faces.copy()
        testw = np.absolute(flux_w_faces) < k0
        testo = np.absolute(flux_o_faces) < k0
        test3 = np.absolute(flux_total) < k0
        flux_w_faces[testw] = 0.0
        flux_o_faces[testo] = 0.0
        flux_total[test3] = 0.0

        tw = flux_w_faces > 0
        to = flux_o_faces > 0
        ttotal = flux_total > 0

        upwind_w = np.full((len(internal_faces), 2), False, dtype=bool)
        upwind_o = upwind_w.copy()

        upt = upwind_o.copy()
        upt[:,0] = True

        if gravity is True:
            # n_zero_condition = (~test3) & (~testw) & (~testo)
            # pdb.set_trace()
            #
            # # internal_faces = self.elements_lv0['internal_faces']
            #
            # v1 = ttotal & tw & n_zero_condition
            # v2 = ttotal & to & n_zero_condition
            #
            # upwind_w[v1, 0] = True
            # upwind_o[v2, 0] = True
            #
            # v1 = ttotal & (~tw) & n_zero_condition
            # v2 = ttotal & (~to) & n_zero_condition
            #
            # upwind_w[v1, 1] = True
            # upwind_o[v2, 1] = True
            #
            # v1 = (~ttotal) & (~tw) & n_zero_condition
            # v2 = (~ttotal) & (~to) & n_zero_condition
            #
            # upwind_w[v1, 1] = True
            # upwind_o[v2, 1] = True
            #
            # v1 = (~ttotal) & (tw) & n_zero_condition
            # v2 = (~ttotal) & (to) & n_zero_condition
            #
            # upwind_w[v1, 0] = True
            # upwind_o[v2, 0] = True

            # pdb.set_trace()

            v1 = (~test3) & (testw) & to
            v2 = (~test3) & (testo) & tw

            upwind_w[v1, 0] = True
            upwind_o[v2, 0] = True

            v1 = (~test3) & (testw) & (~to)
            v2 = (~test3) & (testo) & (~tw)

            upwind_w[v1, 1] = True
            upwind_o[v2, 1] = True

            v1 = (~test3) & (~testw) & tw
            v2 = (~test3) & (~testo) & to

            upwind_w[v1, 0] = True
            upwind_o[v2, 0] = True

            v1 = (~test3) & (~testw) & (~tw)
            v2 = (~test3) & (~testo) & (~to)

            upwind_w[v1, 1] = True
            upwind_o[v2, 1] = True

            tt = test3
            if tt.sum() > 0:

                upwind_w[tt] = False
                upwind_o[tt] = False

                centroids = centroid_volumes[v0[tt]]
                delta_z = centroids[:,1] - centroids[:,0]
                delta_z = delta_z[:,2]

                upgw = np.full((tt.sum(), 2), False, dtype=bool)
                upgo = np.full((tt.sum(), 2), False, dtype=bool)

                tz = delta_z >= 0

                upgw[tz, 1] = True
                upgo[tz, 0] = True

                tz = ~tz
                upgw[tz, 0] = True
                upgo[tz, 1] = True

                upwind_w[tt] = upgw
                upwind_o[tt] = upgo

            ##############################################
            ## antigo
            # upwind_w[tw, 0] = True
            # upwind_o[to, 0] = True
            #
            # tw = ~tw
            # to = ~to
            # upwind_w[tw, 1] = True
            # upwind_o[to, 1] = True
            #
            # tt = np.absolute(flux_total) < k0
            # upwind_w[tt] = False
            # upwind_o[tt] = False
            #
            # centroids = centroid_volumes[v0[tt]]
            # delta_z = centroids[:,1] - centroids[:,0]
            # delta_z = delta_z[:,2]
            #
            # upgw = np.full((tt.sum(), 2), False, dtype=bool)
            # upgo = np.full((tt.sum(), 2), False, dtype=bool)
            #
            # tz = delta_z >= 0
            #
            # upgw[tz, 1] = True
            # upgo[tz, 0] = True
            #
            # tz = ~tz
            # upgw[tz, 0] = True
            # upgo[tz, 1] = True
            #
            # upwind_w[tt] = upgw
            # upwind_o[tt] = upgo
            #########################################



        else:
            tw = flux_w_faces >= 0
            to = flux_o_faces >= 0
            upwind_w[tw, 0] = True
            upwind_o[to, 0] = True
            tw = ~tw
            to = ~to
            upwind_w[tw, 1] = True
            upwind_o[to, 1] = True

        # verif1 = upwind_w[:,0] ^ upwind_w[:,1]
        # verif2 = upwind_o[:,0] ^ upwind_o[:,1]
        # verif1 = ~verif1
        # verif2 = ~verif2

        # pdb.set_trace()

        self.test_upwind_dup(upwind_w, upwind_o)
        return upwind_w, upwind_o

    def visualize_upwind_vec(self, internal_faces, abs_u_normal_internal_faces, centroid_volumes, volumes_adj_internal_faces, upwind_w, upwind_o):

        u_normal_internal_faces = abs_u_normal_internal_faces
        v0 = volumes_adj_internal_faces
        c0 = centroid_volumes[v0[:,0]]
        c1 = centroid_volumes[v0[:,1]]
        dc = c1-c0
        norm = np.linalg.norm(dc, axis=1).reshape(dc.shape[0], 1)
        dcu = dc/norm
        upwind_w_out = np.zeros((len(internal_faces), 3), dtype=int)
        upwind_w_out[upwind_w[:,0]] = dcu[upwind_w[:,0]].astype(int)
        upwind_w_out[upwind_w[:,1]] = -dcu[upwind_w[:,1]].astype(int)

        upwind_o_out = np.zeros((len(internal_faces), 3), dtype=int)
        upwind_o_out[upwind_o[:,0]] = dcu[upwind_o[:,0]].astype(int)
        upwind_o_out[upwind_o[:,1]] = -dcu[upwind_o[:,1]].astype(int)

        return upwind_w_out, upwind_o_out

    def get_g_source_w_o_internal_faces_dep(self, nkga_internal_faces, mobility_w_internal_faces, mobility_o_internal_faces, rho_w, rho_o, hi):

        assert len(nkga_internal_faces) ==  len(mobility_w_internal_faces) == len(mobility_o_internal_faces)

        g_source_w_internal_faces = nkga_internal_faces*mobility_w_internal_faces*rho_w*hi.sum(axis=1)
        g_source_o_internal_faces = nkga_internal_faces*mobility_o_internal_faces*rho_o*hi.sum(axis=1)

        return g_source_w_internal_faces, g_source_o_internal_faces

    def get_g_source_w_o_internal_faces(self, areas_internal_faces, u_direction_internal_faces, g_velocity_w_internal_faces, g_velocity_o_internal_faces):

        gw = (g_velocity_w_internal_faces*u_direction_internal_faces).sum(axis=1)
        gw = gw*areas_internal_faces
        go = (g_velocity_o_internal_faces*u_direction_internal_faces).sum(axis=1)
        go = go*areas_internal_faces

        return gw, go

    def get_g_velocity_w_o_internal_faces(self, gravity_vector, mobility_w_internal_faces, mobility_o_internal_faces, rho_w, rho_o, hi, keq_internal_faces):

        ni = len(keq_internal_faces)
        mob_w_int_f = mobility_w_internal_faces.reshape(ni, 1)
        mob_o_int_f = mobility_o_internal_faces.reshape(ni, 1)
        keq = keq_internal_faces.reshape(ni, 1)
        hi2 = hi.sum(axis=1).reshape(ni, 1)

        v_gw = rho_w*mob_w_int_f*gravity_vector*keq*hi2
        v_go = rho_o*mob_o_int_f*gravity_vector*keq*hi2

        return v_gw, v_go


