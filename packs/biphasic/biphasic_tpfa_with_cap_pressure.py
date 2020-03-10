from .biphasic_tpfa import BiphasicTpfa
from ..utils.capillary_pressure import capillaryPressureBiphasic
import numpy as np
import scipy.sparse as sp
import time

class biphasicTpfaCapPressure(BiphasicTpfa):

    def __init__(self, M, data_impress, elements_lv0, wells, data_name: str='BiphasicTpfaCapPressure.npz'):
        super().__init__(M, data_impress, elements_lv0, wells, data_name=data_name)
        self.capillary_model = capillaryPressureBiphasic()
        # self.remap_gids()

        if not self.load:
            self.data_impress['capillary_pressure'] = self.capillary_model.get_pcow_from_sat(self.data_impress['saturation'])

    def get_cap_flux(self):

        # import pdb; pdb.set_trace()

        internal_faces = self.elements_lv0['internal_faces']
        v0 = self.elements_lv0['neig_internal_faces']
        areas = self.data_impress['area'][internal_faces]
        k_harm = self.data_impress['k_harm'][internal_faces]
        dh = self.data_impress['dist_cent'][internal_faces]
        # lambda_w = self.data_impress['lambda_w'][v0[self._data['upwind_identificate']]]
        lambda_w = self.lambda_w_internal_faces
        cap_pressure = self.data_impress['capillary_pressure'][v0]

        flux_cap_faces, flux_cap_volumes = self.capillary_model.get_capillary_pressure_flux(
            internal_faces, v0, areas, k_harm, dh, lambda_w, cap_pressure
            )

        self.data_impress['flux_cap_faces'] = np.zeros(len(self.data_impress['flux_cap_faces']))
        self.data_impress['flux_cap_faces'][internal_faces] = flux_cap_faces
        self.data_impress['flux_cap_volumes'] = flux_cap_volumes

    def get_b_cap(self):
        b2 = self.data_impress['flux_cap_volumes'].copy()
        b2[self.wells['ws_p']] = np.zeros(len(self.wells['ws_p']))
        return b2

    def update_flux_w_and_o_volumes(self) -> None:

        ## testando sem gravidade

        # import pdb; pdb.set_trace()

        vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        v0 = vols_viz_internal_faces
        internal_faces = self.elements_lv0['internal_faces']
        # total_flux_faces = self.data_impress['flux_faces']
        # fw_faces = self.data_impress['fw_faces']
        # fw_vol = self.data_impress['fw_vol']
        # flux_volumes = self.data_impress['flux_volumes']
        ws_prod = self.wells['ws_prod']
        ws_inj = self.wells['ws_inj']
        grav_source_term_water_volumes = self._data['grav_source_term_water_volumes']
        lambda_w_internal_faces = self.lambda_w_internal_faces
        lambda_o_internal_faces = self.lambda_o_internal_faces
        # fw_internal_faces = self.data_impress['fw_vol'][v0[self._data['upwind_identificate']]]
        # fw_internal_faces = self.fw_internal_faces

        areas_internal_faces = self.data_impress['area'][internal_faces]
        k_harm_internal_faces = self.data_impress['k_harm'][internal_faces]
        dh_internal_faces = self.data_impress['dist_cent'][internal_faces]

        x = self.data_impress['water_pressure']
        ps1 = x[v0[:, 1]]
        ps0 = x[v0[:, 0]]

        # k0 = 1
        k0 = -1
        # k1 = 1
        k1 = -1
        # k2 = 1
        k2 = -1

        # import pdb; pdb.set_trace()

        # flux_w_faces = fw_faces*total_flux_faces
        # flux_w_internal_faces2 = flux_w_faces[internal_faces]
        flux_w_internal_faces = k0*((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_w_internal_faces/dh_internal_faces + k1*self._data['grav_source_term_water_faces'][internal_faces])

        # x = self.data_impress['pressure']
        # ps1 = x[v0[:, 1]]
        # ps0 = x[v0[:, 0]]
        # flux_o_internal_faces = k0*((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_o_internal_faces/dh_internal_faces + k1*(self.data_impress['flux_grav_faces'][internal_faces] - self._data['grav_source_term_water_faces'][internal_faces]))
        # flux_w_internal_faces = (total_flux_faces[internal_faces] + k_harm_internal_faces*areas_internal_faces*self.lambda_o_internal_faces*(ps1 - ps0)/dh_internal_fa
        # flux_w_internal_faces = fw_internal_faces*(total_flux_faces[internal_faces] - k_harm_internal_faces*areas_internal_faces*self.lambda_o_internal_faces*self.capillary_model.get_gradient_cap_presure(dh_internal_faces, cap_pressure))
        # flux_w_internal_faces = (total_flux_faces[internal_faces] + k_harm_internal_faces*areas_internal_faces*self.lambda_o_internal_faces*(ps1 - ps0)/dh_internal_faces)

        x = self.data_impress['pressure']
        ps1 = x[v0[:, 1]]
        ps0 = x[v0[:, 0]]

        # total_flux_internal_faces = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*self.lambda_t_internal_faces/dh_internal_faces - self.data_impress['flux_grav_faces'][internal_faces] - self._data_impress['flux_cap_faces'][internal_faces])
        total_flux_internal_faces = k0*((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*self.lambda_t_internal_faces/dh_internal_faces + k1*self.data_impress['flux_grav_faces'][internal_faces] + k2*self.data_impress['flux_cap_faces'][internal_faces])
        # lambda_o_internal_faces = self.data_impress['lambda_o'][v0[self._data['upwind_identificate']]]
        # flux_o_internal_faces = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_o_internal_faces/dh_internal_faces - (self.data_impress['flux_grav_faces'][internal_faces] - self._data['grav_source_term_water_internal_faces']))

        flux_o_internal_faces = total_flux_internal_faces - flux_w_internal_faces

        # soma = flux_w_internal_faces + flux_o_internal_faces
        # vv = abs(soma - total_flux_faces[internal_faces])


        # verif = np.allclose(flux_w_internal_faces, flux_w_internal_faces2)
        # import pdb; pdb.set_trace()

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_w_internal_faces, -flux_w_internal_faces]).flatten()
        flux_w_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        data = np.array([total_flux_internal_faces, -total_flux_internal_faces]).flatten()
        total_flux_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
        flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        # self.print_test()
        # import pdb; pdb.set_trace()
        # import pdb; pdb.set_trace()

        # flux_w_volumes[ws_prod] -= flux_volumes[ws_prod]*fw_vol[ws_prod]
        # flux_w_volumes[ws_inj] -= flux_volumes[ws_inj]*fw_vol[ws_inj]
        all_wells = self.wells['all_wells']
        flux_w_volumes[all_wells] -= total_flux_volumes[all_wells]*self.fw_volumes[all_wells]

        # flux_o_internal_faces = total_flux_faces[internal_faces] - flux_w_internal_faces

        # lambda_o_internal_faces = self.data_impress['lambda_t'][v0[self._data['upwind_identificate']]] - self.data_impress['lambda_w'][v0[self._data['upwind_identificate']]]
        # flux_o_internal_faces2 = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_o_internal_faces/dh_internal_faces - (self.data_impress['flux_grav_faces'][internal_faces] - self._data['grav_source_term_water_internal_faces']))

        # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        # cols = np.repeat(0, len(lines))
        # data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
        # flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        flux_o_volumes[all_wells] -= total_flux_volumes[all_wells]*(1-self.fw_volumes[all_wells])
        #
        # flux_w_faces = np.zeros(len(self.data_impress['flux_w_faces']))
        # flux_w_faces[internal_faces] = flux_w_internal_faces

        import pdb; pdb.set_trace()


        self.data_impress['flux_w_faces'] = np.zeros(len(self.data_impress['flux_w_faces']))
        self.data_impress['flux_w_faces'][internal_faces] = flux_w_internal_faces
        self.data_impress['flux_o_volumes'] = flux_o_volumes
        self.data_impress['flux_w_volumes'] = flux_w_volumes
        self.data_impress['flux_volumes'] = total_flux_volumes

    def remap_gids(self):
        import pdb; pdb.set_trace()
        self._data['remaped_gids'] = -1*np.ones(len(self.elements_lv0['volumes']), dtype=int)
        active_vols = self.active_vols()
        self._data['remaped_gids'][active_vols] = np.arange(len(active_vols))

    def run(self):

        T, b = self.get_T_and_b()
        self.get_cap_flux()
        b += self.get_b_cap()
        # import pdb; pdb.set_trace()
        self.data_impress['pressure'] = self.solver.direct_solver(T, b)
        self.data_impress['water_pressure'] = self.data_impress['pressure'] - self.data_impress['capillary_pressure']
        # self.update_flux_w_and_o_volumes()
        # import pdb; pdb.set_trace()
        # flux_test = -(self['Tini']*self.data_impress['pressure'] - self.data_impress['flux_grav_volumes'] - self.data_impress['flux_cap_volumes'])
        # self.data_impress['flux_volumes_test'] = flux_test
        # self.data_impress['verif_po'] = np.absolute(self.data_impress['flux_volumes'] - self.data_impress['flux_volumes_test'])
        # self.print_test()
        # import pdb; pdb.set_trace()

        # self.get_flux_faces_and_volumes()

        # self.data_impress['water_pressure'] = self.data_impress['pressure'] - self.data_impress['capillary_pressure']
        # flux_total = T*self.data_impress['pressure']
        # flux_total -= self.data_impress['flux_cap_volumes']
        # flux_total -= self.data_impress['flux_grav_volumes']
        # self.data_impress['flux_volumes'] = flux_total
        # ft = flux_total
        # rr0 = ft[self.wells['ws_p']]
        # rr1 = ft[self.wells['ws_q']]
        # hh = self.data_impress['flux_cap_volumes'][self.wells['ws_p']]
        # import pdb; pdb.set_trace()

    def run_2(self, save=False):

        t0 = time.time()
        self.update_flux_w_and_o_volumes()
        self.update_delta_t()
        self.update_saturation()
        self.update_t()
        self.update_vpi()
        self.update_relative_permeability()
        self.update_mobilities()
        self.update_transmissibility()
        self.update_loop()
        t1 = time.time()
        dt = t1-t0
        self.update_current_biphasic_results(dt)

        if save:
            self.save_infos()
