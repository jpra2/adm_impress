from .biphasic_tpfa import BiphasicTpfa
from ..utils.capillary_pressure import capillaryPressureBiphasic
import numpy as np
import scipy.sparse as sp

class biphasicTpfaCapPressure(BiphasicTpfa):

    def __init__(self, M, data_impress, elements_lv0, wells, data_name: str='BiphasicTpfaCapPressure.npz'):
        super().__init__(M, data_impress, elements_lv0, wells, data_name=data_name)
        self.capillary_model = capillaryPressureBiphasic()

        if not self.load:
            self.data_impress['capillary_pressure'] = self.capillary_model.get_pcow_from_sat(self.data_impress['saturation'])

    def get_cap_flux(self):

        internal_faces = self.elements_lv0['internal_faces']
        v0 = self.elements_lv0['neig_internal_faces']
        areas = self.data_impress['area'][internal_faces]
        k_harm = self.data_impress['k_harm'][internal_faces]
        dh = self.data_impress['dist_cent']
        lambda_w = self.data_impress['lambda_w'][v0[self._data['upwind_identificate']]]
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

        vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        v0 = vols_viz_internal_faces
        internal_faces = self.elements_lv0['internal_faces']
        total_flux_faces = self.data_impress['flux_faces']
        fw_faces = self.data_impress['fw_faces']
        fw_vol = self.data_impress['fw_vol']
        flux_volumes = self.data_impress['flux_volumes']
        ws_prod = self.wells['ws_prod']
        ws_inj = self.wells['ws_inj']
        grav_source_term_water_volumes = self._data['grav_source_term_water_volumes']
        lambda_w_internal_faces = self.data_impress['lambda_w'][v0[self._data['upwind_identificate']]]

        areas_internal_faces = self.data_impress['area'][internal_faces]
        k_harm_internal_faces = self.data_impress['k_harm'][internal_faces]
        dh_internal_faces = self.data_impress['dist_cent'][internal_faces]

        x = self.data_impress['water_pressure']
        ps0 = x[v0[:, 0]]
        ps1 = x[v0[:, 1]]

        # flux_w_faces = fw_faces*total_flux_faces
        # flux_w_internal_faces2 = flux_w_faces[internal_faces]
        flux_w_internal_faces = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_w_internal_faces/dh_internal_faces - self._data['grav_source_term_water_internal_faces'])

        # lambda_o_internal_faces = self.data_impress['lambda_o'][v0[self._data['upwind_identificate']]]
        # flux_o_internal_faces = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_o_internal_faces/dh_internal_faces - (self.data_impress['flux_grav_faces'][internal_faces] - self._data['grav_source_term_water_internal_faces']))

        # flux_o_internal_faces = total_flux_faces[internal_faces] - flux_w_internal_faces

        # soma = flux_w_internal_faces + flux_o_internal_faces
        # vv = abs(soma - total_flux_faces[internal_faces])


        # verif = np.allclose(flux_w_internal_faces, flux_w_internal_faces2)
        # import pdb; pdb.set_trace()

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_w_internal_faces, -flux_w_internal_faces]).flatten()
        flux_w_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        self.data_impress['flux_w_volumes'] = flux_w_volumes
        self.print_test()
        import pdb; pdb.set_trace()
        # import pdb; pdb.set_trace()

        flux_w_volumes[ws_prod] -= flux_volumes[ws_prod]*fw_vol[ws_prod]
        flux_w_volumes[ws_inj] -= flux_volumes[ws_inj]*fw_vol[ws_inj]

        # flux_o_internal_faces = total_flux_faces[internal_faces] - flux_w_internal_faces

        # lambda_o_internal_faces = self.data_impress['lambda_t'][v0[self._data['upwind_identificate']]] - self.data_impress['lambda_w'][v0[self._data['upwind_identificate']]]
        # flux_o_internal_faces2 = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_o_internal_faces/dh_internal_faces - (self.data_impress['flux_grav_faces'][internal_faces] - self._data['grav_source_term_water_internal_faces']))

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
        flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        flux_o_volumes[ws_prod] -= flux_volumes[ws_prod]*(1 - fw_vol[ws_prod])

        flux_w_faces = np.zeros(len(self.data_impress['flux_w_faces']))
        flux_w_faces[internal_faces] = flux_w_internal_faces

        self.data_impress['flux_w_faces'] = flux_w_faces
        self.data_impress['flux_w_volumes'] = flux_w_volumes
        self.data_impress['flux_o_volumes'] = flux_o_volumes

    def run(self):

        T, b = self.get_T_and_b()
        self.get_cap_flux()
        b += self.get_b_cap()
        p = self.solver.direct_solver(T, b)
        self.data_impress['pressure'] = p
        self.data_impress['water_pressure'] = p - self.data_impress['capillary_pressure']
        self.update_flux_w_and_o_volumes()
