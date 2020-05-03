import numpy as np
import scipy.sparse as sp
from ..biphasic_tpfa import BiphasicTpfa
from .biphasic_multiscale_properties import bipahsicMultiscaleProperties

class BiphasicTpfaMultiscale(BiphasicTpfa, bipahsicMultiscaleProperties):

    def __init__(self, M, data_impress, elements_lv0, wells, data_name: str='BiphasicTpfa.npz'):
        super().__init__(M, data_impress, elements_lv0, wells, data_name=data_name)


    def update_flux_w_and_o_volumes(self):

        import pdb; pdb.set_trace()

        # flux_faces = self.data_impress['flux_faces']
        # fw_faces = self.data_impress['fw_faces']
        # internal_faces = self.elements_lv0['internal_faces']
        # flux_w_internal_faces = fw_faces[internal_faces]*flux_faces[internal_faces]
        # flux_o_internal_faces = (1 - fw_faces[internal_faces])*flux_faces[internal_faces]
        # ws_prod = self.wells['ws_prod']
        # ws_inj = self.wells['ws_inj']
        # fw_vol = self.data_impress['fw_vol']
        # v0 = self.elements_lv0['neig_internal_faces']
        # flux_volumes = self.data_impress['flux_volumes']

        # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        # cols = np.repeat(0, len(lines))
        # data = np.array([flux_w_internal_faces, -flux_w_internal_faces]).flatten()
        # flux_w_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        # flux_w_volumes[ws_prod] -= flux_volumes[ws_prod]*fw_vol[ws_prod]
        # flux_w_volumes[ws_inj] -= flux_volumes[ws_inj]*fw_vol[ws_inj]
        #
        # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        # cols = np.repeat(0, len(lines))
        # data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
        # flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        # flux_o_volumes[ws_prod] -= flux_volumes[ws_prod]*(1 - fw_vol[ws_prod])

        self.data_impress['flux_w_faces'] = np.zeros(len(self.data_impress['flux_w_faces']))
        self.data_impress['flux_w_faces'][internal_faces] = self.flux_w_internal_faces_ms
        self.data_impress['flux_w_volumes'] = self.flux_w_volumes_ms
        self.data_impress['flux_o_volumes'] = self.flux_o_volumes_ms

        flux_o_internal_faces = self.flux_o_internal_faces_ms

        u_normal = self.data_impress['u_normal']
        flux_w_vec_internal_faces = u_normal[internal_faces]*self.data_impress['flux_w_faces'][internal_faces].reshape([len(internal_faces), 1])
        flux_o_vec_internal_faces = u_normal[internal_faces]*flux_o_internal_faces.reshape([len(internal_faces), 1])
        self.data_impress['flux_w_faces_vec'] = np.zeros(self.data_impress['flux_w_faces_vec'].shape)
        self.data_impress['flux_o_faces_vec'] = np.zeros(self.data_impress['flux_w_faces_vec'].shape)
        self.data_impress['flux_w_faces_vec'][internal_faces] = flux_w_vec_internal_faces
        self.data_impress['flux_o_faces_vec'][internal_faces] = flux_o_vec_internal_faces
