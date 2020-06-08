import numpy as np
import scipy.sparse as sp

class bipahsicMultiscaleProperties:

    @property
    def flux_internal_faces_ms(self):
        internal_faces = self.elements_lv0['internal_faces']
        return self.data_impress['flux_faces'][internal_faces]

    @property
    def flux_w_internal_faces_ms(self):
        flux_internal_faces = self.flux_internal_faces_ms
        fw_faces = self.data_impress['fw_faces']
        internal_faces = self.elements_lv0['internal_faces']
        flux_w_internal_faces = fw_faces[internal_faces]*flux_internal_faces
        return flux_w_internal_faces

    @property
    def flux_o_internal_faces_ms(self):
        flux_internal_faces = self.flux_w_internal_faces_ms
        fw_faces = self.data_impress['fw_faces']
        internal_faces = self.elements_lv0['internal_faces']
        flux_o_internal_faces = (1 - fw_faces[internal_faces])*flux_internal_faces
        return flux_o_internal_faces

    @property
    def flux_w_volumes_ms(self):
        flux_w_internal_faces = self.flux_w_internal_faces_ms
        fw_faces = self.data_impress['fw_faces']
        ws_prod = self.wells['ws_prod']
        ws_inj = self.wells['ws_inj']
        fw_vol = self.data_impress['fw_vol']
        v0 = self.elements_lv0['neig_internal_faces']
        flux_volumes = self.data_impress['flux_volumes']

        # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        lines=np.concatenate(v0.T)
        # cols = np.repeat(0, len(lines))
        # data = np.array([flux_w_internal_faces, -flux_w_internal_faces]).flatten()
        data=np.concatenate([flux_w_internal_faces,-flux_w_internal_faces])

        # flux_w_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        flux_w_volumes=np.bincount(lines,weights=data)
        # import pdb; pdb.set_trace()
        flux_w_volumes[ws_prod] -= flux_volumes[ws_prod]*fw_vol[ws_prod]
        flux_w_volumes[ws_inj] -= flux_volumes[ws_inj]*fw_vol[ws_inj]

        # self.data_impress['verif_po'][:] = 0.0
        # self.data_impress['verif_po'][ws_prod] = flux_volumes[ws_prod]*fw_vol[ws_prod]
        # self.data_impress['verif_po'][ws_inj] = flux_volumes[ws_inj]*fw_vol[ws_inj]

        return flux_w_volumes

    @property
    def flux_o_volumes_ms(self):
        flux_o_internal_faces = self.flux_o_internal_faces_ms
        fw_faces = self.data_impress['fw_faces']
        ws_prod = self.wells['ws_prod']
        ws_inj = self.wells['ws_inj']
        fw_vol = self.data_impress['fw_vol']
        v0 = self.elements_lv0['neig_internal_faces']
        flux_volumes = self.data_impress['flux_volumes']

        # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        lines=np.concatenate(v0.T)
        # cols = np.repeat(0, len(lines))
        # data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
        data = np.concatenate([flux_o_internal_faces, -flux_o_internal_faces])
        # flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        flux_o_volumes = np.bincount(lines,weights=data)
        flux_o_volumes[ws_prod] -= flux_volumes[ws_prod]*(1 - fw_vol[ws_prod])

        return flux_o_volumes
