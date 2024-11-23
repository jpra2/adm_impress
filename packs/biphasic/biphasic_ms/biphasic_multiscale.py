import numpy as np
import scipy.sparse as sp
from ..biphasic_tpfa import BiphasicTpfa
from .biphasic_multiscale_properties import bipahsicMultiscaleProperties
import time

class BiphasicTpfaMultiscale(BiphasicTpfa, bipahsicMultiscaleProperties):

    def __init__(self, M, data_impress, elements_lv0, wells, data_name: str='BiphasicTpfaMultiscale.npz'):
        super().__init__(M, data_impress, elements_lv0, wells, data_name=data_name)
        
    def update_flux_w_and_o_volumes(self):
        internal_faces = self.elements_lv0['internal_faces']
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
