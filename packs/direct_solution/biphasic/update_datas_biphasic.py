import pdb
from ... import directories as direc
from . import relative_permeability
import numpy as np



class updateDatas:

    def __init__(self, M):
        name_relative_permeability = direc.data_loaded['biphasic_data']['relative_permeability']
        self.relative_permeability = getattr(relative_permeability, name_relative_permeability)
        self.mesh = M
        self.mi_w = direc.data_loaded['biphasic_data']['mi_w']
        self.mi_o = direc.data_loaded['biphasic_data']['mi_o']

    def update_mobility(self, krws, kros):
        M = self.mesh
        n = len(volumes)
        lambda_w = krws/self.mi_w
        lambda_o = kros/self.mi_o
        lambda_t = lamb_w + lamb_o
        fw_vol = lamb_w/lamb_t

        M.data.variables[M.data.variables_impress['lambda_w']] = lambda_w
        M.data.variables[M.data.variables_impress['lambda_o']] = lambda_o
        M.data.variables[M.data.variables_impress['lambda_t']] = lambda_t
        M.data.variables[M.data.variables_impress['fw_vol']] = fw_vol

    def update_upwind(self):
        M = self.mesh
        n_internal_faces = len(flux_faces)
        ids = np.arange(n_faces)

        internal_faces = M.data.elements_lv0[direc.entities_lv0_0[0]]
        transmissibility_faces = M.data.variables[M.data.variables_impress['transmissibility']]
        t0 = transmissibility_internal_faces
        pretransmissibility_faces = M.data.variables[M.data.variables_impress['pretransmissibility']]

        fw_faces = np.zeros(n_internal_faces)
        lambda_t_faces = fw_faces.copy()
        transmissibility_internal_faces = fw_faces.copy()

        positive_ids = ids[flux_faces >= 0]
        negative_ids = np.setdiff1d(ids, positive_ids)

        fw_faces[positive_ids] = fw_vols[adjs_internal_faces[positive_ids, np.repeat(0, len(positive_ids), dtype=np.int64)]]
        fw_faces[negative_ids] = fw_vols[adjs_internal_faces[negative_ids, np.repeat(1, len(negative_ids), dtype=np.int64)]]

        lambda_t_faces[positive_ids] = total_mobilities_vols[adjs_internal_faces[positive_ids, np.repeat(0, len(positive_ids), dtype=np.int64)]]
        lambda_t_faces[negative_ids] = total_mobilities_vols[adjs_internal_faces[negative_ids, np.repeat(1, len(negative_ids), dtype=np.int64)]]

        transmissibility_internal_faces = pretransmissibility_faces[internal_faces]*lambda_t_faces

        M.data.variables[M.data.variables_impress['lambda_t_faces']][internal_faces] = lambda_t_faces
        M.data.variables[M.data.variables_impress['fw_face']][internal_faces] = fw_faces
        M.data.variables[M.data.variables_impress['transmissibility']][internal_faces] = transmissibility_internal_faces

    def update_saturation_upwind(self)
