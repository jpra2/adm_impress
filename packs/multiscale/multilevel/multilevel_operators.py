from ...data_class.data_manager import DataManager
from ..operators.prolongation.AMS.ams_tpfa import AMSTpfa
from ..operators.prolongation.AMS.ams_mpfa import AMSMpfa
import numpy as np
import scipy.sparse as sp

def get_gids_primalids_dualids(gids, primal_ids, dual_ids):

    gids2 = np.unique(gids)
    ids = np.arange(len(gids))
    primal_ids2 = []
    dual_ids2 = []
    for i in gids2:
        test = ids[gids==i]
        primal_id = np.unique(primal_ids[test])
        dual_id = np.unique(dual_ids[test])
        if len(primal_id) > 1 or len(dual_id) > 1:
            raise ValueError('erro get_gids_primalids_dualids')
        primal_ids2.append(primal_id[0])
        dual_ids2.append(dual_id[0])

    primal_ids2 = np.array(primal_ids2)
    dual_id = np.array(dual_ids2)

    return gids2, primal_ids2, dual_ids2

class MultilevelOperators(DataManager):
    def __init__(self,
                n_levels,
                data_impress,
                data_name='MultilevelOperators.npz',
                load=False):

        super().__init__(data_name=data_name, load=load)

        self.n_levels = n_levels
        self.data_impress = data_impress

        self.restriction = 'restriction_level_'
        self.prolongation = 'prolongation_level_'
        self.infos_level = 'infos_level_'
        self.gid_n = 'gid'
        self.primal_id_n = 'primal_id'
        self.dual_id_n = 'dual_id'
        self.operators = dict()

        if load == False:
            self.get_initial_infos()

        self.get_operators(load=load)

        if load == False:
            self.export_to_npz()

    def get_initial_infos(self):
        t1 = np.dtype(int)
        dt = [(self.gid_n, t1),
              (self.primal_id_n, t1),
              (self.dual_id_n, t1)]

        for n in range(self.n_levels):

            gids = self.data_impress['GID_' + str(n)]
            primal_ids = self.data_impress['GID_' + str(n+1)]
            dual_ids = self.data_impress['DUAL_' + str(n+1)]
            if n > 0:
                gids, primal_ids, dual_ids = get_gids_primalids_dualids(gids, primal_ids, dual_ids)
            tam = len(gids)
            tam2 = len(np.unique(primal_ids))
            structured_array = np.zeros(tam, dtype=dt)
            structured_array[self.gid_n] = gids
            structured_array[self.primal_id_n] = primal_ids
            structured_array[self.dual_id_n] = dual_ids
            self._data[self.infos_level + str(n+1)] = structured_array
            data_r = np.ones(tam)
            lines = primal_ids
            cols = gids
            OR = sp.csc_matrix((data_r,(lines,cols)), shape=(tam2, tam))
            self._data[self.restriction + str(n+1)] = OR

    def get_operators(self, load=False):

        for n in range(self.n_levels):
            level = n+1
            infos = self._data[self.infos_level + str(level)]
            gid = infos[self.gid_n]
            primal_id = infos[self.primal_id_n]
            dual_id = infos[self.dual_id_n]
            interns = gid[dual_id==0]
            faces = gid[dual_id==1]
            edges = gid[dual_id==2]
            vertices = gid[dual_id==3]
            if level == 1:
                operator = AMSTpfa
            else:
                operator = AMSMpfa
            self.operators[str(level)] = operator(interns,
                                                 faces,
                                                 edges,
                                                 vertices,
                                                 gid,
                                                 primal_id,
                                                 load=load)

    def run(self, T: 'fine transmissibility'):
        T_ant = T.copy()
        for n in range(self.n_levels):
            level = n+1
            OP = self.operators[str(level)].run(T_ant)
            self._data[self.prolongation + str(level)] = OP
            OR = self._data[self.restriction + str(level)]
            T_ant = OR*T_ant*OP
