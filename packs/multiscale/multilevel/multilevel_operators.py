from ...data_class.data_manager import DataManager
from ..operators.prolongation.AMS.ams_tpfa import AMSTpfa
from ..operators.prolongation.AMS.ams_mpfa import AMSMpfa
import numpy as np
import scipy.sparse as sp
import os

def get_gids_primalids_dualids(gids, primal_ids, dual_ids):

    gids2 = np.unique(gids)
    # ids = np.arange(len(gids))
    primal_ids2 = []
    dual_ids2 = []
    for i in gids2:
        # test = ids[gids==i]
        test = gids==i
        primal_id = np.unique(primal_ids[test])
        dual_id = np.unique(dual_ids[test])
        if len(primal_id) > 1 or len(dual_id) > 1:
            raise ValueError('erro get_gids_primalids_dualids')
        primal_ids2.append(primal_id[0])
        dual_ids2.append(dual_id[0])

    primal_ids2 = np.array(primal_ids2)
    dual_id = np.array(dual_ids2)

    return gids2, primal_ids2, dual_ids2

def manter_vizinhos_de_face(T, ids, neigh_ids):
    T2 = T.copy()
    nl = T2.shape[0]
    map_neig = dict(zip(ids, neigh_ids))

    for i in range(nl):
        neigs = set(map_neig[i])
        data = sp.find(T[i])
        cols = data[1]
        vals = data[2]
        # id_line = cols==i
        viz_rem = []

        for viz, val in zip(cols, vals):
            if set([viz]) & neigs or viz == i:
                pass
            else:
                viz_rem.append(viz)
                T2[i, i] += val
        viz_rem = np.array(viz_rem)
        T2[np.repeat(i, len(viz_rem)), viz_rem] = np.zeros(len(viz_rem))

    data_t2 = sp.find(T2)
    T2 = sp.csc_matrix((data_t2[2], (data_t2[0], data_t2[1])), shape=(nl, nl))

    return T2

class MultilevelOperators(DataManager):
    def __init__(self,
                n_levels,
                data_impress,
                ml_data,
                data_name='MultilevelOperators.npz',
                load=False):

        super().__init__(data_name=data_name, load=load)
        self.load = load
        self.ml_data = ml_data

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

        if load == True:
            for n in range(n_levels):
                prol_name = self.prolongation + str(n+1)
                rest_name = self.restriction + str(n+1)
                OP = sp.load_npz(os.path.join('flying', prol_name + '.npz'))
                OR = sp.load_npz(os.path.join('flying', rest_name + '.npz'))
                self._data[prol_name] = OP
                self._data[rest_name] = OR

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

    def run(self, T: 'fine transmissibility without boundary conditions'):

        T_ant = T.copy()
        for n in range(self.n_levels):
            level = n+1
            OP = self.operators[str(level)].run(T_ant)
            self._data[self.prolongation + str(level)] = OP
            OR = self._data[self.restriction + str(level)]

            sp.save_npz(os.path.join('flying', self.prolongation + str(level) + '.npz'), OP)
            sp.save_npz(os.path.join('flying', self.restriction + str(level) + '.npz'), OR)

            if level == self.n_levels:
                continue
            T_ant = OR*T_ant*OP
            cids_neigh = self.ml_data['coarse_id_neig_face_level_'+str(level)]
            cids_level = self.ml_data['coarse_primal_id_level_'+str(level)]
            # T_ant = manter_vizinhos_de_face(T_ant, cids_level, cids_neigh)
