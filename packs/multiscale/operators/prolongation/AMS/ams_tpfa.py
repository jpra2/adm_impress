from .....flux_schemes.tpfa_scheme import TpfaScheme
import scipy.sparse as sp
import numpy as np


class AMSTpfa(TpfaScheme):
    name = 'AMSTpfa_'
    id = 1

    def __init__(self,
        T: 'transmissibility_matrix',
        gids: 'global_ids',
        reordered_ids: 'ids_for_wirebasket_pattern',
        ni: 'number_of_interns',
        nf: 'number of faces',
        ne: 'number of edges',
        nv: 'number_of_vertices',
        name: 'unique name for data'='',
        load=False):

        if name == '':
            data_name = AMSTpfa.name + str(AMSTpfa.id) + '.npz'
        else:
            data_name = AMSTpfa.name + name + '.npz'
        super().__init__(data_name=data_name, load=load)
        AMSTpfa.id += 1
        self['gids'] = gids.copy()
        self['reordered_ids'] = reordered_ids.copy()
        self.nv = nv
        self.ne = ne
        self.nf = nf
        self.ni = ni

    def get_permutation_matrix(self):
        lines = self['reordered_ids']
        cols = self['gids']
        data = np.ones(len(lines))
        self['permutation'] = sp.csc_matrix((data, (lines, cols)), shape=(self.nv, self.nv))

    def get_prolongation_operator(self):
        pass

    def get_AS(self):

        Tmod = self['Tini'].copy().tolil()
        ni = self.ni
        nf = self.nf
        ne = self.ne
        nv = self.nv

        nni = ni
        nnf = nf + nni
        nne = ne + nnf
        nnv = nv + nne

        #internos
        Aii = Tmod[0:nni, 0:nni]
        Aif = Tmod[0:nni, nni:nnf]

        #faces
        Aff = Tmod[nni:nnf, nni:nnf]
        Afe = Tmod[nni:nnf, nnf:nne]
        soma = Aif.transpose().sum(axis=1)
        d1 = np.matrix(Aff.diagonal()).reshape([nf, 1])
        d1 += soma
        Aff.setdiag(d1)

        #arestas
        Aee = Tmod[nnf:nne, nnf:nne]
        Aev = Tmod[nnf:nne, nne:nnv]
        soma = Afe.transpose().sum(axis=1)
        d1 = np.matrix(Aee.diagonal()).reshape([ne, 1])
        d1 += soma
        Aee.setdiag(d1)
        Ivv = sp.identity(nv)

        As = {}
        As['Aii'] = Aii
        As['Aif'] = Aif
        As['Aff'] = Aff
        As['Afe'] = Afe
        As['Aee'] = Aee
        As['Aev'] = Aev
        As['Ivv'] = Ivv

        return As

    self.get_transmissibility_matrix_without_boundary_conditions()
