from .....data_class.data_manager import DataManager
import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg
import time

class AMSTpfa:
    # name = 'AMSTpfa_'
    id = 1

    def __init__(self,
        internals,
        faces,
        edges,
        vertices,
        gids: 'global_ids',
        primal_ids: 'primal_ids',
        data_name='AMSTpfa_',
        load=False):

        # data_name = AMSTpfa.name + str(AMSTpfa.id) + '.npz'
        # data_name = data_name + str(AMSTpfa.id) + '.npz'
        # self.id = AMSTpfa.id
        # super().__init__(data_name=data_name, load=load)
        # AMSTpfa.id += 1

        # self.T = T
        self.wirebasket_elements = np.array([internals, faces, edges, vertices])
        self.wirebasket_numbers = np.array([len(internals), len(faces), len(edges), len(vertices)])
        self.nv = self.wirebasket_numbers[-1]
        self.ns_sum = [self.wirebasket_numbers[0]]
        for i in range(3):
            self.ns_sum.append(self.ns_sum[i] + self.wirebasket_numbers[i+1])
        self.ns_sum = np.array(self.ns_sum)

        n_reord = 0
        self.wirebasket_ids = []
        for i in range(4):
            n2 = len(self.wirebasket_elements[i])
            self.wirebasket_ids.append(np.arange(n_reord, n_reord + n2))
            n_reord += n2

        self.wirebasket_ids = np.array(self.wirebasket_ids)
        self.get_G()
        # dt = [('gid', np.dtype(int)), ('primal_id', np.dtype(int))]
        # self.gid_to_primal = np.zeros(len(gids), dtype=dt)
        # self.gid_to_primal['gid'] = gids
        # self.gid_to_primal['primal_id'] = primal_ids
        cols = primal_ids[vertices]
        lines = np.arange(len(cols))
        data = np.ones(len(lines))
        self.G2 = sp.csc_matrix((data,(lines,cols)), shape=(self.nv, self.nv))
        # self.T_wire = self.G*self.T*self.GT
        # # self.T_wire = self.GT*self.T*self.G
        # self.get_as()
        # self._data['OP_AMS_' + str(self.id)] = self.get_OP_AMS_TPFA_by_AS()

        # self.run()

    def get_G(self):
        cols = np.concatenate(self.wirebasket_elements)
        n = len(cols)
        lines = np.concatenate(self.wirebasket_ids)
        data = np.ones(n)
        self.G = sp.csc_matrix((data,(lines,cols)), shape=(n, n))
        self.GT = self.G.copy()
        self.GT = self.GT.transpose()

    def get_as(self, T_wire):

        As = dict()

        Tmod = T_wire.copy().tolil()
        # As['Tf'] = Tmod
        ni = self.wirebasket_numbers[0]
        nf = self.wirebasket_numbers[1]
        ne = self.wirebasket_numbers[2]
        nv = self.wirebasket_numbers[3]

        nni = self.ns_sum[0]
        nnf = self.ns_sum[1]
        nne = self.ns_sum[2]
        nnv = self.ns_sum[3]

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

        As['Aii'] = Aii
        As['Aif'] = Aif
        As['Aff'] = Aff
        As['Afe'] = Afe
        As['Aee'] = Aee
        As['Aev'] = Aev
        As['Ivv'] = Ivv

        return As

    def get_OP_AMS_TPFA_by_AS(self, As):

        ni = self.wirebasket_numbers[0]
        nf = self.wirebasket_numbers[1]
        ne = self.wirebasket_numbers[2]
        nv = self.wirebasket_numbers[3]

        nni = self.ns_sum[0]
        nnf = self.ns_sum[1]
        nne = self.ns_sum[2]
        nnv = self.ns_sum[3]

        lines = np.arange(nne, nnv).astype(np.int32)
        ntot = (self.wirebasket_numbers.sum())
        op = sp.lil_matrix((ntot, nv))
        op[lines] = As['Ivv'].tolil()

        M = As['Aee']
        M = linalg.spsolve(M.tocsc(), sp.identity(ne).tocsc())
        M = M.dot(-1*As['Aev'])
        op[nnf:nne] = M.tolil()

        M2 = As['Aff']
        M2 = linalg.spsolve(M2.tocsc(), sp.identity(nf).tocsc())
        M2 = M2.dot(-1*As['Afe'])
        M = M2.dot(M)
        op[nni:nnf] = M.tolil()

        M2 = As['Aii']
        M2 = linalg.spsolve(M2.tocsc(), sp.identity(ni).tocsc())
        M2 = M2.dot(-1*As['Aif'])
        M = M2.dot(M)
        op[0:nni] = M.tolil()

        return self.GT*op*self.G2

    def run(self, T: 'transmissibility matrix'):

        T_wire = self.G*T*self.GT
        # self._data['T_wire'] = T_wire
        As = self.get_as(T_wire)
        OP = self.get_OP_AMS_TPFA_by_AS(As)
        return OP
