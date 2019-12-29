from .ams_tpfa import AMSTpfa
import time
import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg

class AMSMpfa(AMSTpfa):

    def __init__(self,
        internals,
        faces,
        edges,
        vertices,
        gids: 'global_ids',
        primal_ids: 'primal_ids',
        load=False,
        data_name='AMSMpfa_'):

        super().__init__(
            internals,
            faces,
            edges,
            vertices,
            gids,
            primal_ids,
            load=load,
            data_name=data_name
        )

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
        Aie = Tmod[0:nni, nnf:nne]
        Aiv = Tmod[0:nni, nne:nnv]

        #faces
        Aff = Tmod[nni:nnf, nni:nnf]
        Afe = Tmod[nni:nnf, nnf:nne]
        Afv = Tmod[nni:nnf, nne:nnv]

        Aant = Tmod[nni:nnf, 0:nni]
        soma = Aant.sum(axis=1)
        d1 = np.matrix(Aff.diagonal()).reshape([nf, 1])
        d1 += soma
        Aff.setdiag(d1)

        #arestas
        Aee = Tmod[nnf:nne, nnf:nne]
        Aev = Tmod[nnf:nne, nne:nnv]

        Aant = Tmod[nnf:nne, 0:nnf]
        soma = Aant.sum(axis=1)
        d1 = np.matrix(Aee.diagonal()).reshape([ne, 1])
        d1 += soma
        Aee.setdiag(d1)

        Ivv = sp.identity(nv)

        As['Aii'] = Aii
        As['Aif'] = Aif
        As['Aie'] = Aie
        As['Aiv'] = Aiv
        As['Aff'] = Aff
        As['Afe'] = Afe
        As['Afv'] = Afv
        As['Aee'] = Aee
        As['Aev'] = Aev
        As['Ivv'] = Ivv

        return As

    def get_OP_AMS_MPFA_by_AS(self, As):

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

        Meeinv = As['Aee']
        Meeinv = linalg.spsolve(Meeinv.tocsc(), sp.identity(ne).tocsc())
        Pe = Meeinv.dot(-1*As['Aev'])
        op[nnf:nne] = Pe.tolil()

        Mffinv = As['Aff']
        Mffinv = linalg.spsolve(Mffinv.tocsc(), sp.identity(nf).tocsc())
        Pf = As['Afe'].dot(-1*Pe) - As['Afv']
        Pf = Mffinv.dot(Pf)
        op[nni:nnf] = Pf.tolil()

        Miiinv = As['Aii']
        Miiinv = linalg.spsolve(Miiinv.tocsc(), sp.identity(ni).tocsc())
        Pi = As['Aif'].dot(-1*Pf) + As['Aie'].dot(-1*Pe) - As['Aiv']
        Pi = Miiinv.dot(Pi)
        op[0:nni] = Pi.tolil()

        return self.GT*op*self.G2

    def run(self, T: 'transmissibility matrix'):

        T_wire = self.G*T*self.GT
        As = self.get_as(T_wire)
        OP = self.get_OP_AMS_MPFA_by_AS(As)
        return OP
