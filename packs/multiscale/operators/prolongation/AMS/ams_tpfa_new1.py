from .....data_class.data_manager import DataManager
import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg
import time
import pdb

class AMSTpfa:

    def __init__(self,
        gids,
        primal_ids,
        dual_flags,
        proc=-1,
        locker=False):

        self.proc = proc
        self.wirebasket_elements, self.wirebasket_numbers, self.nv, self.ns_sum, self.wirebasket_ids = get_wirebasket_elements(gids, dual_flags)
        self.locker = locker
        self.get_G()
        self.G2 = get_G2(self.wirebasket_elements[3], primal_ids)

    def run(self, T: 'transmissibility matrix', total_source_term=np.array([])):

        assert self.G.shape == T.shape

        try:
            T_wire = self.get_Twire(T)
        except Exception as e:
            with self.locker:
                print()
                print(self.G.shape)
                print(T.shape)
                print(self.GT.shape)
                raise e
        As = self.get_as(T_wire)
        del T_wire

        OP = self.get_OP_AMS_TPFA_by_AS(As)
        if len(total_source_term) == 0:
            pcorr = np.zeros(OP.shape[0])
        else:
            total_source_term_wire = self.G*total_source_term
            pcorr = self.get_pcorr(As, total_source_term_wire)

        return OP, pcorr

    def get_G(self):
        cols = np.concatenate(self.wirebasket_elements)
        n = len(cols)
        lines = np.concatenate(self.wirebasket_ids)
        data = np.ones(n)

        try:
            self.G = sp.csc_matrix((data,(lines,cols)), shape=(n, n))
        except Exception as e:
            if self.locker == False:
                raise e
            else:
                with self.locker:
                    print()
                    print(cols)
                    print(lines)
                    print(len(lines))
                    print(len(cols))
                    print(n)
                    print()
                raise e

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
        soma = Tmod[nni:nnf, 0:nni].sum(axis=1)
        d1 = np.matrix(Aff.diagonal()).reshape([nf, 1])
        d1 += soma
        Aff.setdiag(d1)

        #arestas
        Aee = Tmod[nnf:nne, nnf:nne]
        Aev = Tmod[nnf:nne, nne:nnv]
        soma = Tmod[nnf:nne, nni:nnf].sum(axis=1)
        d1 = np.matrix(Aee.diagonal()).reshape([ne, 1])
        d1 += soma
        Aee.setdiag(d1)
        Ivv = sp.identity(nv)

        As['Aii'] = Aii.tocsc()
        As['Aif'] = Aif.tocsc()
        As['Aff'] = Aff.tocsc()
        As['Afe'] = Afe.tocsc()
        As['Aee'] = Aee.tocsc()
        As['Aev'] = Aev.tocsc()
        As['Ivv'] = Ivv.tocsc()

        return As

    def get_as_off_diagonal(self, T_wire):

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
        # Aii = Tmod[0:nni, 0:nni]
        Aif = Tmod[0:nni, nni:nnf]
        Aie = Tmod[0:nni, nnf:nne]

        #faces
        # Aff = Tmod[nni:nnf, nni:nnf]
        Afe = Tmod[nni:nnf, nnf:nne]
        # soma = Tmod[nni:nnf, 0:nni].sum(axis=1)
        # d1 = np.matrix(Aff.diagonal()).reshape([nf, 1])
        # d1 += soma
        # Aff.setdiag(d1)

        #arestas
        # Aee = Tmod[nnf:nne, nnf:nne]
        # Aev = Tmod[nnf:nne, nne:nnv]
        # soma = Tmod[nnf:nne, nni:nnf].sum(axis=1)
        # d1 = np.matrix(Aee.diagonal()).reshape([ne, 1])
        # d1 += soma
        # Aee.setdiag(d1)
        # Ivv = sp.identity(nv)

        # As['Aii'] = Aii.tocsc()
        As['Aif'] = Aif.tocsc()
        As['Aie'] = Aie.tocsc()
        # As['Aff'] = Aff.tocsc()
        As['Afe'] = Afe.tocsc()
        # As['Aee'] = Aee.tocsc()
        # As['Aev'] = Aev.tocsc()
        # As['Ivv'] = Ivv.tocsc()

        return As

    def get_OP_AMS_TPFA_by_AS(self, As):

        nv = self.wirebasket_numbers[3]

        Pv = sp.identity(nv)
        Pe = -linalg.spsolve(As['Aee'].tocsc(),(As['Aev']*Pv).tocsc())
        Pf = -linalg.spsolve(As['Aff'].tocsc(),(As['Afe']*Pe).tocsc())
        Pi = -linalg.spsolve(As['Aii'].tocsc(),(As['Aif']*Pf).tocsc())
        op = sp.vstack([Pi,Pf,Pe,Pv])

        return self.GT*op*self.G2

    def get_pcorr(self, As, total_source_term_wire):

        '''
        obtem a correcao da pressao
        '''

        ni = self.wirebasket_numbers[0]
        nf = self.wirebasket_numbers[1]
        ne = self.wirebasket_numbers[2]
        nv = self.wirebasket_numbers[3]

        nni = self.ns_sum[0]
        nnf = self.ns_sum[1]
        nne = self.ns_sum[2]
        nnv = self.ns_sum[3]

        q2 = total_source_term_wire
        pcorr = np.zeros(len(q2), dtype=float)

        pcorr_ee = linalg.spsolve(As['Aee'], q2[nnf:nne])
        pcorr_fe = -linalg.spsolve(As['Aff'], As['Afe']*pcorr_ee)
        pcorr_ie = -linalg.spsolve(As['Aii'], As['Aif']*pcorr_fe)

        pcorr_ff = linalg.spsolve(As['Aff'], q2[nni:nnf])
        pcorr_if = -linalg.spsolve(As['Aii'], As['Aif']*pcorr_ff)

        pcorr_ii = linalg.spsolve(As['Aii'], q2[0:nni])

        pcorr[0:nni] = pcorr_ii + pcorr_if + pcorr_ie
        pcorr[nni:nnf] = pcorr_ff + pcorr_fe
        pcorr[nnf:nne] = pcorr_ee
        pcorr = self.GT*pcorr

        return pcorr

    def get_Twire(self, T):
        T_wire = self.G*T*self.GT
        return T_wire



def get_wirebasket_elements(gids, dual_flags):

    internals = gids[dual_flags==0]
    faces = gids[dual_flags==1]
    edges = gids[dual_flags==2]
    vertices = gids[dual_flags==3]

    wirebasket_elements = np.array([internals, faces, edges, vertices])
    wirebasket_numbers = np.array([len(internals), len(faces), len(edges), len(vertices)])
    nv = wirebasket_numbers[-1]
    ns_sum = [wirebasket_numbers[0]]
    for i in range(3):
        ns_sum.append(ns_sum[i] + wirebasket_numbers[i+1])
    ns_sum = np.array(ns_sum)

    n_reord = 0
    wirebasket_ids = []
    for i in range(4):
        n2 = wirebasket_numbers[i]
        wirebasket_ids.append(np.arange(n_reord, n_reord + n2))
        n_reord += n2

    wirebasket_ids = np.array(wirebasket_ids)

    return wirebasket_elements, wirebasket_numbers, nv, ns_sum, wirebasket_ids

def get_G2(vertices, primal_ids):

    nv = len(vertices)
    cols = primal_ids[vertices]
    lines = np.arange(len(cols))
    data = np.ones(len(lines))
    G2 = sp.csc_matrix((data,(lines,cols)), shape=(nv, nv))

    return G2
