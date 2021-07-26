from .....data_class.data_manager import DataManager
import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg
from scipy.sparse.linalg import inv
import time

class AMSTpfa:
    # name = 'AMSTpfa_'
    id = 1

    def __init__(self,
        gids,
        dual_flags,
        primal_ids,
        load=False,
        tpfalizar=False,
        get_prolongation_operator=True,
        get_correction_term=False):

        self.get_prolongation_operator = get_prolongation_operator
        self.get_correction_term = get_correction_term
        self.tpfalizar = tpfalizar
        self.wirebasket_elements, self.wirebasket_numbers, self.nv, self.ns_sum, self.wirebasket_ids = get_wirebasket_elements(gids, dual_flags)
        self.get_G()
        self.G2 = get_G2(self.wirebasket_elements[3], primal_ids)

    def run(self, T: 'transmissibility matrix', total_source_term=None, B_matrix=None, Eps_matrix=None, return_correction_matrix=False):

        if self.get_correction_term:
            B_wire = self.G*B_matrix*self.GT
            Eps_wire = self.G*Eps_matrix*self.GT
            self.I = sp.identity(len(np.concatenate(self.wirebasket_elements))).tocsc()
            self.E_wire = Eps_wire*B_wire + self.I - B_wire
            self.it = 0
            del B_wire, Eps_wire

        T_wire = self.G*T*self.GT
        if self.tpfalizar:
            T_wire = self.tpfalize(T_wire)
        # self._data['T_wire'] = T_wire
        As = self.get_as(T_wire)
        del T_wire

        if self.get_prolongation_operator:
            OP = self.get_OP_AMS_TPFA_by_AS(As)
        else:
            OP = sp.lil_matrix((len(np.concatenate(self.wirebasket_elements)), self.G2.shape[0])).tocsc()

        if self.get_correction_term:
            total_source_term_wire = self.G*total_source_term
            pcorr = self.get_pcorr(As, total_source_term_wire)
        else:
            pcorr = np.zeros(len(np.concatenate(self.wirebasket_elements)))

        Cmatrix = self.get_correction_matrix(As, return_correction_matrix=return_correction_matrix)


        return OP.tocsc(), pcorr, Cmatrix

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

    def get_OP_AMS_TPFA_by_AS(self, As):

        nv = self.wirebasket_numbers[3]

        Pv = sp.identity(nv)
        Pe = -linalg.spsolve(As['Aee'].tocsc(),(As['Aev']*Pv).tocsc())
        Pf = -linalg.spsolve(As['Aff'].tocsc(),(As['Afe']*Pe).tocsc())
        Pi = -linalg.spsolve(As['Aii'].tocsc(),(As['Aif']*Pf).tocsc())
        op = sp.vstack([Pi,Pf,Pe,Pv])

        return self.GT*op*self.G2

    def tpfalize(self, T_wire):

        ni = self.wirebasket_numbers[0]
        nf = self.wirebasket_numbers[1]
        ne = self.wirebasket_numbers[2]
        nv = self.wirebasket_numbers[3]

        nni = self.ns_sum[0]
        nnf = self.ns_sum[1]
        nne = self.ns_sum[2]
        nnv = self.ns_sum[3]

        T_wire2 = T_wire.copy().tolil()

        rr = np.array(T_wire2[0:nni, nnf:nnv].sum(axis=1).transpose())[0]
        T_wire2[0:nni,nnf:nnv] = 0
        inds = np.arange(nni)
        T_wire2[inds,inds] += rr

        rr = np.array(T_wire2[nni:nnf, nne:nnv].sum(axis=1).transpose())[0]
        T_wire2[nni:nnf, nne:nnv] = 0
        inds = np.arange(nni, nnf)
        T_wire2[inds,inds] += rr

        rr = np.array(T_wire2[nnf:nne, 0:nni].sum(axis=1).transpose())[0]
        T_wire2[nnf:nne, 0:nni] = 0
        inds = np.arange(nnf, nne)
        T_wire2[inds,inds] += rr

        rr = np.array(T_wire2[nne:nnv, 0:nnf].sum(axis=1).transpose())[0]
        T_wire2[nne:nnv, 0:nnf] = 0
        inds = np.arange(nne, nnv)
        T_wire2[inds,inds] += rr

        return T_wire2

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

        if self.it > 0:
            self.E_wire = self.I

        # q2 = self.E_wire*total_source_term_wire
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

        self.it += 1

        return pcorr

    def get_correction_matrix(self, As, return_correction_matrix=False):
        nf = sum(self.wirebasket_numbers)
        Cmatrix = sp.lil_matrix((nf, nf))

        if not return_correction_matrix:
            return Cmatrix.tocsc()

        ni = self.wirebasket_numbers[0]
        nf = self.wirebasket_numbers[1]
        ne = self.wirebasket_numbers[2]
        nv = self.wirebasket_numbers[3]

        nni = self.ns_sum[0]
        nnf = self.ns_sum[1]
        nne = self.ns_sum[2]

        # Aee_inv = linalg.spsolve(As['Aee'], sp.identity(As['Aee'].shape[0]).tocsc())
        # Aff_inv = linalg.spsolve(As['Aff'], sp.identity(As['Aff'].shape[0]).tocsc())
        # Aii_inv = linalg.spsolve(As['Aii'], sp.identity(As['Aff'].shape[0]).tocsc())

        # Cmatrix[0:nni, 0:nni] = Aii_inv
        # Cmatrix[nni:nnf, nni:nnf] = Aff_inv
        # Cmatrix[nnf:nne, nnf:nne] = Aee_inv
        # Cmatrix[nni:nnf, nnf:nne] = -Aff_inv * As['Afe'] * Aee_inv
        # Cmatrix[0:nni, nni:nnf] = -Aii_inv * As['Aif'] * Aff_inv
        # Cmatrix[0:nni, nnf:nne] = Aii_inv * As['Aif'] * (-Cmatrix[nni:nnf, nnf:nne])

        # inversa = linalg.spsolve(As['Aee'], sp.identity(As['Aee'].shape[0]).tocsc())
        inversa = inv(As['Aee'])
        Cmatrix[nnf:nne, nnf:nne] = inversa
        if nf > 0:
            # inversa = linalg.spsolve(As['Aff'], sp.identity(As['Aff'].shape[0]).tocsc())
            inversa = inv(As['Aff'])
            Cmatrix[nni:nnf, nni:nnf] = inversa
            Cmatrix[nni:nnf, nnf:nne] = -Cmatrix[nni:nnf, nni:nnf] * As['Afe'] * Cmatrix[nnf:nne, nnf:nne]

        if ni > 0:
            # inversa = linalg.spsolve(As['Aii'], sp.identity(As['Aff'].shape[0]).tocsc())
            inversa = inv(As['Aii'])
            Cmatrix[0:nni, 0:nni] = inversa
            Cmatrix[0:nni, nni:nnf] = -Cmatrix[0:nni, 0:nni] * As['Aif'] * Cmatrix[nni:nnf, nni:nnf]
            Cmatrix[0:nni, nnf:nne] = Cmatrix[0:nni, 0:nni] * As['Aif'] * (-Cmatrix[nni:nnf, nnf:nne])

        return self.GT*Cmatrix.tocsc()*self.G

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
