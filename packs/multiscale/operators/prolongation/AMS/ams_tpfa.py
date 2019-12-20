from .....data_class.data_manager import DataManager
import numpy as np
import scipy.sparse as sp

class AMSTpfa:
    name = 'AMSTpfa_'
    id = 1

    def __init__(self, T: 'transmissibility_matrix',
        internals,
        faces,
        edges,
        vertices,
        load=False):

        self.T = T
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

        # data_name = AMSTpfa.name + str(AMSTpfa.id)
        # super().__init__(data_name=data_name, load=load)
        # AMSTpfa.id += 1

        self.run()

    def get_G(self):
        cols = self.wirebasket_elements.flatten()
        lines = self.wirebasket_ids.flatten()
        data = np.ones(len(cols))
        self.G = sp.csc_matrix((data,(lines,cols)),shape=(self.nv, self.nv))
        self.GT = self.G.copy()
        self.GT = self.GT.transpose()

    def get_as(self):

        Tmod = self.T_wire.copy().tolil()
        # As['Tf'] = Tmod
        ni = self.wirebasket_numbers[0]
        nf = self.wirebasket_numbers[1]
        ne = self.wirebasket_numbers[2]
        nv = self.wirebasket_numbers[3]

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

        self.As['Aii'] = Aii
        self.As['Aif'] = Aif
        self.As['Aff'] = Aff
        self.As['Afe'] = Afe
        self.As['Aee'] = Aee
        self.As['Aev'] = Aev
        self.As['Ivv'] = Ivv

    def run(self):

        self.T_wire = self.G*self.T*self.GT
        self.get_as()
