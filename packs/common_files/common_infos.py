import copy
import numpy as np

class CommonInfos:

    def get_local_t(self, T, volumes):
        '''
            Retorna a Matriz de transmissibilidade local
        '''

        T2 = T[volumes][:,volumes]
        data = np.array(T2.sum(axis=1).transpose())[0]
        data2 = T2.diagonal()
        data2 -= data
        T2.setdiag(data2)

        return T2

    def get_local_matrix(self, T, volumes):

        T2 = T[volumes][:,volumes]
        return T2

    def copy(self):
        return copy.deepcopy(self)
