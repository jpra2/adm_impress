import scipy.sparse as sp
import numpy as np

class MatricesForCorrection:

    @staticmethod
    def get_B_matrix(total_source_term, gravity_source_term):

        ###################
        ## teste
        rr = total_source_term==0
        total_source_term2 = total_source_term.copy()
        total_source_term2[rr] = np.ones(len(total_source_term[rr]))
        ###################

        diagg = gravity_source_term/total_source_term2
        n = len(total_source_term)

        B_matrix = sp.lil_matrix((n ,n))
        B_matrix.setdiag(diagg)

        return B_matrix.tocsc()

    @staticmethod
    def get_Eps_matrix(gids, volumes_without_gravity_source_term):

        n = len(gids)
        n2 = len(volumes_without_gravity_source_term)
        diagg = np.ones(n, dtype=int)
        diagg[volumes_without_gravity_source_term] = np.zeros(n2, dtype=int)

        Eps_matrix = sp.lil_matrix((n, n))
        Eps_matrix.setdiag(diagg)

        return Eps_matrix.tocsc()
