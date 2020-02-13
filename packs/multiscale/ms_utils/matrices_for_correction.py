import scipy.sparse as sp
import numpy as np

class MatricesForCorrection:

    @staticmethod
    def get_B_matrix(total_source_term, gravity_source_term):

        ###################
        ## teste
        # rr = total_source_term==0
        # total_source_term[rr] = np.ones(len(total_source_term[rr]))
        ###################
        diagg=np.ones(len(total_source_term))
        non_null=abs(total_source_term)>1e-10
        if len(non_null) > 0:
            # import pdb; pdb.set_trace()
            diagg[non_null]=gravity_source_term[non_null]/total_source_term[non_null]

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
