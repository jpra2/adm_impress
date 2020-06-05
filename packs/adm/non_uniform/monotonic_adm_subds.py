import numpy as np
import scipy.sparse as sp

def get_preprossed_monotonic_primal_objects(data_impress, T, OP_AMS):
    lcd=sp.find(T)
    lines=lcd[0]
    cols=lcd[1]
    data=lcd[2]
    OP_copy=OP_AMS.tolil()
    gids0=data_impress['GID_0']
    gids1=data_impress['GID_1']
    lines_1=gids1[lines]
    cols_1=gids1[cols]
    preprocessed_primal_objects=[]
    for gid1 in np.unique(gids1):
        pos_lines1 = lines_1 == gid1
        pos_cols1 = cols_1 == gid1
        l0=lines[pos_lines1]
        c0=cols[pos_cols1]
        T_l=T_copy[pos_coarse][pos_coarse]
        OP_l=OP_copy[pos_coarse]
        preprocessed_primal_objects.append(PrepMonotonicPrimal())


class PrepMonotonicPrimal:
    def __init__(self, data_impress, OP_AMS):
        self.primal_monotonic_subds(data_impress,OP_AMS)

    def create_primal_subds(data_impress, OP_AMS):

        import pdb; pdb.set_trace()

class MonotonicPrimal:
    def __init__(self,prep_monotonic_primals):
        a=1
