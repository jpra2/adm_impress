from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu as lu
import numpy as np
import pdb

def get_local_lu_and_global_ids(separated_dual_structures, ts):
    local_lu_and_gids=[]
    
    for separated_dual_structure in separated_dual_structures:
        if len(separated_dual_structure[0])>0:
            local_lu=[]
            global_ids=[]
            ss=separated_dual_structure
            lines, cols, faces, vols =ss#ss[0], ss[1], ss[2], ss[3]
            for l, c, f, v in zip(lines, cols, faces, vols):
                t=ts[f]
                d=np.zeros_like(f,dtype=float)
                d[f>=0]=ts[f[f>=0]]
                d[f<0]=-ts[-f[f<0]-1] # if f<0, the corresponding face is -f-1 (convention for unicity)
                mat=csc_matrix((d,(l,c)),shape=(len(v),len(v)))
                local_lu.append(lu(mat))

                global_ids.append(v)

            local_lu_and_gids.append([local_lu, global_ids])
        else:
            local_lu_and_gids.append([])
    return local_lu_and_gids

class LocalLU:
    def update_lu_objects(self, separated_dual_structures, ts):
        '''
        separated_dual_structure: preprocessed dual clusters from each dual class
        ts: transmissibility of all faces
        '''
        self.local_lu_and_global_ids=get_local_lu_and_global_ids(separated_dual_structures, ts)
