from scipy.sparse import csc_matrix
class CorrectionFunction:
    def __init__(self,q):
        self.cf = calculate_correction_function(q)
        matrices_and_gids = self.get_local_matrices_and_global_ids(separated_dual_structure, ks)

    # def calculate_correction_function():
    #     pass

    def get_local_matrices_and_global_ids(separated_dual_structure, ks):
        conjs_matrices={}
        for dual_type in separated_dual_structure.keys()
            adjs = separated_dual_structure[dual_type][0]
            faces = separated_dual_structure[dual_type][1]
            local_matrices=[]
            global_ids=[]
            for adjs, faces in zip(edge_adjs, edge_faces):
                vols=np.unique(np.concatenate(adjs))
                mape=np.arange(vols.max()+1)
                mape[vols]=np.arange(len(vols))
                ad=mape[adjs]
                lines=np.concatenate(ad[:,0], ad[:,1], ad[:,0], ad[:,1])
                cols=np.concatenate(ad[:,1], ad[:,0], ad[:,0], ad[:,1])
                data=np.concatenate(ks, ks, -ks, -ks)
                mat=csc_matrix((data,(lines,cols)),shape=(len(vols),len(vols)))
                local_matrices.append(mat)
                global_ids.append(vols)
            conjs_matrices.append([local_matrices, global_gids])
        return conjs_matrices
