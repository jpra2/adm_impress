def get_local_matrices_and_global_ids(separated_dual_structures, ks):
    conjs_matrices=[]

    for separated_dual_structure in separated_dual_structures:
        if len(separated_dual_structure[0])>0:
            
            ss=separated_dual_structure
            line, cols, faces, vols =ss[0], ss[1], ss[2], ss[3]
            for line, cols, faces, vols in zip(line, cols, faces, vols):
                data=np.concatenate(ks, ks, -ks, -ks)
                mat=csc_matrix((data,(lines,cols)),shape=(len(vols),len(vols)))
                local_matrices.append(mat)
                global_ids.append(vols)

            conjs_matrices.append([local_matrices, global_gids])
        else:
            conjs_matrices.append([])
    return cnp.array(conjs_matrices, global_ids)
