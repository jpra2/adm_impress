from packs import directories as direc
import numpy as np

def get_reservoir_partitions(coord_nodes, external_vertex_on_boundary, uniform_dual=False):
    cr1 = direc.data_loaded['Crs']['Cr1']
    cr2 = direc.data_loaded['Crs']['Cr2']
    # cr1 = np.load('flying/crs.npy').tolist()
    # cr1 = [5, 5, 1]
    crs = [cr1, cr2]

    # import pdb; pdb.set_trace()
    if np.array(cr1).max()<5:
        uniform_dual=False
    Lx, Ly, Lz = coord_nodes.max(axis=0)-coord_nodes.min(axis=0)
    min_x, min_y, min_z = coord_nodes.min(axis=0)
    max_x, max_y, max_z = min_x+Lx, min_y+Ly, min_z+Lz
    min_j=[min_x, min_y, min_z]
    max_j=[max_x, max_y, max_z]

    d_x=Lx/(len(np.unique(np.round(coord_nodes[:,0])))-1)
    d_y=Ly/(len(np.unique(np.round(coord_nodes[:,1])))-1)
    d_z=Lz/(len(np.unique(np.round(coord_nodes[:,2])))-1)
    d_j=[d_x, d_y, d_z]
    P=[]
    D=[]
    for i in range(len(crs)):
        P_i = []
        D_i = []
        for j in range(3):
            if uniform_dual:
                if (max_j[j]-min_j[j])<=crs[i][j]*d_j[j]+1e-10:
                    Pij = np.arange(min_j[j],round(max_j[j])+d_j[j],crs[i][j]*d_j[j])
                else:
                    n_homog_prim=round((max_j[j]-min_j[j])/(crs[i][j]*d_j[j]))-1
                    length_non_homog_prim=max_j[j]-min_j[j]-n_homog_prim*crs[i][j]*d_j[j]
                    initial_homog_primal=min_j[j]+int((length_non_homog_prim/d_j[j])/2)*d_j[j]
                    Pij=np.array([min_j[j]])
                    p_homog=initial_homog_primal+np.cumsum(np.repeat(crs[i][j]*d_j[j],n_homog_prim+2))-crs[i][j]*d_j[j]
                    Pij=np.append(Pij,p_homog)
            else:
                
                Pij = np.arange(min_j[j],round(max_j[j])+d_j[j],crs[i][j]*d_j[j])

            Pij[-1] = max_j[j]

            Dij = (Pij[1:]+Pij[:-1])/2

            # if j==0:
            #     import pdb; pdb.set_trace()
            if external_vertex_on_boundary and len(Dij)>1:
                Dij[0] = min_j[j]+d_j[j]/2
                Dij[-1] = max_j[j]-d_j[j]/2
            # import pdb; pdb.set_trace()
            P_i.append(Pij)
            D_i.append(Dij)
        P.append(P_i)
        D.append(D_i)
    # import pdb; pdb.set_trace()
    return P, D, min_j, max_j, d_j

def distribute_reservoir_partitions(P_all, D_all, nworker):
    P = P_all.copy()
    D = D_all.copy()
    coarser_primal = P[-1]
    max_workers = 0
    partitioning_direction = 0
    for j in range(3):
        if len(coarser_primal[j])-1 > max_workers:
            max_workers = len(coarser_primal[j])-1
            partitioning_direction = j
    if nworker > max_workers:
        nworker = max_workers
        print('More workers than partitions, working now with: {} processes'.format(nworker))

    vector_to_partitionate = P[-1][partitioning_direction]
    number_of_intervals = len(vector_to_partitionate)-1
    intervals_by_worker = int(number_of_intervals/nworker)
    i0=vector_to_partitionate[0:number_of_intervals:intervals_by_worker]
    i1=vector_to_partitionate[intervals_by_worker:number_of_intervals+1:intervals_by_worker]
    intervals=np.array([i0,i1]).T
    intervals[-1][-1]=vector_to_partitionate[-1]

    subP = []
    subD = []
    for interval in intervals:
        sP=[]
        sD=[]
        for l in range(len(P)):
            sPl=[]
            sDl=[]
            for j in range(3):
                if j==partitioning_direction:
                    prim = P[l][j]
                    dual = D[l][j]
                    sPl.append(prim[(prim>=interval[0]) & (prim<=interval[1])])
                    sDl.append(dual[(dual>=interval[0]) & (dual<=interval[1])])
                else:
                    sPl.append(P[l][j])
                    sDl.append(D[l][j])
            sP.append(sPl)
            sD.append(sDl)
        subP.append(sP)
        subD.append(sD)
    return subP, subD

def get_box(centroids, limites):
    l0 = centroids[:,0] > limites[0,0]
    l1 = centroids[:,1] > limites[0,1]
    l2 = centroids[:,2] > limites[0,2]
    lc1 = l0 & l1 & l2
    l0 = centroids[:,0] < limites[1,0]
    l1 = centroids[:,1] < limites[1,1]
    l2 = centroids[:,2] < limites[1,2]
    lc2 = l0 & l1 & l2
    conts = lc1 & lc2
    return conts

def create_dual_and_primal(subP, subD, min_j, max_j, d_j, cent_volumes):
    all_volumes=np.arange(len(cent_volumes))
    min_x, min_y, min_z = min_j
    max_x, max_y, max_z = max_j
    Lx1, Ly1, Lz1 = subP[0]
    Lx2, Ly2, Lz2 = subP[1]
    Lxd1, Lyd1, Lzd1 = subD[0]
    Lxd2, Lyd2, Lzd2 = subD[1]
    primal_2 = []
    primal_1 = []
    dual_flag_1 = []
    dual_flag_2 = []
    for i in range(len(Lx2)-1):
        bx = np.array([[Lx2[i], min_y, min_z],[Lx2[i+1], max_y, max_z]])
        indx=get_box(cent_volumes[all_volumes], bx)
        vx = all_volumes[indx]
        for j in range(len(Ly2)-1):
            by = np.array([[Lx2[i], Ly2[j], min_z],[Lx2[i+1], Ly2[j+1], max_z]])
            indy = get_box(cent_volumes[vx], by)
            vy = vx[indy]
            for k in range(len(Lz2)-1):
                bz = np.array([[Lx2[i], Ly2[j], Lz2[k]],[Lx2[i+1], Ly2[j+1], Lz2[k+1]]])
                indz = get_box(cent_volumes[vy], bz)
                vz = vy[indz]
                primal_2.append(vz)

                ################### Creates dual flags lv2###############
                lxd2 = Lxd2[(Lxd2>=Lx2[i]) & (Lxd2<Lx2[i+1])]
                lyd2 = Lyd2[(Lyd2>=Ly2[j]) & (Lyd2<Ly2[j+1])]
                lzd2 = Lzd2[(Lzd2>=Lz2[k]) & (Lzd2<Lz2[k+1])]
                lzd2=np.unique(lzd2)

                ld2=[lxd2, lyd2, lzd2]
                dual_flag=np.zeros(len(vz))
                for d in range(3):
                    try:
                        dual_flag += (cent_volumes[vz,d] > ld2[d]-d_j[d]/2) & (cent_volumes[vz,d] < ld2[d]+d_j[d]/2)
                    except:
                        import pdb; pdb.set_trace()
                dual_flag_2.append(dual_flag)
                ########################################################

                lx1 = Lx1[(Lx1>=Lx2[i]) & (Lx1<=Lx2[i+1])]
                ly1 = Ly1[(Ly1>=Ly2[j]) & (Ly1<=Ly2[j+1])]
                lz1 = Lz1[(Lz1>=Lz2[k]) & (Lz1<=Lz2[k+1])]

                for l in range(len(lx1)-1):
                    bx = np.array([[lx1[l], Ly2[j], Lz2[k]], [lx1[l+1], Ly2[j+1], Lz2[k+1]]])
                    indx=get_box(cent_volumes[vz], bx)
                    vx1 = vz[indx]
                    for m in range(len(ly1)-1):
                        by = np.array([[lx1[l], ly1[m], Lz2[k]],[lx1[l+1], ly1[m+1], Lz2[k+1]]])
                        indy = get_box(cent_volumes[vx1], by)
                        vy1 = vx1[indy]

                        for n in range(len(lz1)-1):
                            bz = np.array([[lx1[l], ly1[m], lz1[n]],[lx1[l+1], ly1[m+1], lz1[n+1]]])
                            indz = get_box(cent_volumes[vy1], bz)
                            vz1 = vy1[indz]
                            primal_1.append(vz1)

                            lxd1 = Lxd1[(Lxd1>=lx1[l]) & (Lxd1<=lx1[l+1])]
                            lyd1 = Lyd1[(Lyd1>=ly1[m]) & (Lyd1<=ly1[m+1])]
                            lzd1 = Lzd1[(Lzd1>=lz1[n]) & (Lzd1<=lz1[n+1])]
                            lzd1 = np.unique(lzd1)
                            ld1=[lxd1, lyd1, lzd1]
                            dual_flag=np.zeros(len(vz1))
                            for d in range(3):
                                dual_flag += (cent_volumes[vz1,d] > ld1[d]-d_j[d]/2) & (cent_volumes[vz1,d] < ld1[d]+d_j[d]/2)

                            dual_flag_1.append(dual_flag)
    return primal_1, primal_2, dual_flag_1, dual_flag_2

def set_tags(M1, primal_1, primal_2, dual_flag_1, dual_flag_2):
    volumes=np.array(M1.all_volumes)
    p1=np.concatenate(primal_1)
    p2=np.concatenate(primal_2)
    id1=[np.repeat(i,len(primal_1[i])) for i in range(len(primal_1))]
    id2=[np.repeat(i,len(primal_2[i])) for i in range(len(primal_2))]
    id1=np.concatenate(id1)
    id2=np.concatenate(id2)
    M1.mb.tag_set_data(M1.fine_to_primal1_classic_tag,volumes[p1], id1)
    M1.mb.tag_set_data(M1.fine_to_primal2_classic_tag,volumes[p2], id2)
    d1=np.concatenate(dual_flag_1).astype(int)
    d2=np.concatenate(dual_flag_2).astype(int)
    M1.mb.tag_set_data(M1.D1_tag,volumes[p1], d1)
    M1.mb.tag_set_data(M1.D2_tag,volumes[p2], d2)

    for i in range(len(primal_1)):
        ms=M1.mb.create_meshset()
        M1.mb.add_entities(ms,volumes[primal_1[i]])
        M1.mb.tag_set_data(M1.primal_id_tag1,ms,i)

    for i in range(len(primal_2)):
        ms=M1.mb.create_meshset()
        M1.mb.add_entities(ms,volumes[primal_2[i]])
        M1.mb.tag_set_data(M1.primal_id_tag2,ms,i)

    # ms=M1.mb.create_meshset()
    # M1.mb.add_entities(ms,M1.all_volumes)
    # M1.mb.write_file("results/dual_test.vtk",[ms])
    # import pdb; pdb.set_trace()
class DualPrimal:
    def __init__(self, M1, coord_nodes, cent_volumes, external_vertex_on_boundary=True):

        P, D, min_j, max_j, d_j = get_reservoir_partitions(coord_nodes, external_vertex_on_boundary, uniform_dual=True)

        # subP, subD = distribute_reservoir_partitions(P, D, nworker=3)
        primal_1, primal_2, dual_flag_1, dual_flag_2 = create_dual_and_primal(P, D, min_j, max_j, d_j, cent_volumes)
        set_tags(M1, primal_1, primal_2, dual_flag_1, dual_flag_2)
        # m=M1.mb.create_meshset()
        # M1.mb.add_entities(m,M1.all_volumes)
        # M1.mb.write_file('results/trashs.vtk',[m])
        # import pdb; pdb.set_trace()
