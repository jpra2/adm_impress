import numpy as np
from pymoab import core, types, rng, topo_util
import multiprocessing as mp
import time

def get_box(conjunto, all_centroids, limites, return_inds):
    l0=all_centroids[:,0] > limites[0,0]
    l1=all_centroids[:,1] > limites[0,1]
    l2=all_centroids[:,2] > limites[0,2]
    lc1=l0 & l1 & l2
    l0=all_centroids[:,0] < limites[1,0]
    l1=all_centroids[:,1] < limites[1,1]
    l2=all_centroids[:,2] < limites[1,2]
    lc2=l0 & l1 & l2
    conts= lc1 & lc2
    if return_inds:
        return (rng.Range(np.array(conjunto)[conts]),np.arange(len(conts))[conts].tolist())
    else:
        return rng.Range(np.array(conjunto)[conts])

class dual_and_primal:
    def __init__(self, subDomain):

        M1 = subDomain.M1
        input_dual = subDomain.input_dual_and_primal
        partition = subDomain.partition
        try:
            master = subDomain.master
            paralell = True
            id = subDomain.id
            print("Process {} started".format(id))
        except:
            paralell = False
            master = None
            print("Running serial")

        t0=time.time()
        self.create(M1, input_dual, partition, master, paralell)
        if paralell:
            print("Process {} finished after {} seconds".format(id, time.time()-t0))
        # print("Tempo total para criação da primal e dual: {} segundos".format(time.time()-t0))

    def create(self, M1, input_dual_and_primal, partition, master, paralell):

        lx1, ly1, lz1 = partition.l1

        lxd1, lyd1, lzd1 = partition.ld1
        lx2, ly2, lz2 = partition.l2
        lxd2, lyd2, lzd2 = partition.ld2
        dx0, dy0, dz0 = input_dual_and_primal.d0
        l1, l2 = input_dual_and_primal.l
        D_x, D_y, D_z = input_dual_and_primal.D
        Lx, Ly, Lz = input_dual_and_primal.L
        xmin, xmax, ymin, ymax, zmin, zmax = input_dual_and_primal.lims
        nc1, nc2 = input_dual_and_primal.ncs
        L2_meshset = input_dual_and_primal.meshsets
        v1 = []
        v2 = []

        volumes_1 = np.array([],dtype=np.uint64)
        volumes_2 = np.array([],dtype=np.uint64)
        dual_1 = np.array([],dtype=int)
        dual_2 = np.array([],dtype=int)
        primal_1 = np.array([],dtype=int)
        primal_2 = np.array([],dtype=int)
        print(lx2)
        xref=lx2[0]
        yref=ly2[0]
        zref=lz2[0]
        for i in range(len(lx2)-1):
            t1=time.time()
            if i==len(lx2)-2:
                sx=D_x
            sy=0

            #################################################
            x0=lx2[i]
            x1=lx2[i+1]
            box_x=np.array([[x0-0.01,ymin,zmin],[x1+0.01,ymax,zmax]])
            vols_x=get_box(M1.all_volumes, M1.all_centroids, box_x, False)
            x_centroids=np.array([M1.mtu.get_average_position([v]) for v in vols_x])
            ######################################

            for j in range(len(ly2)-1):
                if j==len(ly2)-2:
                    sy=D_y
                sz=0
                #########################
                y0=ly2[j]
                y1=ly2[j+1]
                box_y=np.array([[x0-0.01,y0-0.01,zmin],[x1+0.01,y1+0.01,zmax]])
                vols_y=get_box(vols_x, x_centroids, box_y, False)
                y_centroids=np.array([M1.mtu.get_average_position([v]) for v in vols_y])
                ###############
                for k in range(len(lz2)-1):
                    if k==len(lz2)-2:
                        sz=D_z
                    ########################################
                    z0=lz2[k]
                    z1=lz2[k+1]
                    tb=time.time()
                    box_dual_1=np.array([[x0-0.01,y0-0.01,z0-0.01],[x1+0.01,y1+0.01,z1+0.01]])
                    vols=get_box(vols_y, y_centroids, box_dual_1, False)
                    ####################
                    l2_meshset=M1.mb.create_meshset()
                    cont=0
                    elem_por_L2=vols
                    v2.append(vols)
                    M1.mb.add_entities(l2_meshset,elem_por_L2)
                    centroid_p2=np.array([M1.mtu.get_average_position([np.uint64(v)]) for v in elem_por_L2])

                    cx,cy,cz=centroid_p2[:,0],centroid_p2[:,1],centroid_p2[:,2]
                    posx=(abs(cx-lxd2[i])<=l1[0]/1.9)
                    posy=(abs(cy-lyd2[j])<=l1[1]/1.9)
                    posz=(abs(cz-lzd2[k])<=l1[2]/1.9)
                    f1a2v3=np.zeros(len(elem_por_L2),dtype=int)
                    f1a2v3[posx]+=1
                    f1a2v3[posy]+=1
                    f1a2v3[posz]+=1

                    volumes_2 = np.append(volumes_2,elem_por_L2)
                    dual_2 = np.append(dual_2,f1a2v3)
                    primal_2 = np.append(primal_2, np.repeat(nc2,len(elem_por_L2)))

                    # M1.mb.tag_set_data(M1.D2_tag, elem_por_L2, f1a2v3)
                    # M1.mb.tag_set_data(M1.fine_to_primal2_classic_tag, elem_por_L2, np.repeat(nc2,len(elem_por_L2)))
                    M1.mb.add_child_meshset(L2_meshset,l2_meshset)

                    sg=M1.mb.get_entities_by_handle(l2_meshset)
                    # print(k, len(sg), time.time()-t1)
                    t1=time.time()
                    M1.mb.tag_set_data(M1.primal_id_tag2, l2_meshset, nc2)

                    centroids_primal2=np.array([M1.mtu.get_average_position([np.uint64(v)]) for v in elem_por_L2])
                    nc2+=1
                    s1x=0
                    vs1=[]
                    for m in range(len(lx1)):
                        a=int(l2[0]/l1[0])*i+m
                        if Lx-D_x==lx2[i]+lx1[m]+l1[0]:# and D_x==Lx-int(Lx/l1[0])*l1[0]:
                            s1x=D_x
                        s1y=0
                        for n in range(len(ly1)):
                            b=int(l2[1]/l1[1])*j+n
                            if Ly-D_y==ly2[j]+ly1[n]+l1[1]:# and D_y==Ly-int(Ly/l1[1])*l1[1]:
                                s1y=D_y
                            s1z=0

                            for o in range(len(lz1)):
                                c=int(l2[2]/l1[2])*k+o
                                if Lz-D_z==lz2[k]+lz1[o]+l1[2]:
                                    s1z=D_z
                                l1_meshset=M1.mb.create_meshset()
                                box_primal1 = np.array([np.array([lx2[i]+lx1[m], ly2[j]+ly1[n], lz2[k]+lz1[o]]), np.array([lx2[i]+lx1[m]+l1[0]+s1x, ly2[j]+ly1[n]+l1[1]+s1y, lz2[k]+lz1[o]+l1[2]+s1z])])
                                elem_por_L1 = get_box(elem_por_L2, centroids_primal2, box_primal1, False)
                                vs1.append(np.array(elem_por_L1))
                                M1.mb.add_entities(l1_meshset,elem_por_L1)

                                centroid_p1=np.array([M1.mtu.get_average_position([np.uint64(v)]) for v in elem_por_L1])
                                cx,cy,cz=centroid_p1[:,0],centroid_p1[:,1],centroid_p1[:,2]
                                try:
                                    posx=(abs(cx-lxd1[a])<=dx0/2)
                                except:
                                    import pdb; pdb.set_trace()
                                posy=(abs(cy-lyd1[b])<=dy0/2)
                                posz=(abs(cz-lzd1[c])<=dz0/2)
                                f1a2v3=np.zeros(len(elem_por_L1),dtype=int)
                                f1a2v3[posx]+=1
                                f1a2v3[posy]+=1
                                f1a2v3[posz]+=1
                                volumes_1 = np.append(volumes_1,elem_por_L1)
                                dual_1 = np.append(dual_1,f1a2v3)
                                primal_1 = np.append(primal_1, np.repeat(nc1,len(elem_por_L1)))

                                M1.mb.tag_set_data(M1.primal_id_tag1, l1_meshset, nc1)
                                nc1+=1
                                M1.mb.add_child_meshset(l2_meshset,l1_meshset)
                    v1.append(np.array(vs1))

        if paralell:
            master.send([volumes_1,dual_1,primal_1, volumes_2, dual_2, primal_2, np.array(v1), np.array(v2)])
        else:
            M1.mb.tag_set_data(M1.D1_tag, volumes_1,dual_1)
            M1.mb.tag_set_data(M1.fine_to_primal1_classic_tag, volumes_1, primal_1)
            M1.mb.tag_set_data(M1.D2_tag, volumes_2,dual_2)
            M1.mb.tag_set_data(M1.fine_to_primal2_classic_tag, volumes_2, primal_2)

class SubD:
    def __init__(self, l1, l2, ld1, ld2):
        self.l1 = l1
        self.ld1 = ld1
        self.l2 = l2
        self.ld2 = ld2

class SubDomain():
    def __init__(self, input_dual, partition, M1, master, id):
        self.input_dual_and_primal = input_dual
        self.partition = partition
        self.M1 = M1
        self.master = master
        self.id = id

class paralell_dual_and_primal:
    def __init__(self, subDomain,nworker, first):
        self.v1, self.v2 = self.run_paralell(subDomain,nworker, first)

    def run_paralell(self, subDomain, nworker, first):
        input_dual=subDomain.input_dual_and_primal
        partition=subDomain.partition
        M1=subDomain.M1

        if first:
            first=False
            print("Main program started")

            l1 = partition.l1
            ld1 = partition.ld1
            lx2, ly2, lz2 = partition.l2
            lxd2, lyd2, lzd2 = partition.ld2
            lxd1, lyd1, lzd1 = partition.ld1

            if nworker>len(lxd2):
                nworker=len(lxd2)
                print("Mais processos que partições, usando agora {} processos".format(len(lxd2)))

            nn = int(len(lxd2)/nworker)
            master2worker = [mp.Pipe() for p in range(nworker)]

            m2w, w2m = list(zip(*master2worker))


            subDomains = []
            for i in range(nworker):
                if i==nworker-1:
                    xd1 = lxd1[i*nn*3:]
                    xd2 = lxd2[i*nn:]
                    x2 = lx2[i*nn:]
                else:
                    xd1 = lxd1[i*nn*3:(i+1)*nn*3]
                    xd2 = lxd2[i*nn:(i+1)*nn]
                    x2 = lx2[i*nn:(i+1)*nn+1]
                ld1 = xd1, lyd1, lzd1
                l2 = x2, ly2, lz2
                ld2 = xd2, lyd2, lzd2
                sub_partition = SubD(l1,l2,ld1,ld2)
                subDomains.append(SubDomain(input_dual,sub_partition, M1, w2m[i], i))
            
            procs = [mp.Process(target=dual_and_primal, args=[s]) for s in subDomains]
            tp=time.time()
            for p in procs:
                p.start()

            volumes_1 = np.array([],dtype=np.uint64)
            volumes_2 = np.array([],dtype=np.uint64)
            dual_1 = np.array([],dtype=int)
            dual_2 = np.array([],dtype=int)
            primal_1 = np.array([],dtype=int)
            primal_2 = np.array([],dtype=int)
            L2_meshsets = np.array([],dtype=np.uint64)
            p1=0
            p2=0
            v1 = []
            v2 = []
            for m in m2w:
                msg=m.recv()

                ''' Trata exceção devido ao tipo de dado: apenas uma dual por processo
                gera numpy array, mais que uma dual gera lista de numpy arrays'''

                try:
                    volumes_1 = np.append(volumes_1, np.concatenate(msg[0]))
                    dual_1 = np.append(dual_1, np.concatenate(msg[1]))
                    primal_1= np.append(primal_1, np.concatenate(msg[2])+p1)

                    volumes_2 = np.append(volumes_2, np.concatenate(msg[3]))
                    dual_2 = np.append(dual_2, np.concatenate(msg[4]))
                    primal_2 = np.append(primal_2, np.concatenate(msg[5])+p2)

                except(ValueError):
                    volumes_1 = np.append(volumes_1, np.concatenate([msg[0]]))
                    dual_1 = np.append(dual_1, np.concatenate([msg[1]]))
                    primal_1= np.append(primal_1, np.concatenate([msg[2]])+p1)

                    volumes_2 = np.append(volumes_2, np.concatenate([msg[3]]))
                    dual_2 = np.append(dual_2, np.concatenate([msg[4]]))
                    primal_2 = np.append(primal_2, np.concatenate([msg[5]])+p2)
                    v1.append(msg[6])
                    v2.append(msg[7])

                try:
                    p1=primal_1.max()+1
                    p2=primal_2.max()+1

                except:
                    p1=np.array([])
                    p2=np.array([])


            M1.mb.tag_set_data(M1.D1_tag, volumes_1,dual_1)
            M1.mb.tag_set_data(M1.fine_to_primal1_classic_tag, volumes_1, primal_1)
            M1.mb.tag_set_data(M1.D2_tag, volumes_2,dual_2)
            M1.mb.tag_set_data(M1.fine_to_primal2_classic_tag, volumes_2, primal_2)

            for p in procs:
                p.join()
            print("Total time with {} processes, {} seconds".format(nworker,time.time()-tp))
            return v1, v2
