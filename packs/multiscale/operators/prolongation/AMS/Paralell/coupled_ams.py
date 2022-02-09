from pymoab import rng, types
import scipy
import numpy as np
import multiprocessing as mp
from scipy.sparse import csc_matrix,csr_matrix, linalg, vstack, find
import time
from packs.multiscale.operators.prolongation.AMS.Paralell.partitionating_parameters import calibrate_partitioning_parameters
import pdb

class DualDomain:
    def __init__(self, data_impress, elements_lv0, volumes, local_couple=0, couple_bound=False):
        self.local_couple=local_couple
        self.couple_bound=couple_bound
        self.adjs, self.ks, self.ids_globais_vols, self.ns =  self.get_local_informations(data_impress, elements_lv0, volumes, local_couple=local_couple, couple_bound=couple_bound)
        self.coarse_ids = data_impress['GID_1'][self.vertices]
        self.A_b_t=[]

    def get_local_informations(self, data_impress, elements_lv0, volumes, local_couple=0, couple_bound=False):
        viz=np.unique(np.concatenate(elements_lv0['volumes_face_volumes'][volumes]))
        nv=len(volumes)
        faces_entities = np.unique(np.concatenate(elements_lv0['volumes_face_faces'][volumes]))
        int_facs=np.setdiff1d(faces_entities, elements_lv0['boundary_faces'])

        so_viz=np.setdiff1d(viz,volumes)
        if len(so_viz) > 0:
            so_viz_faces=np.unique(np.concatenate(elements_lv0['volumes_face_faces'][so_viz]))
            int_facs=np.setdiff1d(int_facs,so_viz_faces)

        dual_id_volumes = data_impress['DUAL_1'][volumes]
        if local_couple>0:
            dual_flags=np.repeat(-1,len(elements_lv0["volumes"]))
            dual_flags[volumes]=dual_id_volumes
            if couple_bound:
                reduce_flag=volumes
                volumes_red=volumes
                dual_flags_red=dual_flags[reduce_flag]
                dual_flags_red[dual_flags_red==2]=1
                if local_couple==2:
                    dual_flags_red[dual_flags_red==1]=0
            else:
                try:
                    reduce_flag = np.setdiff1d(volumes, np.concatenate(elements_lv0['volumes_face_volumes'][so_viz]))
                except:
                    reduce_flag = volumes
                volumes_red = reduce_flag

                dual_flags_red=dual_flags[reduce_flag]
                dual_flags_red[dual_flags_red==2]-=1
                if local_couple==2:
                    dual_flags_red[dual_flags_red==1]=0


            # M.mb.tag_set_data(M.D1_tag,volumes_red,dual_flags_red)
            data_impress['DUAL_1'][volumes_red] = dual_flags_red
            dual_id_volumes = data_impress['DUAL_1'][volumes]

        vertices = volumes[dual_id_volumes==3]
        edges = volumes[dual_id_volumes==2]
        faces = volumes[dual_id_volumes==1]
        internals = volumes[dual_id_volumes==0]

        nv=len(vertices)
        ne=len(edges)
        nf=len(faces)
        ni=len(internals)
        ns=[nv,ne,nf,ni]

        self.vertices = vertices

        adjs = elements_lv0['faces_face_volumes'][int_facs]
        adjs = np.concatenate(adjs).reshape((adjs.shape[0], 2))

        map_l=np.zeros(adjs.max()+1)
        map_l[internals]=np.arange(ni)
        map_l[faces]=np.arange(ni,ni+nf)
        map_l[edges]=np.arange(ni+nf,ni+nf+ne)
        map_l[vertices]=np.arange(ni+nf+ne, ni+nf+ne+nv)

        adjs_l0=map_l[adjs[:,0]]
        adjs_l1=map_l[adjs[:,1]]

        adjs=np.array([adjs_l0, adjs_l1]).T
        # ks=M.mb.tag_get_data(M.k_eq_tag,np.uint64(int_facs)[vv],flat=True)
        ks=data_impress['transmissibility'][int_facs]
        # ids_globais_vols=M.mb.tag_get_data(M.ID_reordenado_tag,np.concatenate([np.uint64(internals),np.uint64(faces), np.uint64(edges),vertices]),flat=True)
        ids_globais_vols=np.concatenate([np.uint64(internals),np.uint64(faces), np.uint64(edges),vertices])

        if nv>4:
            print(nv)
            data_impress['coupled_flag'][volumes]=np.repeat(1.0,len(volumes))
        return adjs, ks, ids_globais_vols, ns

class OP_local:
    def __init__(self, sub_d):
        # self.Nvols=sub_d.Nvols
        # self.Nverts=sub_d.Nvert
        self.OP = self.get_OP(sub_d)

    def get_submatrix(self,id0, id1, ks, slice):
        id0, id1, ks= np.concatenate([id0,id1]), np.concatenate([id1,id0]), np.concatenate([ks, ks])
        (xi, xs, yi, ys)=slice
        inds =(id0>=xi) & (id0<xs) & (id1>=yi) & (id1<ys)

        l1=id0[inds]-xi
        c1=id1[inds]-yi
        d1=ks[inds]
        if xi==yi:
            inds_sup0=((id0>=xi) & (id0<xs) & (id1>=ys))
            ls0=id0[inds_sup0]-xi
            cs0=ls0
            ds0=-ks[inds_sup0]

            l=np.concatenate([l1,l1, ls0])
            c=np.concatenate([c1,l1, cs0])
            d=np.concatenate([d1,-d1, ds0])
        else:
            l=l1
            c=c1
            d=d1
        submatrix=csc_matrix((d,(l,c)),shape=(xs-xi,ys-yi))
        return(submatrix)

    def get_OP(self, sub_d):
        adjs, ks, ids_globais_vols, ns = sub_d.adjs, sub_d.ks, sub_d.ids_globais_vols, sub_d.ns
        nv=ns[0]
        ne=ns[1]
        nf=ns[2]
        ni=ns[3]

        adjs0=adjs[:,0]
        adjs1=adjs[:,1]

        II=self.get_submatrix(adjs0, adjs1, ks, (0, ni, 0, ni))
        IF=self.get_submatrix(adjs0, adjs1, ks, (0, ni, ni, ni+nf))

        if sub_d.local_couple>0:
            IE=self.get_submatrix(adjs0, adjs1, ks, (0, ni, ni+nf, ni+nf+ne))
            IV=self.get_submatrix(adjs0, adjs1, ks, (0, ni, ni+nf+ne, ni+nf+ne+nv))

        FF=self.get_submatrix(adjs0, adjs1, ks, (ni, ni+nf, ni, ni+nf))
        FE=self.get_submatrix(adjs0, adjs1, ks, (ni,ni+nf, ni+nf,ni+nf+ne))

        if sub_d.local_couple>0:
            FV=self.get_submatrix(adjs0, adjs1, ks, (ni, ni+nf, ni+nf+ne, ni+nf+ne+nv))

        EE=self.get_submatrix(adjs0, adjs1, ks, (ni+nf, ni+nf+ne, ni+nf, ni+nf+ne))
        EV=self.get_submatrix(adjs0, adjs1, ks, (ni+nf, ni+nf+ne, ni+nf+ne, ni+nf+ne+nv))

        Pv=scipy.sparse.identity(nv)
        t0=time.time()
        Pe=-linalg.spsolve(EE,EV*Pv)
        sub_d.A_b_t.append([EE.shape[0], nv,time.time()-t0])
        if sub_d.local_couple==0:
            t0=time.time()
            Pf=-linalg.spsolve(FF,FE*Pe)
            sub_d.A_b_t.append([FF.shape[0], nv,time.time()-t0])
            t0=time.time()
            Pi=-linalg.spsolve(II,IF*Pf)
            sub_d.A_b_t.append([II.shape[0], nv,time.time()-t0])
        else:
            t0=time.time()
            Pf=-linalg.spsolve(FF,FE*Pe+FV)
            sub_d.A_b_t.append([FF.shape[0], nv,time.time()-t0])
            t0=time.time()
            Pi=-linalg.spsolve(II,IF*Pf+IE*Pe+IV)
            sub_d.A_b_t.append([II.shape[0], nv,time.time()-t0])
        try:
            OP=vstack([Pi,Pf,Pe,Pv])
        except:
            import pdb; pdb.set_trace()

        lcd=scipy.sparse.find(OP)
        lines=ids_globais_vols[np.array(lcd[0])].astype(int)
        cols = sub_d.coarse_ids[lcd[1]].astype(int)
        data=np.array(lcd[2])
        return lines, cols, data

class Partitioner:
    def __init__(self,all_subds, nworker, regression_degree):
        if len(all_subds)>0:
            estimated_time_by_subd = self.get_estimated_time_by_subd(all_subds,regression_degree)
            partitioned_subds = self.balance_processes(all_subds, estimated_time_by_subd, nworker=nworker)
            self.partitioned_subds = partitioned_subds
        else:
            self.partitioned_subds=[]
        # A_b_t=np.zeros((1,3))
        # for subd in all_subds:
        #     A_b_t=np.vstack([A_b_t,np.array(subd.A_b_t)])
        # A_b_t=A_b_t[1:,:]
        # try:
        #     Abt=np.load("flying/A_b_t.npy")
        #     A_b_t=np.vstack([A_b_t,Abt])
        #     np.save("flying/A_b_t.npy",A_b_t)
        # except:
        #     np.save("flying/A_b_t.npy",A_b_t)


    def get_estimated_time_by_subd(self, all_subds, regression_degree = 2):
        n_A = [np.array(subd.ns) for subd in all_subds]

        n_b = [subd.ns[0] for subd in all_subds]
        try:
            n_A=np.array(n_A)[:,1:]
        except:
            import pdb; pdb.set_trace()
        n_b=np.array(n_b)
        if regression_degree==1:
            print("linear")
            cx, cy, intercept = np.load("flying/partitioning_coeffitients_cx_cy_intercept.npy")
            cx2, cxy, cy2 = 0, 0, 0
        else:
            print("quadratico")
            try:
                cx, cy, cx2, cxy, cy2, intercept= np.load("flying/partitioning_coeffitients_bx_cy_dx2_exy_fy2_intercept.npy")
            except:
                calibrate_partitioning_parameters()
                cx, cy, cx2, cxy, cy2, intercept= np.load("flying/partitioning_coeffitients_bx_cy_dx2_exy_fy2_intercept.npy")

        x=n_A
        y=np.array([n_b]).T
        estimated_time_by_subd=(cx*x+cy*y+cx2*x*x+cxy*x*y+cy2*y*y+intercept).sum(axis=1)
        return estimated_time_by_subd

    def balance_processes(self, all_subds, estimated_time_by_subd, nworker=1):
        if nworker>len(all_subds):
            print("more workers than subdomains, working with {} processes".format(len(all_subds)))
            nworker=len(all_subds)

        parts = np.zeros((nworker,len(all_subds)))
        u_vals=-np.sort(np.unique(-estimated_time_by_subd))
        for u in u_vals:
            posics=np.arange(len(estimated_time_by_subd))[estimated_time_by_subd==u]
            for p in posics:
                worker_id=np.arange(nworker)[parts.sum(axis=1)==parts.sum(axis=1).min()][0]
                parts[worker_id,p]=estimated_time_by_subd[p]
        if (parts!=0).sum(axis=0).min()!=1 or (parts>0).sum(axis=0).min()!=1:
            print("verificar particionamento")
            import pdb; pdb.set_trace()
        print(parts.sum(axis=1), (parts>0).sum(axis=1),len(all_subds),"aqui")

        partitioned_subds=[]
        for i in range(nworker):
            partitioned_subds.append(np.array(all_subds)[parts[i]>0])

        return partitioned_subds

class OP_AMS:
    def __init__(self, data_impress, elements_lv0, all_conjs_duais, local_couple=0, couple_bound=False):
        t0=time.time()

        print("Time to calibrate partitioning parameters: {} segundos".format(time.time()-t0))
        t0=time.time()
        all_subds = [DualDomain(data_impress, elements_lv0, all_conjs_duais[i], local_couple=local_couple, \
        couple_bound = couple_bound) for i in range(len(all_conjs_duais))]
        print("Time to partitionate subdomains: {} segundos".format(time.time()-t0))
        Nvols=len(elements_lv0['volumes'])
        Nverts = (data_impress['DUAL_1']==3).sum()
        regression_degree=2
        nworker=3

        partitioned_subds=Partitioner(all_subds, nworker, regression_degree).partitioned_subds

        lines, cols, data = self.get_OP_paralell(partitioned_subds)
        self.OP=csc_matrix((data,(lines,cols)),shape=(Nvols,Nverts))


        # # To test bugs on serial, use this###############################
        # lines, cols, data = self.get_OP(all_subds, paralell=False)
        # self.OP=csc_matrix((data,(lines,cols)),shape=(Nvols,Nverts))
        #######################################



    def get_OP(self,partitioned_subd, paralell=True):
        if paralell:
            print("process {} started".format(partitioned_subd[-1].id))
        t0=time.time()
        # Processes imputs
        ################################
        lcd=np.zeros((3,1))
        for dual_d in partitioned_subd:
            lcd=np.hstack([lcd,OP_local(dual_d).OP])
        ###################################

        # Send results to master process
        ###################################################
        if paralell:
            print("process {} finished after {}".format(partitioned_subd[-1].id, time.time()-t0))
            master=dual_d.master
            master.send(lcd)
            #############################################
        else:
            return lcd

    def get_OP_paralell(self, partitioned_subds):
        nworker = len(partitioned_subds)
        print("calculating prolongation operator with {} processes".format(nworker))

        # Setup communication structure
        #########################################
        master2worker = [mp.Pipe() for p in range(nworker)]
        m2w, w2m = list(zip(*master2worker))
        for i in range(len(partitioned_subds)):
            partitioned_subds[i][-1].master = w2m[i]
            partitioned_subds[i][-1].id = i
        ########################################

        # Creates & start processes
        #########################################
        procs = [mp.Process(target=self.get_OP, args=[s]) for s in partitioned_subds]
        for p in procs:
            p.start()
        #########################################

        # Get processed output & kill subprocesses
        #################################
        l=[]
        c=[]
        d=[]
        for m in m2w:
            msg=m.recv()
            l.append(msg[0])
            c.append(msg[1])
            d.append(msg[2])

        l=np.concatenate(l)
        c=np.concatenate(c)
        d=np.concatenate(d)
        for p in procs:
            p.join()
        ###############################################
        lines=l.astype(int)
        cols=c.astype(int)
        data=d
        return lines, cols, data
