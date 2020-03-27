from pymoab import rng, types
import scipy
import numpy as np
import multiprocessing as mp
from scipy.sparse import csc_matrix,csr_matrix, linalg, vstack, find
import time

class DualDomain:
    def __init__(self, data_impress, elements_lv0, volumes, local_couple=0, couple_bound=True):
        self.local_couple=local_couple
        self.couple_bound=couple_bound

        self.adjs, self.ks, self.ids_globais_vols, self.ns =  self.get_local_informations(data_impress, elements_lv0, volumes, local_couple=0, couple_bound=True)
        self.Nvols=len(elements_lv0['volumes'])
        self.Nvert = self.Nvert
        self.coarse_ids = data_impress['GID_1'][self.vertices]

    def get_local_informations(self, data_impress, elements_lv0, volumes, local_couple=0,           couple_bound=True):
        # viz=M.mtu.get_bridge_adjacencies(volumes,2,3)
        viz=np.unique(np.concatenate(elements_lv0['volumes_face_volumes'][volumes]))
        nv=len(volumes)
        # M.mb.tag_set_data(M.local_id_dual_tag,viz, np.repeat(nv+5,len(viz)))
        # faces_entities=M.mtu.get_bridge_adjacencies(volumes,2,2)
        faces_entities = np.unique(np.concatenate(elements_lv0['volumes_face_faces'][volumes]))
        int_facs=np.setdiff1d(faces_entities, elements_lv0['boundary_faces'])

        so_viz=np.setdiff1d(viz,volumes)
        if len(so_viz) > 0:
            so_viz_faces=np.unique(np.concatenate(elements_lv0['volumes_face_faces'][so_viz]))
            int_facs=np.setdiff1d(int_facs,so_viz_faces)

        if local_couple>0:
            dual_flags=M.mb.tag_get_data(M.D1_tag, np.array(volumes),flat=True)
            if couple_bound:
                reduce_flag=np.repeat(True,len(volumes))
                volumes_red=np.array(volumes)
                volumes_red=np.array(volumes)[reduce_flag]
                dual_flags_red=dual_flags[reduce_flag]
                dual_flags_red[dual_flags_red==2]-=1
                if local_couple==2:
                    dual_flags_red[dual_flags_red==1]=2
                    dual_flags_red[dual_flags_red==0]=2

                # dual_flags[dual_flags!=3]=2
                # dual_flags_red=dual_flags

            else:
                M.mb.tag_set_data(M.local_id_dual_tag,np.array(volumes), np.arange(len(volumes)))
                adjs=[M.mb.get_adjacencies(f,3) for f in int_facs]
                adjs=np.array(adjs)
                adjs=M.mb.tag_get_data(M.local_id_dual_tag,np.concatenate(adjs),flat=True).reshape(len(adjs),2)

                v0=adjs[:,0]<len(volumes)
                v1=adjs[:,1]<len(volumes)
                vv=v0&v1

                adjs=adjs[vv]
                adjs0=np.concatenate([adjs[:,0],adjs[:,1]])
                adjs1=np.concatenate([adjs[:,1],adjs[:,0]])

                T=csc_matrix((np.ones(len(adjs0)),(adjs0,adjs1)),shape=(len(volumes),len(volumes)))
                reduce_flag=np.array(T.sum(axis=1)==T.sum(axis=1).max()).T[0] & (dual_flags!=3)

                volumes_red=np.array(volumes)[reduce_flag]
                dual_flags_red=dual_flags[reduce_flag]
                dual_flags_red[dual_flags_red==2]-=1
                if local_couple==2:
                    dual_flags_red[dual_flags_red==1]=2
                    dual_flags_red[dual_flags_red==0]=2


            M.mb.tag_set_data(M.D1_tag,volumes_red,dual_flags_red)

        # ms=M.mb.create_meshset()
        # M.mb.add_entities(ms,volumes)
        # M.mb.write_file("results/test_couple.vtk",[ms])

        dual_id_volumes = data_impress['DUAL_1'][volumes]

        # vertices=M.mb.get_entities_by_type_and_tag(ms, types.MBHEX, np.array([M.D1_tag]), np.array([3]))
        # edges=M.mb.get_entities_by_type_and_tag(ms, types.MBHEX, np.array([M.D1_tag]), np.array([2]))
        # faces=M.mb.get_entities_by_type_and_tag(ms, types.MBHEX, np.array([M.D1_tag]), np.array([1]))
        # internals=M.mb.get_entities_by_type_and_tag(ms, types.MBHEX, np.array([M.D1_tag]), np.array([0]))

        vertices = volumes[dual_id_volumes==3]
        edges = volumes[dual_id_volumes==2]
        faces = volumes[dual_id_volumes==1]
        internals = volumes[dual_id_volumes==0]

        nv=len(vertices)
        ne=len(edges)
        nf=len(faces)
        ni=len(internals)
        ns=[nv,ne,nf,ni]

        self.Nvert = nv
        self.vertices = vertices

        # M.mb.tag_set_data(M.local_id_dual_tag,internals, np.arange(ni))
        # M.mb.tag_set_data(M.local_id_dual_tag,faces, np.arange(ni,ni+nf))
        # M.mb.tag_set_data(M.local_id_dual_tag,edges, np.arange(ni+nf,ni+nf+ne))
        # M.mb.tag_set_data(M.local_id_dual_tag,vertices, np.arange(ni+nf+ne,ni+nf+ne+nv))

        # adjs=[M.mb.get_adjacencies(f,3) for f in int_facs]
        # adjs=np.array(adjs)
        #
        # adjs=M.mb.tag_get_data(M.local_id_dual_tag,np.concatenate(adjs),flat=True).reshape(len(adjs),2)

        adjs = elements_lv0['faces_face_volumes'][int_facs]
        adjs = np.concatenate(adjs).reshape((adjs.shape[0], 2))

        map_l=np.zeros(adjs.max()+1)
        map_l[internals]=np.arange(ni)
        map_l[faces]=np.arange(ni,ni+nf)
        map_l[edges]=np.arange(ni+nf,ni+nf+ne)
        map_l[vertices]=np.arange(ni+nf+ne, ni+nf+ne+nv)

        adjs_l0=map_l[adjs[:,0]]
        adjs_l1=map_l[adjs[:,1]]
        # import pdb; pdb.set_trace()
        # adjs0=adjs[:,0]
        # adjs1=adjs[:,1]
        #
        # v0=adjs0<len(volumes)
        # v1=adjs1<len(volumes)
        # vv=v0&v1
        adjs=np.array([adjs_l0, adjs_l1]).T
        # ks=M.mb.tag_get_data(M.k_eq_tag,np.uint64(int_facs)[vv],flat=True)
        ks=data_impress['transmissibility'][int_facs]
        # ids_globais_vols=M.mb.tag_get_data(M.ID_reordenado_tag,np.concatenate([np.uint64(internals),np.uint64(faces), np.uint64(edges),vertices]),flat=True)
        ids_globais_vols=np.concatenate([np.uint64(internals),np.uint64(faces), np.uint64(edges),vertices])
        return adjs, ks, ids_globais_vols, ns

class OP_local:
    def __init__(self, sub_d):
        self.Nvols=sub_d.Nvols
        self.Nverts=sub_d.Nvert
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
        # vertex_global_ids=ids_globais_vols[sub_d.vertices]
        # vertex_global_ids=sub_d.vertices
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
        Pe=-linalg.spsolve(EE,EV*Pv)
        if sub_d.local_couple==0:
            Pf=-linalg.spsolve(FF,FE*Pe)
            Pi=-linalg.spsolve(II,IF*Pf)
        else:
            Pf=-linalg.spsolve(FF,FE*Pe+FV)
            Pi=-linalg.spsolve(II,IF*Pf+IE*Pe+IV)

        OP=vstack([Pi,Pf,Pe,Pv])

        lcd=scipy.sparse.find(OP)
        lines=ids_globais_vols[np.array(lcd[0])]
        # cols=vertex_global_ids[np.array(lcd[1])]-ni-nf-ne
        cols = sub_d.coarse_ids[lcd[1]]
        data=np.array(lcd[2])
        OP=scipy.sparse.csc_matrix((data,(lines,cols)), shape=(self.Nvols,self.Nverts))
        print(OP.sum(axis=1).max(),OP.sum(axis=1).min())

        return OP

class OP_AMS:
    def __init__(self, data_impress, elements_lv0, volumes, local_couple=0, couple_bound=True):
        dual_d=DualDomain(data_impress, elements_lv0, volumes, local_couple=0, couple_bound=True)
        self.OP=OP_local(dual_d).OP
