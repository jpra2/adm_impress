import scipy
import numpy as np
import multiprocessing as mp
from scipy.sparse import csc_matrix,csr_matrix, linalg, vstack, find
import time

class prolongation_operator():
    def __init__(self, global_ids_internal_faces, coupled_edges_local_ids=None):
        self.OP_part = self.get_OP_part(global_ids_internal_faces,coupled_edges_local_ids=None)

    def get_submatrix(self,id0, id1, ks, xi, xs, yi, ys,sum_diag=None,just_some_lines=None):
        if just_some_lines:
            just_some_lines=np.array(just_some_lines)[0]
            s0=(id0==just_some_lines[0]) & (id1>=yi) & (id1<ys)

            for j in just_some_lines[1:]:
                s0+=(id0==j) & (id1>=yi) & (id1<ys)
            p0=s0>0
            l=id0[p0]-xi
            c=id1[p0]-yi
            d=ks[p0]

        else:
            inds =(id0>=xi) & (id0<xs) & (id1>=yi) & (id1<ys)
            l=id0[inds]-xi
            c=id1[inds]-yi
            d=ks[inds]

        if xi==yi:
            c=np.concatenate([c,l])
            l=np.concatenate([l,l])
            d=np.concatenate([d,-d])

        if sum_diag:
            #Add to diagonal the terms in other submatrices
            for diagon in sum_diag:
                dia=-np.array(diagon.sum(axis=1).transpose())[0]
                n_lines=len(dia)
                l=np.concatenate([l,range(n_lines)])
                c=np.concatenate([c,range(n_lines)])
                d=np.concatenate([d,dia])
        try:
            submatrix=csc_matrix((d,(l,c)),shape=(xs-xi,ys-yi))
        except:
            pass
        return(submatrix)

    def map_local_properties(self,global_ids_internal_faces):
        local_adjs_with_global_id=all_adjacent_volumes[global_ids_internal_faces]
        local_transmissibility=kst[global_ids_internal_faces]
        ids_globais_vols=np.unique(local_adjs_with_global_id.flatten())

        internal_global_ids=ids_globais_vols[ids_globais_vols<ni]
        faces_global_ids=ids_globais_vols[(ids_globais_vols>=ni) & (ids_globais_vols<ni+nf)]
        edge_global_ids=ids_globais_vols[(ids_globais_vols>=ni+nf) & (ids_globais_vols<ni+nf+na)]
        vertex_global_ids=ids_globais_vols[(ids_globais_vols>=ni+nf+na) & (ids_globais_vols<ni+nf+na+nv)]

        ids_globais_vols=np.concatenate([internal_global_ids, faces_global_ids, edge_global_ids, vertex_global_ids])

        n_intern=len(internal_global_ids)
        n_faces=len(faces_global_ids)
        n_arests=len(edge_global_ids)
        n_verts=len(vertex_global_ids)
        n_vols=len(ids_globais_vols)

        ids_locais_vols=range(n_vols)

        all_ids_globais=np.arange(n_volumes)
        all_ids_globais[ids_globais_vols]=ids_locais_vols

        ids0=all_ids_globais[local_adjs_with_global_id[:,0]]  #Esses são os ids locais
        ids1=all_ids_globais[local_adjs_with_global_id[:,1]]

        i0=ids0.copy()
        ids0=np.concatenate([ids0,ids1])
        ids1=np.concatenate([ids1,i0])

        ks=local_transmissibility
        ks=np.concatenate([ks,ks])
        return(ids_globais_vols, vertex_global_ids,ids0, ids1, ks, n_intern, n_faces, n_arests, n_verts)

    def lu_inv3(self,M,lines):
        lines=np.array(lines)
        L=len(lines)
        s=1000
        n=int(L/s)
        r=int(L-int(L/s)*s)
        tinv=time.time()
        LU=linalg.splu(M)
        if L<s:
            l=lines
            c=range(len(l))
            d=np.repeat(1,L)
            B=csr_matrix((d,(l,c)),shape=(M.shape[0],L))
            B=B.toarray()
            inversa=csc_matrix(LU.solve(B,'T')).transpose()
        else:
            c=range(s)
            d=np.repeat(1,s)
            for i in range(n):
                l=lines[s*i:s*(i+1)]
                B=csc_matrix((d,(l,c)),shape=(M.shape[0],s))
                B=B.toarray()
                if i==0:
                    inversa=csc_matrix(LU.solve(B,'T')).transpose()
                else:
                    inversa=csc_matrix(vstack([inversa,csc_matrix(LU.solve(B,'T')).transpose()]))

            if r>0:
                l=lines[s*n:L]
                c=range(r)
                d=np.repeat(1,r)
                B=csc_matrix((d,(l,c)),shape=(M.shape[0],r))
                B=B.toarray()
                inversa=csc_matrix(vstack([inversa,csc_matrix(LU.solve(B,'T')).transpose()]))
        f=find(inversa)
        ll=f[0]
        c=f[1]
        d=f[2]
        pos_to_line=dict(zip(range(len(lines)),lines))
        lg=[pos_to_line[l] for l in ll]
        inversa=csc_matrix((d,(lg,c)),shape=(M.shape[0],M.shape[0]))
        #print(time.time()-tinv,L,"tempo de inversão")
        return inversa

    def get_OP_part(self,global_ids_internal_faces, pid, coupled_edges_local_ids=None):
        t7=time.time()
        ids_globais_vols, vertex_global_ids, ids0, ids1, ks, n_intern, n_facs, n_arests, n_verts= self.map_local_properties(global_ids_internal_faces)

        t8=time.time()
        Aif=self.get_submatrix(ids0, ids1, ks, 0,n_intern, n_intern, n_intern+n_facs)
        Aii=self.get_submatrix(ids0, ids1, ks, 0,n_intern, 0, n_intern, [Aif])

        Afe=self.get_submatrix(ids0, ids1, ks, n_intern,n_intern+n_facs, n_intern+n_facs, n_intern+n_facs+n_arests)
        Aff=self.get_submatrix(ids0, ids1, ks, n_intern,n_intern+n_facs, n_intern, n_intern+n_facs,[Afe])
        Aev=self.get_submatrix(ids0, ids1, ks, n_intern+n_facs,n_intern+n_facs+n_arests, n_intern+n_facs+n_arests, n_intern+n_facs+n_arests+n_verts)
        try:
            couple=len(coupled_edges_local_ids)
        except:
            couple=False

        if couple:
            ef=self.get_submatrix(ids0, ids1, ks, n_intern+n_facs,n_intern+n_facs+n_arests, n_intern, n_intern+n_facs,just_some_lines=[coupled_edges_local_ids])
            Aee=self.get_submatrix(ids0, ids1, ks, n_intern+n_facs,n_intern+n_facs+n_arests, n_intern+n_facs, n_intern+n_facs+n_arests,[Aev, ef])

            cnl=np.array(ef.sum(axis=0))[0]!=0
            cols=np.arange(len(cnl))[cnl]
            Aee-=ef*self.lu_inv3(Aff,cols)*Afe #Lu_inv3 inverte apenas as linhas solicitadas

        else:
            Aee=self.get_submatrix(ids0, ids1, ks, n_intern+n_facs,n_intern+n_facs+n_arests, n_intern+n_facs, n_intern+n_facs+n_arests,[Aev])

        Pv=scipy.sparse.identity(n_verts)
        Pe=-linalg.spsolve(Aee,Aev)
        Pf=-linalg.spsolve(Aff,Afe*Pe)
        Pi=-linalg.spsolve(Aii,Aif*Pf)
        OP_part=vstack([Pi,Pf,Pe,Pv])
        lcd=find(OP_part)

        lines=ids_globais_vols[np.array(lcd[0])]
        cols=vertex_global_ids[np.array(lcd[1])]-ni-nf-na
        data=np.array(lcd[2])
        return(lines,cols,data)

class SubDomain():
    id = 0
    def __init__(self,dual_blocks,id,coupled_edges_local_ids, master):
        self.dual_blocks = dual_blocks
        self.coupled_edges_local_ids=coupled_edges_local_ids
        self.id = id
        SubDomain.id += 1
        self.master = master

class acumulate_OP(prolongation_operator):
    def __init__(self,sd):
        ttt=time.time()
        dual_blocks=sd.dual_blocks
        coupled_edges=sd.coupled_edges_local_ids
        id = sd.id
        print("Process {} started".format(id))
        master = sd.master
        lines=[]
        cols=[]
        data=[]
        for d in dual_blocks:
            l, c, d = self.get_OP_part(d,id,coupled_edges)
            lines.append(l)
            cols.append(c)
            data.append(d)
        master.send(np.vstack([lines,cols,data]))
        print("Process {} finished after {} seconds".format(id,time.time()-ttt))

class OP:
    def __init__(self, OP_input):
        nworker=OP_input.nworker
        first=OP_input.first
        ids=OP_input.global_ids_internal_faces
        coupled=OP_input.coupled_edges_local_ids
        coupling=OP_input.realize_coupling
        self.OP = self.get_OP(nworker,first,ids,coupled, coupling)

    def get_OP(self, nworker, first, ids, coupled, coupling):
        if first:
            first=False
            print("Main program started")
            if coupling:
                global_ids_internal_faces = ids
                coupled_edges_local_ids   = np.repeat(None,len(global_ids_internal_faces))
            else:
                global_ids_internal_faces = ids
                coupled_edges_local_ids   = coupled

            global all_adjacent_volumes, kst, n_entities, ni, nf, na, nv, n_volumes

            all_adjacent_volumes      = np.load("flying/all_adjacent_volumes.npy")  #IDs dos volumes adjascentes a cada uma das faces internas
            kst                       = np.load("flying/kst.npy")                   # Permeabilidade equivalente em cada uma das faces internas ao reservatório
            n_entities                = np.load("flying/n_entities.npy")            # Número de entidades de cada tipo [Internal, Face, Edge, Vertex]

            ni = n_entities[0]
            nf = n_entities[1]
            na = n_entities[2]
            nv = n_entities[3]
            n_volumes = ni+nf+na+nv
            # nworker = 6
            t1=time.time()

            n_duais=len(global_ids_internal_faces)
            if nworker>n_duais:
                nworker=n_duais
                print("Mais processos que volumes duais, usando agora {} processos".format(n_duais))

            conjs_duais=[]
            for n in range(nworker):
                if n==nworker-1:
                    conjs_duais.append(global_ids_internal_faces[n*int(n_duais/nworker):])
                else:
                    conjs_duais.append(global_ids_internal_faces[n*int(n_duais/nworker):(n+1)*int(n_duais/nworker)])

            if len(np.concatenate(conjs_duais))!=n_duais:
                print("Verificar as regras de distribuição das tarefas!")
                import pdb; pdb.set_trace()

            master2worker = [mp.Pipe() for p in range(nworker)]
            m2w, w2m = list(zip(*master2worker))
            subDomains = [SubDomain(conjs_duais[id],id,coupled_edges_local_ids[id], w2m[id]) for id in range(nworker)]
            procs = [mp.Process(target=acumulate_OP, args=[s]) for s in subDomains]

            for p in procs:
                p.start()

            l=[]
            c=[]
            d=[]
            for m in m2w:
                msg=m.recv()
                ''' Trata exceção devido ao tipo de dado: apenas uma dual por processo
                gera numpy array, mais que uma dual gera lista de numpy arrays'''
                try:
                    l.append(np.concatenate(msg[0]))
                    c.append(np.concatenate(msg[1]))
                    d.append(np.concatenate(msg[2]))
                except:
                    l.append(np.concatenate([msg[0]]))
                    c.append(np.concatenate([msg[1]]))
                    d.append(np.concatenate([msg[2]]))

            l=np.concatenate(l)
            c=np.concatenate(c)
            d=np.concatenate(d)

            for p in procs:
                p.join()

            t2=time.time()
            print("Total time with {} processes, {} seconds".format(nworker,t2-t1))


            OP=csc_matrix((d,(l,c)),shape=(n_volumes,nv))

            slin=OP.sum(axis=1)
            ops_rep=slin>0.9
            ops_nrep=slin<=0.9
            if ops_nrep.sum()>0:
                print("{} linhas foram calculadas e {} não foram".format(ops_rep.sum(),ops_nrep.sum()))

            den=np.array((np.round(slin)+ops_nrep)).transpose()[0]
            diag=csc_matrix((1/den,(range(n_volumes),range(n_volumes))),shape=(n_volumes,n_volumes))
            OP=diag*OP
            soma_cols=np.array(OP.sum(axis=1)[OP.sum(axis=1)>0])[0]
            # import pdb; pdb.set_trace()
            if soma_cols.max()<soma_cols.min()+0.0000001 and (slin>0.99).sum()==round(OP.sum()):
                print("Teste do Operador indica que seu cálculo está correto")
            else:
                print("verificar o operador")

            print("Main program ended")
            return OP
