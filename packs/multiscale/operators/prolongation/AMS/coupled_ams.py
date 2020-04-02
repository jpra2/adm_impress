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
        self.adjs, self.ks, self.ids_globais_vols, self.ns =  self.get_local_informations(data_impress, elements_lv0, volumes, local_couple=local_couple, couple_bound=couple_bound)
        self.Nvols=len(elements_lv0['volumes'])
        self.Nvert = (data_impress['DUAL_1']==3).sum()
        self.coarse_ids = data_impress['GID_1'][self.vertices]
        self.A_b_t=[]

    def get_local_informations(self, data_impress, elements_lv0, volumes, local_couple=0, couple_bound=True):
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

        dual_id_volumes = data_impress['DUAL_1'][volumes]
        if local_couple>0:
            # dual_flags=M.mb.tag_get_data(M.D1_tag, np.array(volumes),flat=True)
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

                reduce_flag = np.setdiff1d(volumes, np.concatenate(elements_lv0['volumes_face_volumes'][so_viz]))

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
    def __init__(self, data_impress, elements_lv0, all_conjs_duais, local_couple=0, couple_bound=True):
        all_subds = [DualDomain(data_impress, elements_lv0, all_conjs_duais[i], local_couple=local_couple, \
        couple_bound = couple_bound) for i in range(len(all_conjs_duais))]
        partitioned_subds = self.partitionate_subds(all_subds, nworker=1)
        self.OP=self.get_OP(partitioned_subds)
        self.coefs=[]
        A_b_t=np.zeros((1,3))
        for subd in all_subds:
            A_b_t=np.vstack([A_b_t,np.array(subd.A_b_t)])
        A_b_t=A_b_t[1:,:]
        try:
            Abt=np.load("flying/A_b_t.npy")
            A_b_t=np.vstack([A_b_t,Abt])
            np.save("flying/A_b_t.npy",A_b_t)
        except:
            np.save("flying/A_b_t.npy",A_b_t)
        xs=A_b_t[:,0]
        ys=A_b_t[:,1]
        zs=A_b_t[:,2]
        print(len(A_b_t))
        # x1=np.random.rand(200,2)
        # xs=x1[:,0]
        # ys=x1[:,1]
        # zs=1*xs+2*ys+3*xs*xs+4*xs*ys+5*ys*ys+6
        import matplotlib.pyplot as plt

        fig = plt.figure()
        name="A_b_t"
        # self.imprima3d(xs, ys, zs, name)
        self.imprima3d2(xs, ys,zs, fig)
        self.reg3(xs, ys, zs, fig)
        plt.close('all')


    def partitionate_subds(self, all_subds, nworker=1):
        n_A = [np.array(subd.ns) for subd in all_subds]
        n_b = [subd.ns[0] for subd in all_subds]

        partitioned_subds=all_subds
        return partitioned_subds

    def get_OP(self,partitioned_subds):
        for dual_d in partitioned_subds:
            try:
                OP+=OP_local(dual_d).OP
            except:
                OP=OP_local(dual_d).OP
        return OP

    def imprima3d(self,xs, ys, zs, name):
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.ticker import LinearLocator, FormatStrFormatter
        import numpy as np



        ax = fig.gca(projection='3d')


        nn=int(xs.shape[0])

        surf = ax.plot_surface(xs.reshape(nn,1), ys.reshape(nn,1), zs.reshape(nn,1), cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)

        ax.set_zlim(zs.min()-1, zs.max()+1)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))


        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.savefig("results/"+name+".png")

    def imprima3d2(self, xs, ys, zs, fig):
        # This import registers the 3D projection, but is otherwise unused.
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

        import matplotlib.pyplot as plt
        import numpy as np



        ax = fig.add_subplot(111, projection='3d')

        n = 100
        m="o"

        ax.scatter(xs, ys, zs, marker=m)

        ax.set_xlabel('A')
        ax.set_ylabel('b')
        ax.set_zlabel('time_to_solve (s)')
        plt.savefig("results/teste_scatter3d.png")

    def reg3(self, xs, ys, zs, fig):
        from sklearn.preprocessing import PolynomialFeatures
        from sklearn.linear_model import LinearRegression
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
        import matplotlib.pyplot as plt
        import sympy
        # x1=np.random.rand(10,2)
        # z1=1*x1[:,0]+2*x1[:,1]+3*x1[:,0]**2+4*x1[:,0]*x1[:,1]+5*x1[:,1]**2+6
        xx=np.array([xs])
        dx2=(abs(xx.T-xx)**2).sum(axis=1)

        yy=np.array([ys])
        dy2=(abs(yy.T-yy)**2).sum(axis=1)
        d2=dx2+dy2
        d=d2**0.5
        pesos=d/d.sum()


        degree=1
        poly = PolynomialFeatures(degree=degree)
        X_t = poly.fit_transform(np.array([xs,ys]).T,pesos)

        # X_t = poly.fit_transform(x1)
        clf = LinearRegression()
        clf.fit(X_t, zs)
        # clf.fit(X_t, z1)
        print(np.array(clf.coef_)[1:])

        print(clf.intercept_)

        self.coefs = np.concatenate([np.array(clf.coef_)[1:],np.array([clf.intercept_])])
        x, y = sympy.symbols('x y')

        coefs=np.array(clf.coef_)[1:]
        intercept=clf.intercept_
        if degree==2:
            syms=np.array([x, y, x*x, x*y, y*y])
        if degree==1:
            syms=np.array([x, y])
        func=(coefs*syms).sum()+intercept

        ax = fig.add_subplot(111, projection='3d')
        m="o"
        ax.scatter(xs, ys, zs, marker=m)
        ax.set_xlabel('A')
        ax.set_ylabel('b')
        ax.set_zlabel('time_to_solve (s)')
        plt.show()
        sympy.plotting.plot3d(func, (x, xs.min(), xs.max()), (y, ys.min(), ys.max()))
        plt.savefig("results/sym_surface2.png")
