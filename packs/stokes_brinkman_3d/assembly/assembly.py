import numpy as np

class global_assembly:
    def __init__(self, sb, M, dx, dy, dz):
        self.nv=len(M.volumes.all)
        self.nfi=len(M.faces.internal)
        self.mi_l=0.1
        self.M_cont=self.get_continuity_matrix(M, dx, dy, dz)
        self.M_momentum_x=self.get_momentum_matrix(sb.fx, M,dx,dy, dz,0)
        self.M_momentum_y=self.get_momentum_matrix(sb.fy, M,dx,dy, dz,1)
        self.M_momentum_z=self.get_momentum_matrix(sb.fz, M,dx,dy, dz,2)
        self.M=np.vstack([self.M_cont,self.M_momentum_x,self.M_momentum_y,self.M_momentum_z])

    def get_continuity_matrix(self, M, dx, dy, dz):
        M_cont=np.zeros((self.nv,self.nv+self.nfi))
        internal_faces=M.faces.internal
        volumes=M.volumes.all
        faces=M.volumes.bridge_adjacencies(volumes,2,2)
        faces_center_x=M.faces.center(faces.flatten())[:,0].reshape(faces.shape)
        faces_center_y=M.faces.center(faces.flatten())[:,1].reshape(faces.shape)
        faces_center_z=M.faces.center(faces.flatten())[:,2].reshape(faces.shape)

        for i in range(len(faces)):
            facs=faces[i]
            xcents=faces_center_x[i]
            ycents=faces_center_y[i]
            zcents=faces_center_z[i]

            xfe=facs[xcents==xcents.min()]
            xfd=facs[xcents==xcents.max()]
            if xfe[0] in internal_faces:
                ide=M.id_fint[xfe[0]]
                M_cont[volumes[i],self.nv+ide[0]]=-1

            if xfd[0] in internal_faces:
                idd=M.id_fint[xfd[0]]
                M_cont[volumes[i],self.nv+idd[0]]=1

            yfi=facs[ycents==ycents.min()]
            yfs=facs[ycents==ycents.max()]
            if yfi[0] in internal_faces:
                idi=M.id_fint[yfi[0]]
                M_cont[volumes[i],self.nv+idi[0]]=-1

            if yfs[0] in internal_faces:
                ids=M.id_fint[yfs[0]]
                M_cont[volumes[i],self.nv+ids[0]]=1

            zfi=facs[zcents==zcents.min()]
            zfs=facs[zcents==zcents.max()]
            if zfi[0] in internal_faces:
                idi=M.id_fint[zfi[0]]
                M_cont[volumes[i],self.nv+idi[0]]=-1

            if zfs[0] in internal_faces:
                ids=M.id_fint[zfs[0]]
                M_cont[volumes[i],self.nv+ids[0]]=1
        return(M_cont)

    def get_momentum_matrix(self, fd ,M,dx,dy, dz,col):
        T=1
        mi=1
        k=1

        adjs_d=M.faces.bridge_adjacencies(fd,2,3)
        adjs_d0=adjs_d[:,0]
        adjs_d1=adjs_d[:,1]

        M_d=np.zeros((len(fd),self.nv+self.nfi))
        c0=M.volumes.center(adjs_d0)[:,col]
        c1=M.volumes.center(adjs_d1)[:,col]

        higher_coord=np.zeros(len(fd),dtype=int)
        higher_coord[c0>c1]=adjs_d0[c0>c1]
        higher_coord[c0<c1]=adjs_d1[c0<c1]

        lower_coord=np.zeros(len(fd),dtype=int)
        lower_coord[c0<c1]=adjs_d0[c0<c1]
        lower_coord[c0>c1]=adjs_d1[c0>c1]

        fd_l=M.id_fint[fd].T[0]
        fd_l-=fd_l.min()


        M_d[fd_l,higher_coord]=1
        M_d[fd_l,lower_coord]=-1

        fac_viz_ares = M.faces.bridge_adjacencies(fd,1,2)
        laterais=[]
        for fac in fac_viz_ares:
            laterais.append(np.intersect1d(fac,fd))

        laterais_l=[]
        for l in laterais:
            laterais_l.append(M.id_fint[l].T[0])
        faces_viz_adjs_d=M.volumes.bridge_adjacencies(adjs_d.flatten(),2,2).reshape(len(fd),12) #12 é o número de faces 6faces/vol vezes 2 vols por adj

        sup_inf=[]
        i=0
        for fac in faces_viz_adjs_d:
            sup_inf.append(np.setdiff1d(np.intersect1d(fac,fd),fd[i]))
            i+=1
        sup_inf_l=[]
        for si in sup_inf:
            sup_inf_l.append(M.id_fint[si].T[0])
        fd_l_orig=M.id_fint[fd].T[0]
        for i in range(len(fd)):
            M_d[fd_l[i],self.nv+fd_l_orig[i]]=1-self.mi_l/(dy*dy)
            M_d[fd_l[i],self.nv+sup_inf_l[i]]=self.mi_l/(dy*dy)

            if len(laterais[i])==1:
                M_d[fd_l[i],self.nv+fd_l_orig[i]]-=self.mi_l/(dx*dx)
                M_d[fd_l[i],self.nv+laterais_l[i]]+=1*self.mi_l/(dx*dx)
            else:
                M_d[fd_l[i],self.nv+fd_l_orig[i]]-=self.mi_l/(dx*dx)
                M_d[fd_l[i],self.nv+fd_l_orig[i]]+=self.mi_l/(dx*dx)
        return M_d
