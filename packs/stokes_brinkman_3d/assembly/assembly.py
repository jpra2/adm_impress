import numpy as np
from scipy import sparse

class global_assembly:
    def __init__(self, sb, M):
        self.M=self.get_lhs(M, sb)
        self.rhs=self.get_rhs(sb)

    def get_lhs(self, M, sb):
        self.mi_l=0.1
        lcd_cont=self.get_continuity_matrix(M, sb)
        # LHS=M_cont
        lcd=lcd_cont

        if len(sb.fx)>0:
            lcdx=self.get_momentum_matrix(sb.fx, M, sb,0)
            # LHS =np.vstack([LHS,M_momentum_x])
            lcd[0]=np.concatenate([lcd[0],lcdx[0]+sb.nv])
            lcd[1]=np.concatenate([lcd[1],lcdx[1]])
            lcd[2]=np.concatenate([lcd[2],lcdx[2]])
        if len(sb.fy)>0:
            lcdy=self.get_momentum_matrix(sb.fy, M, sb,1)
            # LHS =np.vstack([LHS,M_momentum_y])
            lcd[0]=np.concatenate([lcd[0],lcdy[0]+sb.nv+len(sb.fx)])
            lcd[1]=np.concatenate([lcd[1],lcdy[1]])
            lcd[2]=np.concatenate([lcd[2],lcdy[2]])
        if len(sb.fz)>0:
            lcdz=self.get_momentum_matrix(sb.fz, M, sb,2)
            # LHS =np.vstack([LHS,M_momentum_z])
            lcd[0]=np.concatenate([lcd[0],lcdz[0]+sb.nv+len(sb.fx)+len(sb.fy)])
            lcd[1]=np.concatenate([lcd[1],lcdz[1]])
            lcd[2]=np.concatenate([lcd[2],lcdz[2]])
        # LHS =np.vstack([M_cont,M_momentum_x,M_momentum_y,M_momentum_z])
        # import pdb; pdb.set_trace()
        lhs=sparse.csc_matrix((lcd[2],(lcd[0],lcd[1])),shape=(sb.nv+sb.nfi,sb.nv+sb.nfi))
        # import pdb; pdb.set_trace()
        return lhs

    def get_continuity_matrix(self, M, sb):
        ds=np.array([sb.dx, sb.dy, sb.dz])
        M_cont=np.zeros((sb.nv,sb.nv+sb.nfi))
        internal_faces=M.faces.internal
        volumes=M.volumes.all
        faces=M.volumes.bridge_adjacencies(volumes,2,2)
        faces_center_x=M.faces.center(faces.flatten())[:,0].reshape(faces.shape)
        faces_center_y=M.faces.center(faces.flatten())[:,1].reshape(faces.shape)
        faces_center_z=M.faces.center(faces.flatten())[:,2].reshape(faces.shape)
        lines=[]
        cols=[]
        data=[]
        for i in range(len(faces)):
            facs=faces[i]
            xcents=faces_center_x[i]
            ycents=faces_center_y[i]
            zcents=faces_center_z[i]
            cents=([xcents, ycents, zcents])
            for cent in cents:
                dfe=facs[cent==cent.min()]
                dfd=facs[cent==cent.max()]

                if dfe[0] in internal_faces:
                    ide=M.id_fint[dfe[0]]
                    M_cont[volumes[i],sb.nv+ide[0]]=-1
                    #####################
                    lines.append(volumes[i])
                    cols.append(sb.nv+ide[0])
                    data.append(-1)
                    ################

                if dfd[0] in internal_faces:
                    idd=M.id_fint[dfd[0]]
                    M_cont[volumes[i],sb.nv+idd[0]]=1
                    ###############
                    lines.append(volumes[i])
                    cols.append(sb.nv+idd[0])
                    data.append(1)
                    ##################

        lcd=[np.array(lines),np.concatenate(cols),np.array(data)]

        return lcd

    def get_momentum_matrix(self, fd ,M, sb,col):
        dx, dy, dz= sb.dx, sb.dy, sb.dz
        k_harms=M.k_harm[M.faces.all].T[0]
        # k_harms=np.ones(len(k_harms))
        mi=1
        k=1

        adjs_d=M.faces.bridge_adjacencies(fd,2,3)
        adjs_d0=adjs_d[:,0]
        adjs_d1=adjs_d[:,1]

        M_d=np.zeros((len(fd),sb.nv+sb.nfi))
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
        lines=[]
        cols=[]
        data=[]

        M_d[fd_l,higher_coord]=k_harms[fd]
        M_d[fd_l,lower_coord]=-k_harms[fd]
        ###########################3
        lines.append(fd_l)
        cols.append(higher_coord)
        data.append(k_harms[fd])

        lines.append(fd_l)
        cols.append(lower_coord)
        data.append(-k_harms[fd])
        ##########################
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
            #
            M_d[fd_l[i],sb.nv+fd_l_orig[i]]=1+k_harms[sup_inf[i]].sum()*self.mi_l/(dy*dy)
            M_d[fd_l[i],sb.nv+sup_inf_l[i]]=-k_harms[sup_inf[i]]*self.mi_l/(dy*dy)
            # import pdb; pdb.set_trace()
            ##########################
            lines.append([fd_l[i]])
            cols.append([sb.nv+fd_l_orig[i]])
            data.append([1+k_harms[fd[i]]*self.mi_l/(dy*dy)])

            lines.append(np.repeat(fd_l[i],len(sup_inf_l[i])))
            cols.append(sb.nv+sup_inf_l[i])
            data.append(-k_harms[sup_inf[i]]*self.mi_l/(dy*dy))

            ##############################
            try:
                if len(np.concatenate(cols))!=len(np.concatenate(lines)):
                    print("erro tamanho diferente no momento geral")
                    import pdb; pdb.set_trace()
            except:
                print("erro ao concatenar")
                import pdb; pdb.set_trace()
            if len(laterais[i])==1:
                M_d[fd_l[i],sb.nv+fd_l_orig[i]]+=k_harms[laterais[i]]*self.mi_l/(dx*dx)
                M_d[fd_l[i],sb.nv+laterais_l[i]]-=k_harms[laterais[i]]*self.mi_l/(dx*dx)
                ######################
                lines.append([fd_l[i]])
                cols.append([sb.nv+fd_l_orig[i]])
                data.append(np.array(k_harms[laterais[i]]*self.mi_l/(dx*dx)))

                lines.append(np.array([fd_l[i]]))
                cols.append(np.array(sb.nv+laterais_l[i]))
                data.append(np.array(-k_harms[laterais[i]]*self.mi_l/(dx*dx)))
                try:
                    # print(len(np.concatenate(lines)),len(np.concatenate(cols)),len(np.concatenate(data)))
                    if len(np.concatenate(cols))!=len(np.concatenate(lines)):
                        print("erro tamanho diferente")
                        import pdb; pdb.set_trace()
                except:
                    print("erro ao concatenar em 1")
                    import pdb; pdb.set_trace()
                ##################################
            else:
                M_d[fd_l[i],sb.nv+fd_l_orig[i]]+=(k_harms[laterais[i]]*self.mi_l/(dx*dx)).sum()
                M_d[fd_l[i],sb.nv+laterais_l[i]]-=k_harms[laterais[i]]*self.mi_l/(dx*dx)

                ######################
                lines.append([fd_l[i]])
                cols.append([sb.nv+fd_l_orig[i]])
                data.append([(k_harms[laterais[i]]*self.mi_l/(dx*dx)).sum()])

                lines.append(np.repeat(fd_l[i], len(laterais_l[i])))
                cols.append(sb.nv+laterais_l[i])
                data.append(-k_harms[laterais[i]]*self.mi_l/(dx*dx))

                try:
                    if len(np.concatenate(cols))!=len(np.concatenate(lines)):
                        print("erro tamanho diferente no momento")
                        import pdb; pdb.set_trace()
                except:
                    print("erro ao concatenar")
                    import pdb; pdb.set_trace()

                # import pdb; pdb.set_trace()
                # print(len(np.concatenate(lines)),len(np.concatenate(cols)),len(np.concatenate(data)))
                ######################

        lcd=[np.concatenate(lines),np.concatenate(cols),np.concatenate(data)]
        return lcd

    def get_rhs(self,sb):
        rhs=np.zeros(sb.nv+sb.nfi)
        rhs[0]=1
        rhs[sb.nv-1]=0
        return rhs
