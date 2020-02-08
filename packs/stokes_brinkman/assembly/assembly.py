import numpy as np

class global_assembly:
    def __init__(self, M, dx, dy, dz):
        self.nv=len(M.volumes.all)
        self.nfi=len(M.faces.internal)
        self.mi_l=0.1
        self.M_cont=self.get_continuity_matrix(M, dx, dy, dz)
        self.M_momentum_y=self.get_momentum_y_matrix(M,dx,dy)
        self.M_momentum_x=self.get_momentum_x_matrix(M,dx,dy)
        self.M=np.vstack([self.M_cont,self.M_momentum_x,self.M_momentum_y])


    def get_continuity_matrix(self, M, dx, dy, dz):

        M_cont=np.zeros((self.nv,self.nv+self.nfi))
        '''
        volumes=M.volumes.all
        adjs=M.volumes.bridge_adjacencies(volumes,2,2)
        for f in M.faces.boundary:
            adjs[adjs==f]=-1
        for aa in range(len(adjs)):
            all_faces_vol=adjs[aa][adjs[aa,:]!=-1]
            us=[]
            vs=[]
            for f in all_faces_vol:
                verts=M.faces.bridge_adjacencies(f,0,0)
                coords=M.nodes.coords[verts]
                if abs(coords[:,0].max()-coords[:,0].min())<10**-6:
                    us.append(f)
                else:
                    vs.append(f)
            xs=[]
            for u in us:
                verts=M.faces.bridge_adjacencies(u,0,0)
                xs.append(M.nodes.coords[verts][:,0][0])

            usl=M.id_fint[us].T[0]
            vsl=M.id_fint[vs].T[0]

            if len(us)==2:
                if xs[0]>xs[1]:
                    M_cont[volumes[aa],self.nv+usl[1]]=1
                    M_cont[volumes[aa],self.nv+usl[0]]=-1
                else:
                    M_cont[volumes[aa],self.nv+usl[1]]=-1
                    M_cont[volumes[aa],self.nv+usl[0]]=1
            else:
                if abs(xs[0]-dx)<0.6:
                    M_cont[volumes[aa],self.nv+usl[0]]=1
                else:
                    M_cont[volumes[aa],self.nv+usl[0]]=-1

            ys=[]
            for v in vs:
                verts=M.faces.bridge_adjacencies(v,0,0)
                ys.append(M.nodes.coords[verts][:,1][0])
            if len(vs)==2:
                if ys[0]>ys[1]:
                    M_cont[volumes[aa],self.nv+vsl[1]]=1
                    M_cont[volumes[aa],self.nv+vsl[0]]=-1
                else:
                    M_cont[volumes[aa],self.nv+vsl[1]]=-1
                    M_cont[volumes[aa],self.nv+vsl[0]]=1
            else:

                if abs(ys[0]-dy)<0.6:
                    M_cont[volumes[aa],self.nv+vsl[0]]=1
                else:
                    M_cont[volumes[aa],self.nv+vsl[0]]=-1

        '''
        internal_faces=M.faces.internal
        volumes=M.volumes.all
        faces=M.volumes.bridge_adjacencies(volumes,2,2)
        faces_center_x=M.faces.center(faces.flatten())[:,0].reshape(faces.shape)
        faces_center_y=M.faces.center(faces.flatten())[:,1].reshape(faces.shape)

        for i in range(len(faces)):
            facs=faces[i]
            xcents=faces_center_x[i]
            ycents=faces_center_y[i]

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


        #import pdb; pdb.set_trace()





            # import pdb; pdb.set_trace()
        # internal_faces=M.faces.internal
        # adjs=M.faces.bridge_adjacencies(internal_faces,2,3)
        # adjs0=adjs[:,0]
        # adjs1=adjs[:,1]
        # M_cont[adjs0,nfi+adjs1]=1
        # M_cont[adjs1,nfi+adjs0]=-1
        #
        # M_cont[adjs1,nv+adjs0]=-1
        # M_cont[adjs0,nv+adjs1]=1
        # np.savetxt("results/M_cont.csv",M_cont,delimiter=",")


        return(M_cont)

    def get_momentum_y_matrix(self,M,dx,dy):
        T=1
        mi=1
        k=1

        faces=M.faces.internal
        nodes=M.faces.bridge_adjacencies(faces,0,0).flatten()
        coords=M.nodes.coords[nodes].reshape(len(faces),4,3)
        deltas_x_fac=coords[:,:,0].max(axis=1)-coords[:,:,0].min(axis=1)
        horizontal=faces[deltas_x_fac>10**-8]
        vertical=faces[deltas_x_fac<10**-8]
        adjs_h=M.faces.bridge_adjacencies(vertical,2,3)
        adjs_h0=adjs_h[:,0]
        adjs_h1=adjs_h[:,1]
        M_y=np.zeros((len(vertical),self.nv+self.nfi))
        c0=M.volumes.center(adjs_h0)[:,0]
        c1=M.volumes.center(adjs_h1)[:,0]

        higher_coord=np.zeros(len(vertical),dtype=int)
        higher_coord[c0>c1]=adjs_h0[c0>c1]
        higher_coord[c0<c1]=adjs_h1[c0<c1]

        lower_coord=np.zeros(len(vertical),dtype=int)
        lower_coord[c0<c1]=adjs_h0[c0<c1]
        lower_coord[c0>c1]=adjs_h1[c0>c1]

        h_l=M.id_fint[horizontal].T[0]
        M_y[h_l,higher_coord]=1
        M_y[h_l,lower_coord]=-1


        fac_viz_ares = M.faces.bridge_adjacencies(horizontal,1,2)
        laterais=[]
        for fac in fac_viz_ares:
            laterais.append(np.intersect1d(fac,horizontal))

        laterais_l=[]
        for l in laterais:
            laterais_l.append(M.id_fint[l].T[0])

        faces_viz_adjs_h=M.volumes.bridge_adjacencies(adjs_h.flatten(),2,2).reshape(len(vertical),12)

        sup_inf=[]
        i=0
        for fac in faces_viz_adjs_h:
            # import pdb; pdb.set_trace()
            sup_inf.append(np.setdiff1d(np.intersect1d(fac,vertical),vertical[i]))
            i+=1
        sup_inf_l=[]
        for si in sup_inf:
            sup_inf_l.append(M.id_fint[si].T[0])
        vertical_l=M.id_fint[vertical]-len(horizontal)

        for i in range(len(vertical)):
            M_y[vertical_l[i],self.nv+len(horizontal)+vertical_l[i]]=1-self.mi_l/(dy*dy)
            M_y[vertical_l[i],self.nv+sup_inf_l[i]]=self.mi_l/(dy*dy)

            if len(laterais[i])==1:
                M_y[vertical_l[i],self.nv+len(horizontal)+vertical_l[i]]-=1*self.mi_l/(dx*dx)
                M_y[vertical_l[i],self.nv+len(horizontal)+laterais_l[i]]+=1*self.mi_l/(dx*dx)
            else:
                M_y[vertical_l[i],self.nv+len(horizontal)+vertical_l[i]]-=self.mi_l/(dx*dx)
                M_y[vertical_l[i],self.nv+len(horizontal)+laterais_l[i]]+=self.mi_l/(dx*dx)

        return M_y

    def get_momentum_x_matrix(self,M,dx,dy):
        T=1
        mi=1
        k=1

        faces=M.faces.internal
        nodes=M.faces.bridge_adjacencies(faces,0,0).flatten()
        coords=M.nodes.coords[nodes].reshape(len(faces),4,3)
        deltas_x_fac=coords[:,:,0].max(axis=1)-coords[:,:,0].min(axis=1)
        horizontal=faces[deltas_x_fac>10**-8]
        vertical=faces[deltas_x_fac<10**-8]

        horizontal, vertical = vertical,horizontal
        adjs_h=M.faces.bridge_adjacencies(vertical,2,3)
        adjs_h0=adjs_h[:,0]
        adjs_h1=adjs_h[:,1]
        M_y=np.zeros((len(vertical),self.nv+self.nfi))
        c0=M.volumes.center(adjs_h0)[:,1]
        c1=M.volumes.center(adjs_h1)[:,1]

        higher_coord=np.zeros(len(vertical),dtype=int)
        higher_coord[c0>c1]=adjs_h0[c0>c1]
        higher_coord[c0<c1]=adjs_h1[c0<c1]

        lower_coord=np.zeros(len(vertical),dtype=int)
        lower_coord[c0<c1]=adjs_h0[c0<c1]
        lower_coord[c0>c1]=adjs_h1[c0>c1]

        h_l=M.id_fint[horizontal].T[0]-len(vertical)
        # import pdb; pdb.set_trace()
        M_y[h_l,higher_coord]=1
        M_y[h_l,lower_coord]=-1


        fac_viz_ares = M.faces.bridge_adjacencies(horizontal,1,2)
        laterais=[]
        for fac in fac_viz_ares:
            laterais.append(np.intersect1d(fac,horizontal))

        laterais_l=[]
        for l in laterais:
            laterais_l.append(M.id_fint[l].T[0])

        faces_viz_adjs_h=M.volumes.bridge_adjacencies(adjs_h.flatten(),2,2).reshape(len(vertical),12)

        sup_inf=[]
        i=0
        for fac in faces_viz_adjs_h:
            # import pdb; pdb.set_trace()
            sup_inf.append(np.setdiff1d(np.intersect1d(fac,vertical),vertical[i]))
            i+=1
        sup_inf_l=[]
        for si in sup_inf:
            sup_inf_l.append(M.id_fint[si].T[0])
        vertical_l=M.id_fint[vertical]

        for i in range(len(vertical)):

            M_y[vertical_l[i],self.nv+vertical_l[i]]=1-self.mi_l/(mi*dy*dy)
            M_y[vertical_l[i],self.nv+sup_inf_l[i]]=self.mi_l/(dy*dy)

            if len(laterais[i])==1:
                M_y[vertical_l[i],self.nv+vertical_l[i]]-=self.mi_l/(dx*dx)
                M_y[vertical_l[i],self.nv+laterais_l[i]-len(horizontal)]=self.mi_l/(dx*dx)
            else:
                M_y[vertical_l[i],self.nv+vertical_l[i]]-=self.mi_l/(dx*dx)
                M_y[vertical_l[i],self.nv+laterais_l[i]-len(horizontal)]=self.mi_l/(dx*dx)

        return M_y
