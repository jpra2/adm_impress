import numpy as np

class global_assembly:
    def __init__(self, M, dx, dy, dz):
        self.nv=len(M.volumes.all)
        self.nfi=len(M.faces.internal)
        self.M_cont=self.get_continuity_matrix(M, dx, dy, dz)
        self.M_momentum=self.get_momentum_x_matrix(M,dx,dy)

    def get_continuity_matrix(self, M, dx, dy, dz):
        M_cont=np.zeros((self.nv,self.nv+self.nfi))
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
            # import pdb; pdb.set_trace()
            if len(us)==2:
                if xs[0]>xs[1]:
                    M_cont[volumes[aa],self.nv+us[1]-30]=1
                    M_cont[volumes[aa],self.nv+us[0]-30]=-1
                else:
                    M_cont[volumes[aa],self.nv+us[1]-30]=-1
                    M_cont[volumes[aa],self.nv+us[0]-30]=1
            else:
                if abs(xs[0]-dx)<0.6:
                    M_cont[volumes[aa],self.nv+us[0]-30]=1
                else:
                    M_cont[volumes[aa],self.nv+us[0]-30]=-1

            ys=[]
            for v in vs:
                verts=M.faces.bridge_adjacencies(v,0,0)
                ys.append(M.nodes.coords[verts][:,1][0])
            if len(vs)==2:
                if ys[0]>ys[1]:
                    M_cont[volumes[aa],self.nv+vs[1]-30]=1
                    M_cont[volumes[aa],self.nv+vs[0]-30]=-1
                else:
                    M_cont[volumes[aa],self.nv+vs[1]-30]=-1
                    M_cont[volumes[aa],self.nv+vs[0]-30]=1
            else:

                if abs(ys[0]-dy)<0.6:
                    M_cont[volumes[aa],self.nv+vs[0]-30]=1
                else:
                    M_cont[volumes[aa],self.nv+vs[0]-30]=-1
        return(M_cont)

                # import pdb; pdb.set_trace()


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
        np.savetxt("results/M_cont.csv",M_cont,delimiter=",")

    def get_momentum_x_matrix(self,M,dx,dy):
        T=1
        mi=1

        faces=M.faces.internal
        nodes=M.faces.bridge_adjacencies(faces,0,0).flatten()
        coords=M.nodes.coords[nodes].reshape(len(faces),4,3)
        deltas_x_fac=coords[:,:,0].max(axis=1)-coords[:,:,0].min(axis=1)
        horizontal=faces[deltas_x_fac>10**-8]
        vertical=faces[deltas_x_fac<10**-8]
        adjs_h=M.faces.bridge_adjacencies(vertical,2,3)
        adjs_h0=adjs_h[:,0]
        adjs_h1=adjs_h[:,1]
        M_x=np.zeros((len(vertical),self.nv+self.nfi))
        c0=M.volumes.center(adjs_h0)[:,0]
        c1=M.volumes.center(adjs_h1)[:,0]

        higher_coord=np.zeros(len(vertical),dtype=int)
        higher_coord[c0>c1]=adjs_h0[c0>c1]
        higher_coord[c0<c1]=adjs_h1[c0<c1]

        lower_coord=np.zeros(len(vertical),dtype=int)
        lower_coord[c0<c1]=adjs_h0[c0<c1]
        lower_coord[c0>c1]=adjs_h1[c0>c1]

        import pdb; pdb.set_trace()
        M_x[self.nv+horizontal-30,higher_coord]=1
        M_x[self.nv+horizontal-30,lower_coord]=-1

        import pdb; pdb.set_trace()



        # import pdb; pdb.set_trace()
