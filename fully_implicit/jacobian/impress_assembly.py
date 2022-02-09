import numpy as np
class assembly():
    def __init__(self,M):
        self.get_jacobian_matrix(M)
    def get_jacobian_matrix(self,M):
        # GID_volumes=M.mb.tag_get_data(M.GLOBAL_ID_tag,M.all_volumes, flat=True)
        GID_volumes=M.volumes.all
        n=len(GID_volumes)
        count=0
        J=np.zeros((2*n,2*n),dtype=float)
        q=np.zeros(2*n)
        import pdb; pdb.set_trace()
        Swns=M.mb.tag_get_data(M.Swns_tag,M.all_volumes,flat=True)
        Swn1s=M.mb.tag_get_data(M.Swn1s_tag,M.all_volumes,flat=True)
        vols=np.repeat(1.0,len(M.all_volumes)) #Volume de cada um dos volumes
        n=len(M.all_volumes)

        pv=M.mb.tag_get_data(M.n_pressure_tag,M.all_volumes,flat=True)
        M.internal_faces=np.setdiff1d(M.all_faces, self.boundary_faces)
        ID_vol=GID_volumes
        symbolic_J=F_Jacobian()
        c_o=symbolic_J.c_o
        c_w=symbolic_J.c_w

        lines=[]
        cols=[]
        data=[]
        lines.append(ID_vol)
        cols.append(n+ID_vol)
        data.append(sym.lambdify((Dx,Dy,phi,Dt),c_o)(self.Dx,self.Dy,0.3,np.repeat(self.dt,len(M.all_volumes))))
        # J[ID_vol][n+ID_vol]+=float(F_Jacobian().c_o.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt}))

        lines.append(n+ID_vol)
        cols.append(n+ID_vol)
        data.append(sym.lambdify((Dx,Dy,phi,Dt),c_w)(self.Dx,self.Dy,0.3,np.repeat(self.dt,len(M.all_volumes))))
        # J[n+ID_vol][n+ID_vol]+=float(F_Jacobian().c_w.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt}))

        linesq=[]
        dataq=[]
        linesq.append(ID_vol)
        dataq.append(sym.lambdify((Dx,Dy,phi,Dt,Sw,Swn),F_Jacobian().acum_o)(self.Dx,self.Dy,0.3,self.dt,Swns,Swn1s))
        # q[ID_vol]+=float(F_Jacobian().acum_o.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt, Sw:Swns[count], Swn:Swn1s[count]}))

        linesq.append(n+ID_vol)
        dataq.append(sym.lambdify((Dx,Dy,phi,Dt,Sw,Swn),F_Jacobian().acum_w)(self.Dx,self.Dy,0.3,self.dt,Swns,Swn1s))
        # q[n+ID_vol]+=float(F_Jacobian().acum_w.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt, Sw:Swns[count], Swn:Swn1s[count]}))

        fac=M.internal_faces
        Adjs=[]
        for f in fac: Adjs.append(M.mb.get_adjacencies(f,M.dimension))
        Adjs=np.array(Adjs)
        adj0=np.array(Adjs[:,0])
        adj1=np.array(Adjs[:,1])

        ids0=M.mb.tag_get_data(M.GLOBAL_ID_tag,adj0,flat=True)
        ids1=M.mb.tag_get_data(M.GLOBAL_ID_tag,adj1,flat=True)
        ID_vol=ids0
        id_j=ids1

        swns0=Swns[ids0]
        swns1=Swns[ids1]

        press0=M.mb.tag_get_data(M.n_pressure_tag,np.array(adj0),flat=True)
        press1=M.mb.tag_get_data(M.n_pressure_tag,np.array(adj1),flat=True)


        pf0=press0
        pf1=press1

        up0=pf0>pf1
        up1=pf0<=pf1
        swf=np.zeros(len(M.internal_faces))
        swf[up0]=swns0[up0]
        swf[up1]=swns1[up1]

        id_up=np.zeros(len(M.internal_faces),dtype=np.int32)
        id_up[up0]=ids0[up0]
        id_up[up1]=ids1[up1]

        Ts=np.ones(len(M.internal_faces))

        J00=sym.lambdify((T,Sw),self.F_Jacobian[0][0])(Ts,swf)
        # J00=float(self.F_Jacobian[0][0].subs({T:1, Sw:swf}))

        J01=sym.lambdify((T,Sw, p_i, p_j),self.F_Jacobian[0][1])(Ts,swf, pf0, pf1)
        # J01=float(self.F_Jacobian[0][1].subs({T:1, Sw:swf, p_i:pv, p_j:pj}))

        J10=sym.lambdify((T,Sw),self.F_Jacobian[1][0])(Ts,swf)
        # J10=float(self.F_Jacobian[1][0].subs({T:1, Sw:swf}))

        J11=sym.lambdify((T,Sw, p_i, p_j),self.F_Jacobian[1][1])(Ts,swf, pf0, pf1)
        # J11=float(self.F_Jacobian[1][1].subs({T:1, Sw:swf, p_i:pv, p_j:pj}))

        linesq.append(ID_vol)
        dataq.append(-sym.lambdify((T,Sw, p_i, p_j),F_Jacobian().F_o)(Ts,swns1, pf0, pf1))
        linesq.append(id_j)
        dataq.append(-sym.lambdify((T,Sw, p_i, p_j),F_Jacobian().F_o)(Ts,swns1, pf1, pf0))
        # q[ID_vol]-=float(F_Jacobian().F_o.subs({T:1.0, Sw:Swns1[count_fac], p_i:pv, p_j:pj}))

        linesq.append(n+ID_vol)
        dataq.append(-sym.lambdify((T,Sw, p_i, p_j),F_Jacobian().F_w)(Ts,swns1, pf0, pf1))
        linesq.append(n+id_j)
        dataq.append(-sym.lambdify((T,Sw, p_i, p_j),F_Jacobian().F_w)(Ts,swns1, pf1, pf0))
        # q[n+ID_vol]-=float(F_Jacobian().F_w.subs({T:1.0, Sw:Swns1[count_fac], p_i:pv, p_j:pj}))

        lines.append(ID_vol)
        cols.append(ID_vol)
        data.append(-J00)
        lines.append(id_j)
        cols.append(id_j)
        data.append(-J00)
        # J[ID_vol][ID_vol]-=J00

        lines.append(ID_vol)
        cols.append(id_j)
        data.append(J00)
        lines.append(id_j)
        cols.append(ID_vol)
        data.append(J00)
        # J[ID_vol][id_j]+=J00

        lines.append(n+ID_vol)
        cols.append(ID_vol)
        data.append(-J10)
        lines.append(n+id_j)
        cols.append(id_j)
        data.append(-J10)
        # J[n+ID_vol][ID_vol]-=J10

        lines.append(n+ID_vol)
        cols.append(id_j)
        data.append(J10)
        lines.append(n+id_j)
        cols.append(ID_vol)
        data.append(J10)
        # J[n+ID_vol][id_j]+=J10
        #
        lines.append(ID_vol)
        cols.append(n+id_up)
        data.append(-J01)
        lines.append(id_j)
        cols.append(n+id_up)
        data.append(J01)
        # J[ID_vol][n+id_up]-=J01

        lines.append(n+ID_vol)
        cols.append(n+id_up)
        data.append(-J11)
        lines.append(n+id_j)
        cols.append(n+id_up)
        data.append(J11)
        # J[n+ID_vol][n+id_up]-=J11
        lines=np.concatenate(lines)
        cols=np.concatenate(cols)
        data=np.concatenate(data)

        linesq=np.concatenate(linesq)
        dataq=np.concatenate(dataq)

        J=scipy.sparse.csc_matrix((data,(lines,cols)),shape=(2*len(M.all_volumes),2*len(M.all_volumes)))
        q=scipy.sparse.csc_matrix((dataq,(linesq,np.zeros(len(linesq)))),shape=(2*len(M.all_volumes),1))

        q=q.transpose().toarray()[0]
        self.iterac+=1
        return(J, q)
