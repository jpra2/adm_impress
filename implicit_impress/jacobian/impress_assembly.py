import numpy as np
from implicit_impress.jacobian.symbolic_jacobian import symbolic_J as s_J
import scipy
from scipy.sparse import find
import sympy as sym

T, S_up, Sw, So, Swn, Son, Dt, k, phi, p_i, p_j, Dx, Dy=sym.symbols("T S Sw So Swn Son Dt k phi p_i p_j Dx Dy")
class assembly():
    def __init__(self,M,timestep, wells):
        self.n=len(M.volumes.all)
        self.nfi=len(M.faces.internal)
        self.vazao=10000.0
        self.Dx=1.0
        self.Dy=1.0
        self.Dz=1.0
        # M.swn1s[:]=M.swns[:]
        self.F_Jacobian=s_J()
        self.dt=timestep
        J,q=self.get_jacobian_matrix(M)
        self.J, self.q= self.apply_BC(J, q, wells)

        # self.J=self.apply_dirichlet(J,[0])#,self.n])
        # self.q=q
        # if abs(J2-self.J).sum() >0:
        #     import pdb; pdb.set_trace()
        # if abs(q2-self.q).sum() != abs(self.vazao):
        #     import pdb; pdb.set_trace()

    def apply_BC(self,J3,q3,wells):
        J2=J3.tolil()
        q2=q3.copy()

        J2[wells['ws_p'],:]=0
        J2[wells['ws_p'], wells['ws_p']]=1
        q2[wells['ws_p']]=0
        # J2[wells['ws_p']+self.n,:]=0
        # J2[wells['ws_p']+self.n, wells['ws_p']+self.n]=1
        # q2[wells['ws_p']+self.n]=0
        J2[wells['ws_inj']+self.n,:]=0
        J2[wells['ws_inj']+self.n, wells['ws_inj']+self.n]=1
        q2[wells['ws_inj']+self.n]=0

        q2[self.n+wells['ws_q']]-=wells['values_q']
        J2=J2.tocsc()
        return J2,q2

    def get_jacobian_matrix(self,M):
        GID_volumes=M.volumes.all
        n=len(GID_volumes)
        count=0

        Swns=M.swns[:].transpose()[0]
        Swn1s=M.swn1s[:].transpose()[0]
        Swns[Swns<0]=0
        Swns[Swns>1]=1
        Swn1s[Swn1s<0]=0
        Swn1s[Swn1s>1]=1
        # import pdb; pdb.set_trace()
        # vols=np.repeat(1.0,self.n) #Volume de cada um dos volumes

        ID_vol=GID_volumes
        symbolic_J=s_J()
        c_o=symbolic_J.c_o
        c_w=symbolic_J.c_w

        lines=[]
        cols=[]
        data=[]
        lines.append(ID_vol)
        cols.append(n+ID_vol)
        data.append(sym.lambdify((Dx,Dy,phi,Dt),c_o)(self.Dx,self.Dy,0.3,np.repeat(self.dt,self.n)))
        # J[ID_vol][n+ID_vol]+=float(F_Jacobian().c_o.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt}))

        lines.append(n+ID_vol)
        cols.append(n+ID_vol)
        data.append(sym.lambdify((Dx,Dy,phi,Dt),c_w)(self.Dx,self.Dy,0.3,np.repeat(self.dt,self.n)))
        # J[n+ID_vol][n+ID_vol]+=float(F_Jacobian().c_w.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt}))

        linesq=[]
        dataq=[]
        linesq.append(ID_vol)

        dataq.append(sym.lambdify((Dx,Dy,phi,Dt,Sw,Swn),self.F_Jacobian.acum_o)(self.Dx,self.Dy,0.3,self.dt,Swns,Swn1s))
        # q[ID_vol]+=float(F_Jacobian().acum_o.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt, Sw:Swns[count], Swn:Swn1s[count]}))

        linesq.append(n+ID_vol)
        dataq.append(sym.lambdify((Dx,Dy,phi,Dt,Sw,Swn),self.F_Jacobian.acum_w)(self.Dx,self.Dy,0.3,self.dt,Swns,Swn1s))
        # q[n+ID_vol]+=float(F_Jacobian().acum_w.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt, Sw:Swns[count], Swn:Swn1s[count]}))

        # fac=M.faces.internal
        # M.core
        Adjs=M.faces.bridge_adjacencies(M.faces.internal_elements[:],2,3)
        Adjs=np.array(Adjs)
        adj0=np.array(Adjs[:,0])
        adj1=np.array(Adjs[:,1])

        ids0=M.volumes.all[adj0]
        ids1=M.volumes.all[adj1]
        ID_vol=ids0
        id_j=ids1

        swns0=Swns[ids0]
        swns1=Swns[ids1]

        press0=M.pressure[adj0].transpose()[0]
        press1=M.pressure[adj1].transpose()[0]

        pf0=press0
        pf1=press1

        up0=pf0>pf1
        up1=pf0<=pf1
        swf=np.zeros(self.nfi)

        swf[up0]=swns0[up0]
        swf[up1]=swns1[up1]

        id_up=np.zeros(self.nfi,dtype=np.int32)
        id_up[up0]=ids0[up0]
        id_up[up1]=ids1[up1]

        Ts=M.k_harm[M.faces.internal].T[0]

        J00=sym.lambdify((T,Sw),self.F_Jacobian.J[0][0])(Ts,swf)
        # J00=float(self.F_Jacobian[0][0].subs({T:1, Sw:swf}))

        J01=sym.lambdify((T,Sw, p_i, p_j),self.F_Jacobian.J[0][1])(Ts,swf, pf0, pf1)
        # J01=float(self.F_Jacobian[0][1].subs({T:1, Sw:swf, p_i:pv, p_j:pj}))

        J10=sym.lambdify((T,Sw),self.F_Jacobian.J[1][0])(Ts,swf)
        # J10=float(self.F_Jacobian[1][0].subs({T:1, Sw:swf}))

        J11=sym.lambdify((T,Sw, p_i, p_j),self.F_Jacobian.J[1][1])(Ts,swf, pf0, pf1)
        # J11=float(self.F_Jacobian[1][1].subs({T:1, Sw:swf, p_i:pv, p_j:pj}))

        linesq.append(ID_vol)
        dataq.append(-sym.lambdify((T,Sw, p_i, p_j),self.F_Jacobian.F_o)(Ts,swf, pf0, pf1))
        linesq.append(id_j)
        dataq.append(-sym.lambdify((T,Sw, p_i, p_j),self.F_Jacobian.F_o)(Ts,swf, pf1, pf0))
        # q[ID_vol]-=float(F_Jacobian().F_o.subs({T:1.0, Sw:Swns1[count_fac], p_i:pv, p_j:pj}))

        linesq.append(n+ID_vol)
        dataq.append(-sym.lambdify((T,Sw, p_i, p_j),self.F_Jacobian.F_w)(Ts,swf, pf0, pf1))
        linesq.append(n+id_j)
        dataq.append(-sym.lambdify((T,Sw, p_i, p_j),self.F_Jacobian.F_w)(Ts,swf, pf1, pf0))
        # q[n+ID_vol]-=float(F_Jacobian().F_w.subs({T:1.0, Sw:Swns1[count_fac], p_i:pv, p_j:pj}))
        # import pdb; pdb.set_trace()
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

        J=scipy.sparse.csc_matrix((data,(lines,cols)),shape=(2*self.n,2*self.n))
        q=scipy.sparse.csc_matrix((dataq,(linesq,np.zeros(len(linesq)))),shape=(2*self.n,1))

        q=q.transpose().toarray()[0]

        #
        # q[0]=0.0
        # q[-1]-=self.vazao
        # self.iterac+=1
        return(J, q)

    def apply_dirichlet(self, Mat,l):
        ff=find(Mat)
        lines=ff[0]
        cols=ff[1]
        data=ff[2]
        l=np.unique(l)
        for ls in l:
            data[lines==ls]=0
            data[(lines==ls) & (cols==ls)]=1
        return scipy.sparse.csc_matrix((data, (lines,cols)),shape=Mat.shape)
