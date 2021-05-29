import numpy as np
import scipy
from scipy.sparse import find
import sympy as sym

T, S_up, Sw, So, Swn, Son, Dt, k, phi, p_i, p_j, Dx, Dy=sym.symbols("T S Sw So Swn Son Dt k phi p_i p_j Dx Dy")

class assembly():
    def __init__(self,adjs, Ts, data_impress, timestep, wells, F_Jacobian):
        self.wells = wells
        self.ids=data_impress['GID_0']
        self.n=len(self.ids)
        self.nfi=len(adjs)
        self.Dx=1.0
        self.Dy=1.0
        self.Dz=1.0
        self.F_Jacobian=F_Jacobian
        self.dt=timestep
        self.pressure=data_impress['pressure'].copy()
        self.swns=data_impress['swns'].copy()
        self.swn1s=data_impress['swn1s'].copy()
        self.adjs = adjs
        self.Ts = Ts
        self.J,self.q=self.get_jacobian_matrix()

    def apply_BC(self,lines, cols, data, q):
        n=int(len(q)/2)
        wells=self.wells

        q[wells['ws_p']]=0
        q[wells['ws_inj']+n]=0
        if wells['count']==0:
            q[wells['ws_q']]+=wells['values_q']
        for l in wells['ws_p']:
            data[lines==l]=0
            lines=np.append(lines,l)
            cols=np.append(cols,l)
            data=np.append(data,1)
        for l in np.setdiff1d(wells['ws_inj'],wells['ws_q']):
            data[lines==l+n]=0
            lines=np.append(lines,l+n)
            cols=np.append(cols,l+n)
            data=np.append(data,1)
        return lines, cols, data, q

    def get_jacobian_matrix(self):
        n=len(self.ids)
        count=0
        Swns=self.swns
        Swn1s=self.swn1s
        Swns[Swns<0]=0
        Swns[Swns>1]=1
        Swn1s[Swn1s<0]=0
        Swn1s[Swn1s>1]=1
        ID_vol=self.ids
        lines=[]
        cols=[]
        data=[]
        lines.append(ID_vol)
        cols.append(n+ID_vol)
        data.append(self.F_Jacobian.c_o(self.Dx,self.Dy,0.3,np.repeat(self.dt,n)))
        # J[ID_vol][n+ID_vol]+=float(F_Jacobian().c_o.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt}))
        lines.append(n+ID_vol)
        cols.append(n+ID_vol)
        data.append(self.F_Jacobian.c_w(self.Dx,self.Dy,0.3,np.repeat(self.dt,n)))
        # J[n+ID_vol][n+ID_vol]+=float(F_Jacobian().c_w.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt}))
        linesq=[]
        dataq=[]
        linesq.append(ID_vol)
        dataq.append(self.F_Jacobian.acum_o(self.Dx,self.Dy,0.3,self.dt,Swns,Swn1s))
        # q[ID_vol]+=float(F_Jacobian().acum_o.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt, Sw:Swns[count], Swn:Swn1s[count]}))
        linesq.append(n+ID_vol)
        dataq.append(self.F_Jacobian.acum_w(self.Dx,self.Dy,0.3,self.dt,Swns,Swn1s))
        # q[n+ID_vol]+=float(F_Jacobian().acum_w.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt, Sw:Swns[count], Swn:Swn1s[count]}))
        Adjs=np.array(self.adjs)
        adj0=np.array(Adjs[:,0])
        adj1=np.array(Adjs[:,1])
        ids0=self.ids[adj0]
        ids1=self.ids[adj1]
        ID_vol=ids0
        id_j=ids1
        swns0=Swns[ids0]
        swns1=Swns[ids1]
        press0=self.pressure[adj0]
        press1=self.pressure[adj1]
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
        Ts=self.Ts

        J00=self.F_Jacobian.J[0][0](Ts,swf)
        # J00=float(self.F_Jacobian[0][0].subs({T:1, Sw:swf}))
        J01=self.F_Jacobian.J[0][1](Ts,swf, pf0, pf1)
        # J01=float(self.F_Jacobian[0][1].subs({T:1, Sw:swf, p_i:pv, p_j:pj}))
        J10=self.F_Jacobian.J[1][0](Ts,swf)
        # J10=float(self.F_Jacobian[1][0].subs({T:1, Sw:swf}))
        J11=self.F_Jacobian.J[1][1](Ts,swf, pf0, pf1)
        # J11=float(self.F_Jacobian[1][1].subs({T:1, Sw:swf, p_i:pv, p_j:pj}))
        linesq.append(ID_vol)
        dataq.append(-self.F_Jacobian.F_o(Ts,swf, pf0, pf1))
        linesq.append(id_j)
        dataq.append(-self.F_Jacobian.F_o(Ts,swf, pf1, pf0))
        # q[ID_vol]-=float(F_Jacobian().F_o.subs({T:1.0, Sw:Swns1[count_fac], p_i:pv, p_j:pj}))
        linesq.append(n+ID_vol)
        dataq.append(-self.F_Jacobian.F_w(Ts,swf, pf0, pf1))
        linesq.append(n+id_j)
        dataq.append(-self.F_Jacobian.F_w(Ts,swf, pf1, pf0))
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
        q=np.bincount(linesq, weights=dataq)
        lines, cols, data, q = self.apply_BC(lines, cols, data, q)
        J=scipy.sparse.csc_matrix((data,(lines,cols)),shape=(2*n,2*n))
        return(J, q)
