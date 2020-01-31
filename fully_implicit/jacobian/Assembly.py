from fully_implicit.jacobian.symbolic_jacobian import symbolic_J as F_Jacobian
import scipy
from scipy.sparse import linalg
import sympy as sym
import numpy as np

# import pandas as pd
from sympy import Heaviside as H

T, S_up, Sw, So, Swn, Son, Dt, k, phi, p_i, p_j, Dx, Dy=sym.symbols("T S Sw So Swn Son Dt k phi p_i p_j Dx Dy")
class Ass:
    def __init__(self,M):
        self.internal_faces=[f for f in M.all_faces if len(M.mb.get_adjacencies(f,M.dimension))==2]

        self.boundary_faces=np.setdiff1d(M.all_faces,self.internal_faces)

        self.F_Jacobian=F_Jacobian().J
        # self.set_properties(M)

        self.vazao=1.0
        self.porous_volume=1*0.3  #Volume*porosidade
        self.Dx=1.0/6
        self.Dy=1.0/6
        self.Dx_av=np.repeat(1.0/6,len(M.all_volumes))
        self.Dy_av=np.repeat(1.0/6,len(M.all_volumes))
        self.Dz_av=np.repeat(1.0/6,len(M.all_volumes))
        self.pvi_lim=1.0

        for j in range(8):
            self.pvi_def=0.15-0.02*j
            ni=int(self.pvi_lim/self.pvi_def)
            self.qov=[]
            self.qwv=[]
            self.vpiv=[]
            self.pvi_acum=0.0
            self.pvi=0.02
            self.dt=self.porous_volume*self.pvi/self.vazao
            self.set_properties(M)
            self.i=0
            self.iterac=0
            for i in range(ni):
                self.Newton(M)
                self.i+=1
            filepath = 'results/production_data.xlsx'
            try:
                old=numpy.savetxt(filepath,self.F_Jacobian,delimiter=",")
                data=np.hstack([old,self.data])
                print('Appended to an old file')
            except:
                print('First Data file')
                data=self.data
                self.iterac+=1
            df = pd.DataFrame (data)
            df.to_excel(filepath, index=False)

    def Newton(self,M):
        n=len(M.all_volumes)
        S0=M.mb.tag_get_data(M.Swns_tag,M.all_volumes,flat=True)
        P0=M.mb.tag_get_data(M.n_pressure_tag,M.all_volumes,flat=True)
        SP=np.concatenate([S0,P0])
        gids=M.mb.tag_get_data(M.GLOBAL_ID_tag,M.all_volumes,flat=True)
        max_iter=8
        if self.i==0:
            max_iter*=3
        i=0

        while i<max_iter and (i==0 or max(abs(sol[n:]))>0.001):
            if i>4 or self.i>0:
                self.pvi=self.pvi_def
                self.dt=self.porous_volume*self.pvi/self.vazao
            J, q=self.get_jacobian_matrix(M)
            J[0]=0.0
            J[0][0]=1.0

            J[n,:]=0.0
            J[n,n]=1.0
            q[n]=0.0

            q[0]=0.0
            # q[-1]=q[-1].toarray()[0][0]-self.vazao
            q[-1]-=self.vazao
            # J=scipy.sparse.csc_matrix(J)
            Jpp=J[0:n,0:n]
            Jps=J[0:n,n:]
            Jsp=J[n:,0:n]
            Jss=J[n:,n:]

            Fp=-q[0:n]
            Fs=-q[n:]

            # J=scipy.sparse.csc_matrix(J)
            # Jpp1=J[0:n,0:n]
            # Jps1=J[0:n,n:]
            # Jsp1=J[n:,0:n]
            # Jss1=J[n:,n:]

            # Mp1=Jpp1-Jps1*linalg.spsolve(Jss1,Jsp1)
            Mp=Jpp-np.dot(Jps,np.linalg.solve(Jss,Jsp))

            # qp1=(scipy.sparse.csc_matrix(Fp).T-scipy.sparse.csc_matrix(Jps1*linalg.spsolve(Jss1,Fs)).T).toarray().transpose()[0]
            qp=Fp-np.dot(Jps,np.linalg.solve(Jss,Fs))

            # dp1=linalg.spsolve(Mp1,qp1)
            dp=np.linalg.solve(Mp,qp)

            # ds1=linalg.spsolve(Jss1,Fs-Jsp1*dp1)
            ds=np.linalg.solve(Jss,Fs-np.dot(Jsp,dp))
            # import pdb; pdb.set_trace()
            sol=np.concatenate([dp,ds])

            S0+=sol[n:]
            pos_men=np.where(S0<0.2)[0]
            pos_mai=np.where(S0>0.8)[0]
            S0[pos_men]=0.2
            S0[pos_mai]=0.8
            P0+=sol[0:n]
            print(max(abs(sol[0:n])),max(abs(sol[n:])))

            M.mb.tag_set_data(M.n_pressure_tag,M.all_volumes,P0[gids])
            M.mb.tag_set_data(M.Swns_tag,M.all_volumes,S0[gids])
            i+=1

        # Traçando curva de produção
        v0_facs=M.mb.get_adjacencies(M.all_volumes[0], M.dimension-1)
        v0_internal_faces=np.setdiff1d(v0_facs, self.boundary_faces)
        pv=M.mb.tag_get_data(M.n_pressure_tag,M.all_volumes[0],flat=True)[0]

        qo=0.0
        qw=0.0

        for f in v0_internal_faces:

            adjs=M.mb.get_adjacencies(f,M.dimension)
            Sw_up=max(M.mb.tag_get_data(M.Swns_tag,np.array(adjs),flat=True))
            pj=max(M.mb.tag_get_data(M.n_pressure_tag,adjs,flat=True))
            qo+=float(F_Jacobian().F_o.subs({T:1.0, k:1.0, Sw:Sw_up, p_i:pv, p_j:pj}))/0.8
            qw+=float(F_Jacobian().F_w.subs({T:1.0, k:1.0, Sw:Sw_up, p_i:pv, p_j:pj}))
        self.pvi_acum+=self.pvi
        self.qov.append(qo)
        self.qwv.append(qw)
        self.vpiv.append(self.pvi_acum)

        self.data=np.array([self.vpiv,self.qov,self.qwv])

        print('oleo, agua',qo,qw,qo+qw,'pvi',self.pvi_acum)
        m4 = M.mb.create_meshset()
        M.mb.add_entities(m4, M.all_volumes)
        M.mb.write_file('results/biphasic/fully_implicit'+str(self.i)+'.vtk',[m4])
        M.mb.tag_set_data(M.Swn1s_tag,M.all_volumes,S0[gids])


    def get_jacobian_matrix(self,M):
        GID_volumes=M.mb.tag_get_data(M.GLOBAL_ID_tag,M.all_volumes, flat=True)
        count=0
        J=np.zeros((2*len(M.all_volumes),2*len(M.all_volumes)),dtype=float)
        q=np.zeros(2*len(M.all_volumes))
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
        J=J.toarray()
        # np.savetxt("results/vec.csv",JJ.toarray(),delimiter=",")

        self.iterac+=1
        return(J, q)

    def set_properties(self,M):
        M.mb.tag_set_data(M.GLOBAL_ID_tag,M.all_volumes,range(len(M.all_volumes)))
        M.mb.tag_set_data(M.k_eq_tag,self.internal_faces,np.repeat(1.0,len(self.internal_faces)))
        M.mb.tag_set_data(M.phi_tag,M.all_volumes,np.repeat(0.3,len(M.all_volumes)))
        M.mb.tag_set_data(M.press_value_tag,np.uint64([M.all_volumes[0],M.all_volumes[-1]]),[0.0,15.0])
        M.mb.tag_set_data(M.Swns_tag,M.all_volumes,np.repeat(0.2,len(M.all_volumes)))
        M.mb.tag_set_data(M.Swns_tag,M.all_volumes[-1],0.8)
        M.mb.tag_set_data(M.Swn1s_tag,M.all_volumes,np.repeat(0.2,len(M.all_volumes)))
        M.mb.tag_set_data(M.Swn1s_tag,M.all_volumes[-1],0.8)
        M.mb.tag_set_data(M.n_pressure_tag,M.all_volumes,np.array(range(len(M.all_volumes)),dtype=np.float))
