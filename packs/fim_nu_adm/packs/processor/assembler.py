from packs.preprocessor.symbolic_calculation import symbolic_J as s_J
import scipy.sparse as sp
import numpy as np
import time

class Assembler:
    def __init__(self, wells, faces, volumes, CFL=10):
        self.wells=wells
        self.F_Jacobian=s_J()
        self.adjs=faces['adjacent_volumes']
        self.centroids=volumes['centroids']
        self.Ts=faces['permeabilities']
        self.GID_0=volumes['GID_0']
        self.gh=((self.centroids[self.adjs[:,0]]-self.centroids[self.adjs[:,1]])*self.F_Jacobian.g)[:,np.array(self.F_Jacobian.g)!=0].T[0]

    def get_jacobian_matrix(self, Swns, Swn1s, p, time_step):
        # Ts, adjs, swns, swn1s, time_step, wells, F_Jacobian
        Swns[Swns<0]=0
        Swns[Swns>1]=1
        Swn1s[Swn1s<0]=0
        Swn1s[Swn1s>1]=1
        Adjs=self.adjs
        ##########teste########
        dp=p[Adjs[:,1]]-p[Adjs[:,0]]
        id_up=np.zeros_like(Adjs[:,0]).astype(int)
        id_up[dp>0]=Adjs[:,1][dp>0]
        id_up[dp<=0]=Adjs[:,0][dp<=0]
        sat_up1=Swn1s[id_up].copy()
        lrw1=sat_up1**2/0.01
        lro1=(1-sat_up1)**2/0.01
        lambda_ij=lrw1+lro1
        #################end teste###
        Ts=self.Ts*lambda_ij

        ID_vol=self.GID_0
        n=len(ID_vol)
        count=0
        #####only for alpha###
        la=np.concatenate([Adjs[:,0],Adjs[:,1],Adjs[:,0],Adjs[:,1]])
        ca=np.concatenate([Adjs[:,1],Adjs[:,0],Adjs[:,0],Adjs[:,1]])
        da=np.concatenate([Ts,Ts,-Ts,-Ts])
        self.Ta=sp.csc_matrix((da,(la,ca)),shape=(n,n))
        #####
        lines=[]
        cols=[]
        data=[]
        lines.append(ID_vol)
        cols.append(n+ID_vol)
        # import pdb; pdb.set_trace()
        data.append(self.F_Jacobian.c_o(0.3,np.repeat(time_step,n)))
        # J[ID_vol][n+ID_vol]+=float(F_Jacobian().c_o.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt}))
        lines.append(n+ID_vol)
        cols.append(n+ID_vol)
        data.append(self.F_Jacobian.c_w(0.3,np.repeat(time_step,n)))
        # J[n+ID_vol][n+ID_vol]+=float(F_Jacobian().c_w.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt}))
        linesq=[]
        dataq=[]
        linesq.append(ID_vol)
        dataq.append(self.F_Jacobian.acum_o(0.3,time_step,Swns,Swn1s))
        # q[ID_vol]+=float(F_Jacobian().acum_o.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt, Sw:Swns[count], Swn:Swn1s[count]}))
        linesq.append(n+ID_vol)
        dataq.append(self.F_Jacobian.acum_w(0.3,time_step,Swns,Swn1s))
        # q[n+ID_vol]+=float(F_Jacobian().acum_w.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt, Sw:Swns[count], Swn:Swn1s[count]}))
        # Adjs=np.array(self.adjs)
        adj0=np.array(Adjs[:,0])
        adj1=np.array(Adjs[:,1])
        ids0=ID_vol[adj0]
        ids1=ID_vol[adj1]
        ID_vol=ids0
        id_j=ids1
        swns0=Swns[ids0]
        swns1=Swns[ids1]
        pf0=p[adj0]
        pf1=p[adj1]
        # pf0=press0
        # pf1=press1
        # import pdb; pdb.set_trace()

        fi_o=(pf1-pf0)-self.F_Jacobian.r_o*self.gh
        fi_w=(pf1-pf0)-self.F_Jacobian.r_w*self.gh
        z0=self.centroids[self.adjs[:,0]][:,np.array(self.F_Jacobian.g)!=0].T[0]
        z1=self.centroids[self.adjs[:,1]][:,np.array(self.F_Jacobian.g)!=0].T[0]
        dz=-(z1-z0)
        nfi=len(Adjs)
        # fi_o[fi_w<0]=-fi_o[fi_w<0]
        # import pdb; pdb.set_trace()
        up0=(fi_o<1e-10)#|(swns0>swns1)
        up1=up0==False
        swf=np.zeros(nfi)
        swf[up0]=swns0[up0]
        swf[up1]=swns1[up1]
        id_up=np.zeros(nfi,dtype=np.int32)
        id_up[up0]=ids0[up0]
        id_up[up1]=ids1[up1]

        up0_w=(fi_w<1e-10)#|(swns0>swns1)
        up1_w=up0_w==False
        swf_w=np.zeros(nfi)
        swf_w[up0_w]=swns0[up0_w]
        swf_w[up1_w]=swns1[up1_w]
        id_up_w=np.zeros(nfi,dtype=np.int32)
        id_up_w[up0_w]=ids0[up0_w]
        id_up_w[up1_w]=ids1[up1_w]
        # import pdb; pdb.set_trace()
        gh=-self.gh
        # if self.iteration>0:
        #     gh=np.zeros_like(gh)

        J00=self.F_Jacobian.J[0][0](Ts,swf)
        # J00=float(self.F_Jacobian[0][0].subs({T:1, Sw:swf}))
        J01=self.F_Jacobian.J[0][1](Ts,swf, pf0, pf1,gh)
        # J01=float(self.F_Jacobian[0][1].subs({T:1, Sw:swf, p_i:pv, p_j:pj}))
        J10=self.F_Jacobian.J[1][0](Ts,swf_w)
        # J10=float(self.F_Jacobian[1][0].subs({T:1, Sw:swf}))
        J11=self.F_Jacobian.J[1][1](Ts,swf_w, pf0, pf1,gh)
        # J11=float(self.F_Jacobian[1][1].subs({T:1, Sw:swf, p_i:pv, p_j:pj}))
        linesq.append(ID_vol)
        dataq.append(-self.F_Jacobian.F_o(Ts,swf, pf0, pf1,gh))
        linesq.append(id_j)
        dataq.append(-self.F_Jacobian.F_o(Ts,swf, pf1, pf0,-gh))
        # q[ID_vol]-=float(F_Jacobian().F_o.subs({T:1.0, Sw:Swns1[count_fac], p_i:pv, p_j:pj}))
        linesq.append(n+ID_vol)
        dataq.append(-self.F_Jacobian.F_w(Ts,swf_w, pf0, pf1,gh))
        linesq.append(n+id_j)
        dataq.append(-self.F_Jacobian.F_w(Ts,swf_w, pf1, pf0,-gh))
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
        '''just for timing'''
        t0=time.time()
        l1=np.concatenate([ID_vol,id_j,ID_vol,id_j])
        c1=np.concatenate([ID_vol,id_j,id_j,ID_vol])
        d1=np.concatenate([-J00,-J00,J00,J00])
        self.Jpp=sp.csr_matrix((d1,(l1,c1)),shape=(n,n))
        self.time_Jpp=time.time()-t0
        ''' end just for timing'''
        # print(len(self.Jpp.data),'------------------')
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
        J=sp.csc_matrix((data,(lines,cols)),shape=(2*n,2*n))
        return(J, q)


    def apply_BC(self, lines, cols, data, q):
            n=int(len(q)/2)
            q[self.wells['ws_p']]=0
            q[self.wells['ws_inj']+n]=0
            for l in self.wells['ws_p']:
                data[lines==l]=0
                lines=np.append(lines,l)
                cols=np.append(cols,l)
                data=np.append(data,1)
            for l in np.setdiff1d(self.wells['ws_inj'],self.wells['ws_q']):
                data[lines==l+n]=0
                lines=np.append(lines,l+n)
                cols=np.append(cols,l+n)
                data=np.append(data,1)

            return lines, cols, data, q
