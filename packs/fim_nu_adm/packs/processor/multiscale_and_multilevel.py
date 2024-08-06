import numpy as np
import scipy.sparse as sp
import time
import os
from ..preprocessor.multiscale import get_dual_and_primal_1, get_local_problems_structure
from .assembler import Assembler
from .. import inputs
import pyamg as amg
from scipy.sparse.linalg import spilu
# from packs.postprocessor.exporter import FieldVisualizer
# visualize=FieldVisualizer()

class NewtonIterationMultilevel:
    def __init__(self, wells, faces, volumes, case):
        self.PVI=0
        self.alpha_lim=inputs.multiscale_and_multilevel_inputs['multilevel']['alpha_lim']
        self.beta_lim=inputs.multiscale_and_multilevel_inputs['multilevel']['beta_lim']
        self.alpha_beta_lim=np.infty
        self.GID_0=volumes['GID_0']
        self.wells=wells
        self.swns=np.zeros(len(self.GID_0))
        self.adjs=faces['adjacent_volumes']
        self.prep_time=[]
        t0=time.time()
        self.GID_1, self.DUAL_1 = get_dual_and_primal_1(volumes['centroids'])
        self.prep_time.append(time.time()-t0)#prep1
        t0=time.time()
        # se for usar local op descomentar
        # self.local_problems_structure, self.local_ID = get_local_problems_structure(self.DUAL_1, self.GID_1, faces['adjacent_volumes'],faces['permeabilities'])
        self.OP, self.OP_matrix = self.get_prolongation_operator(faces['permeabilities'])
        self.prep_time.append(time.time()-t0)#prep2
        self.OR_matrix = sp.csc_matrix((np.ones_like(self.GID_0), (self.GID_1, self.GID_0)), shape=(int(self.GID_1.max()+1),int(self.GID_0.max()+1)))

        self.get_op_to_alpha()
        self.get_beta_groups()
        self.Assembler = Assembler(wells, faces, volumes)
        self.proc_cumulative=[]
        self.int_prim_flag=[]
        self.PVI=[]
        self.porosities=volumes['pore_volume']
        # import pdb; pdb.set_trace()
        self.viz=np.load('results/'+case+'/viz.npy')

    def dual_aglomerator(self):
        t0=time.time()
        JP=self.Assembler.Jpp*self.OP_matrix
        RJP=self.OR_matrix*JP
        self.prep_time.append(time.time()-t0)#prep3
        l, c, d = sp.find(RJP)
        off_diags=l!=c
        l, c, d=l[off_diags], c[off_diags], d[off_diags]
        diags=RJP.diagonal()
        d=d/diags[l] #turns dd into neta
        upper=l>c
        d[upper][d[upper]<d[~upper]]=d[~upper][d[upper]<d[~upper]]
        l1, c1, d1 = l[upper], c[upper], d[upper]
        IJ=np.vstack([l1,c1])
        neta_IJ=np.vstack([d1,np.zeros_like(d1)]).max(axis=0)

    def get_operators(self):
        self.get_finescale_vols()
        t0=time.time()
        self.update_NU_ADM_mesh()
        self.update_averager()
        self.proc_temp.append(time.time()-t0) # time3
        t0=time.time()
        self.update_NU_ADM_operators()
        self.proc_temp.append(time.time()-t0) # time4
        self.update_R_and_P()

        rar=self.P.toarray()
        # visualize.plot_labels(self.NU_ADM_ID)
        # visualize.plot_labels(self.GID_0)
        # visualize.plot_labels(rar[:,5][:-36])
        return self.R, self.P

    def get_finescale_vols(self):
        flag=-np.ones_like(self.GID_1)
        swns=self.swns
        ws_p=self.wells['ws_prod']
        ws_i=self.wells['ws_inj']
        adjs=self.adjs
        deltas=abs(swns[adjs][:,0]-swns[adjs][:,1])
        same=(self.GID_1[adjs[:,0]]==self.GID_1[adjs[:,1]])

        # fs=np.arange(len(deltas))[(deltas>0.05)&(same|(swns.sum(axis=0)<.3))]
        # fs=np.arange(len(deltas))[(deltas>0.05)&((swns.sum(axis=0)<.3)|(deltas>.1))]
        fs=np.arange(len(deltas))[(deltas>0.1) & ((swns[adjs][:,0]<0.8) & (swns[adjs][:,1]<0.8))]
        # import pdb; pdb.set_trace()
        # fs=np.arange(len(deltas))[((deltas>0.05))&((swns.sum(axis=0)>0.2)|same)]
        vols=adjs[fs].flatten()
        # vols=vols[swns[vols]<0.3]
        # for i in range(1):
        #     flag[vols]=1
        #     viz=np.unique(adjs[flag[adjs].sum(axis=1)==0])
        #     vols=np.unique(np.concatenate([vols,viz]))

        t0=time.time()
        self.update_alpha()
        self.proc_temp.append((time.time()-t0)/10)#time1
        self.alpha_beta=self.alphas*self.betas
        alpha_vols=self.GID_0[(self.alphas>self.alpha_lim)]# | (self.alpha_beta>self.alpha_beta_lim)]
        # import pdb; pdb.set_trace()
        '''
        # c_ws_p=np.concatenate([self.GID_0[self.GID_1==p] for p in np.unique(self.GID_1[ws_p])])
        # import pdb; pdb.set_trace()
        fs_vs=np.unique(np.concatenate([c_ws_p, ws_i, vols, alpha_vols]))'''
        fs_vs=np.unique(np.concatenate([ws_i, vols, alpha_vols]))
        # import pdb; pdb.set_trace()
        bs=self.beta_ind[fs_vs]
        binds=np.unique(bs[bs>-1])

        if len(binds)>0:
            bvols=np.concatenate(self.beta_groups[binds])
            fs_vs=np.unique(np.concatenate([fs_vs,bvols]))
        self.fs_vols=fs_vs.copy()
        # import pdb; pdb.set_trace()

    def get_finescale_vols_haji(self):
        flag=-np.ones_like(self.GID_1)
        swns=self.swns
        ws_p=self.wells['ws_prod']
        ws_i=self.wells['ws_inj']
        adjs=self.adjs
        deltas=abs(swns[adjs][:,0]-swns[adjs][:,1])
        fs=np.arange(len(deltas))[deltas>0.1]


    def update_alpha_stabls(self):
        np.set_printoptions(5)
        nf=int(len(self.q)/2)
        # import pdb; pdb.set_trace()
        Jpp=self.J[0:nf,0:nf]

        JP=Jpp*self.OP_matrix
        RJP=self.OR_matrix*JP
        self.GID_0
        self.GID_1
        # JP=JP.tocsr().tocsc()
        l, c, d=sp.find(JP)
        same=c==self.GID_1[l]
        diag=RJP.diagonal()[self.GID_1]
        lines=l[~same]
        maxs=np.zeros(nf)
        np.maximum.at(maxs,lines,d[~same])
        # import pdb; pdb.set_trace()
        self.alphas=maxs/abs(diag)
    # @ profile
    def get_op_to_alpha(self):
        self.OP_to_alpha=self.OP_matrix.copy()
        l, c, _ =sp.find(self.OP_matrix)
        ad1s=self.GID_1[self.adjs]
        self.bound_prim=np.repeat(False,len(self.GID_1))
        self.bound_prim[ad1s[ad1s[:,0]!=ad1s[:,1]]]=True
        self.int_prim_flag=c==self.GID_1[l] | self.bound_prim[l]
        self.OP_to_alpha.data[self.int_prim_flag]=0

    def update_alpha(self):
        np.set_printoptions(5)
        nf=int(self.OP_matrix.shape[0]/2)

        JP=self.Assembler.Ta*self.OP_to_alpha
        RJP=self.OR_matrix*JP

        l, c, d=sp.find(JP)
        same=c==self.GID_1[l]
        diag=RJP.diagonal()[self.GID_1]
        lines=l[~same]
        maxs=np.zeros(l.max()+1)
        # import pdb; pdb.set_trace()
        np.maximum.at(maxs,lines,d[~same])
        # import pdb; pdb.set_trace()
        # da=abs(self.alphas-maxs/abs(diag)).max()
        self.alphas=maxs/abs(diag)

    def get_beta_groups(self):
        pos=self.GID_1[self.OP[0]]==self.OP[1]
        phis=self.OP[2][pos][np.argsort(self.OP[0][pos])]
        self.betas=(1-phis)/phis
        beta_facs=self.betas[self.adjs].max(axis=1)
        ads=self.adjs[beta_facs>self.beta_lim]
        map=np.arange(self.adjs.max())
        uads=np.unique(ads)
        n=len(uads)
        map[uads]=np.arange(n)
        adjs=map[ads]
        adjs=np.vstack([adjs,np.array([adjs[:,1],adjs[:,0]]).T])
        graph = sp.csc_matrix((np.ones_like(adjs[:,0]), (adjs[:,0], adjs[:,1])),shape=(n,n))
        n_l,labels=sp.csgraph.connected_components(graph)
        self.beta_ind=-np.ones_like(self.GID_0)
        self.beta_ind[uads]=labels
        self.beta_groups=np.array([uads[labels==l] for l in range(n_l)])

    def get_prolongation_operator_fast(self, ts):
        i=-1
        ops = []
        glines=[]
        gcols=[]
        gdata=[]
        prob=0
        if len(self.local_problems_structure[-1][0][0])==0:
            self.local_problems_structure=self.local_problems_structure[:-1]
        for structure in self.local_problems_structure:
            ops.append([])
            i+=1
            for local_problem in structure:
                internal_matrix = local_problem[0][0]
                off_diagonal_entries = local_problem[0][1]
                diagonal_entries = local_problem[0][2]
                acumulator = local_problem[0][3]
                internal_gids = local_problem[0][4]
                new_data = np.concatenate([ts[off_diagonal_entries], -ts[diagonal_entries]])
                sums=np.bincount(acumulator,weights=new_data)
                internal_matrix=internal_matrix.tocsc()
                # internal_matrix.data=np.arange(len(internal_matrix.data))
                # if prob==1:

                internal_matrix.data=sums[sums!=0]

                # import pdb; pdb.set_trace()


                for local_external_problem in local_problem[1]:
                    external_matrix = local_external_problem[0]
                    entries = local_external_problem[1]
                    external_gids = local_external_problem[2]
                    entity_up_ids = local_external_problem[3]
                    matrix_connection = local_external_problem[4]
                    columns= local_external_problem[5]
                    if entity_up_ids.max()>-1:
                        d=[]
                        for e in entity_up_ids:
                            d.append(ops[i-1][e])
                        d=np.concatenate(d)
                        # import pdb; pdb.set_trace()
                        matrix_connection.data=d[matrix_connection.data-1]
                        external_matrix.data = ts[entries]
                        external_matrix = external_matrix*matrix_connection
                    else:
                        external_matrix.data = ts[entries]
                        entity_up_ids=external_gids
                op=-sp.linalg.spsolve(internal_matrix, external_matrix)
                # if ((op.sum(axis=1)<0.999) | (op.sum(axis=1)>1.001)).sum()!=0:
                #     # import pdb; pdb.set_trace()
                #     prob=1

                # import pdb; pdb.set_trace()
                fop=sp.find(op)
                glines.append(internal_gids[fop[0]])
                gcols.append(columns[fop[1]])
                gdata.append(fop[2])
                data=op.data
                ops[i].append(data)

        all_volumes=np.arange(len(self.DUAL_1))
        vertices=self.DUAL_1==3
        glines.append(all_volumes[vertices])
        gcols.append(all_volumes[vertices])
        gdata.append(np.ones_like(all_volumes[vertices]))

        glines=np.concatenate(glines)
        gcols=np.concatenate(gcols)
        mapg=-np.ones(gcols.max()+1).astype(int)
        mapg[np.unique(gcols)]=np.arange(len(np.unique(gcols)))
        gcols=mapg[gcols]
        gdata=np.concatenate(gdata)
        op1=[glines, gcols, gdata]
        OP_AMS=sp.csc_matrix((gdata, (glines, gcols)),shape=(int(glines.max()+1), int(gcols.max())+1))
        return op1, OP_AMS

    def get_prolongation_operator(self,ts):
        # _,_1=self.get_prolongation_operator_fast(ts )
        # import pdb; pdb.set_trace()
        dual_1=self.DUAL_1
        n=len(dual_1)
        nv=(dual_1==3).sum()
        ne=(dual_1==2).sum()
        nf=(dual_1==1).sum()
        ni=(dual_1==0).sum()
        wire=-np.ones_like(dual_1)
        wire[dual_1==0]=range(ni)
        wire[dual_1==1]=ni+range(nf)
        wire[dual_1==2]=ni+nf+range(ne)
        wire[dual_1==3]=ni+nf+ne+range(nv)
        # G = sp.csc_matrix((np.ones_like(dual_1), (wire, self.GID_0)),shape=(n,n))
        a=wire[self.adjs.T]

        lines=np.concatenate([a[0],a[1],a[0],a[1]])
        cols=np.concatenate([a[1],a[0],a[0],a[1]])
        data=np.concatenate([ts,ts,-ts,-ts])
        T = sp.csc_matrix((data, (lines, cols)),shape=(n,n))
        # W=G*T*G.T
        Tii=T[0:ni,0:ni]
        Tif=T[0:ni,ni:ni+nf]
        Tff=T[ni:ni+nf,ni:ni+nf]
        Tfe=T[ni:ni+nf,ni+nf:ni+nf+ne]
        Tee=T[ni+nf:ni+ne+nf,ni+nf:ni+ne+nf]
        Tev=T[ni+nf:ni+ne+nf,ni+ne+nf:n]
        Tvv=T[ni+ne+nf:n,ni+ne+nf:n]
        Tvv.data=np.ones(nv)

        auxf=sp.csc_matrix((np.array(Tif.sum(axis=0))[0], (range(nf), range(nf))),shape=(nf,nf))
        auxe=sp.csc_matrix((np.array(Tfe.sum(axis=0))[0], (range(ne), range(ne))),shape=(ne,ne))
        Tee+=auxe
        Tff+=auxf

        ope=sp.linalg.spsolve(Tee,-Tev*Tvv)
        opf=sp.linalg.spsolve(Tff,-Tfe*ope)
        opi=sp.linalg.spsolve(Tii,-Tif*opf)
        prol=sp.vstack([opi,opf,ope,Tvv])


        G=sp.csc_matrix((np.ones(n), (wire.astype(int), self.GID_0)),shape=(n, n))
        fp=sp.find(G.T*prol)
        lp=fp[0]
        cp=fp[1]
        dp=fp[2]
        op1=[lp, cp, dp]

        OP_AMS=sp.csc_matrix((dp, (lp, cp)),shape=(int(lp.max()+1), int(cp.max())+1))

        return op1, OP_AMS

    def update_averager(self):
        # self.fs_vols=fs_vols
        levels=np.ones_like(self.GID_1)
        NU_ADM_ID = -self.levels
        levels[self.fs_vols]=0
        gid1_adjs=self.GID_1[self.adjs]
        same_gid=gid1_adjs[:,0]==gid1_adjs[:,1]
        cc_adjs=self.levels[self.adjs].sum(axis=1)==2
        adjs=self.adjs[same_gid & cc_adjs]
        fines=np.tile(self.fs_vols,(2,1)).T
        adjs=np.vstack([fines,adjs])
        adjs=np.tile(adjs,(2,1))
        data = np.ones(len(adjs))
        n=len(self.levels)
        graph = sp.csc_matrix((data, (adjs[:,0], adjs[:,1])),shape=(n,n))
        n,labels=sp.csgraph.connected_components(graph)
        '''
        self.NU_ADM_ID=labels
        gid1=self.GID_1[self.GID_0[self.DUAL_1==3]]
        self.coarse_id_NU_ADM=gid1
        '''
        cols=self.GID_0
        lines=labels
        data=np.ones_like(cols)
        # averager=sp.csc_matrix(())
        # import pdb; pdb.set_trace()
        self.averager=sp.csc_matrix((data, (lines, cols)), shape=(lines.max()+1, cols.max()+1))

    def update_NU_ADM_mesh(self):
        # self.fs_vols=fs_vols
        self.levels=np.ones_like(self.GID_1)
        self.NU_ADM_ID = -self.levels
        self.levels[self.fs_vols]=0
        coarse_volumes =  self.levels==1
        self.NU_ADM_ID[coarse_volumes]=self.GID_1[coarse_volumes]
        all_cvs=np.unique(self.NU_ADM_ID)
        if all_cvs.min()==-1:
            all_cvs=all_cvs[1:]
        remaining_ids=np.setdiff1d(np.unique(self.GID_1), all_cvs)
        nids=len(self.fs_vols)-len(remaining_ids)
        ids=np.concatenate([remaining_ids, self.GID_1.max()+np.arange(nids)+1])
        # import pdb; pdb.set_trace()
        self.NU_ADM_ID[self.fs_vols]=ids

        self.vertices=self.GID_0[self.DUAL_1==3]
        gid1=self.GID_1[self.vertices]
        for rgid in remaining_ids:
            if rgid in gid1:
                gid1[gid1==rgid]=self.NU_ADM_ID[self.vertices[self.GID_1[self.vertices]==rgid]]
        self.coarse_id_NU_ADM=gid1

    def update_NU_ADM_operators(self):
        l, c, d=self.OP
        coarse=self.levels[l]==1
        # import pdb; pdb.set_trace()
        # mapc = self.NU_ADM_ID[self.DUAL_1==3] #aqui trocar por linha abaixo
        mapc = self.coarse_id_NU_ADM
        lines = l[coarse]
        cols = mapc[c[coarse]]
        # import pdb; pdb.set_trace()
        same=self.GID_1[lines]==c[coarse]
        cols[same]=self.NU_ADM_ID[lines[same]]
        # import pdb; pdb.set_trace()
        # cols = self.NU_ADM_ID[lines]
        data = d[coarse]
        ls = self.fs_vols
        cs = self.NU_ADM_ID[self.fs_vols]
        ds = np.ones_like(cs)

        lines = np.concatenate([lines, ls])
        cols = np.concatenate([cols, cs])
        data = np.concatenate([data, ds])

        self.NU_ADM_OP = [lines, cols, data]
        # import pdb; pdb.set_trace()
        # visualize.plot_labels(self.OP[:,4].T.toarray()[0])
        # visualize.plot_labels(self.NU_ADM_OP[:,12].T.toarray()[0])
        # visualize.plot_labels(self.levels)
        # import pdb; pdb.set_trace()
        cols = self.GID_0
        lines = self.NU_ADM_ID
        data = np.ones(len(lines))
        self.NU_ADM_OR = [lines,cols, data]

    def update_R_and_P(self):
        lp, cp, dp = self.NU_ADM_OP
        lr, cr, dr = self.NU_ADM_OR
        n_f, n_ADM=lp.max()+1, cp.max()+1
        lP=np.concatenate([lp, cr+n_f])
        cP=np.concatenate([cp, lr+n_ADM])
        dP=np.concatenate([dp, dr])

        lR=np.concatenate([lr, lr+n_ADM])
        cR=np.concatenate([cr, cr+n_f])
        dR=np.concatenate([dr, dr])
        self.R=sp.csc_matrix((dR, (lR, cR)), shape=(2*n_ADM, 2*n_f))
        self.P=sp.csc_matrix((dP, (lP, cP)), shape=(2*n_f, 2*n_ADM))

    def tams(self, x0, R, P, tol=1e-4):
        # pv=P*sp.linalg.spsolve(P.T*self.J*P,-P.T*self.q)
        pv=x0
        r0=np.linalg.norm(self.q)
        r1=r0.copy()
        ct=0
        while r1/r0>tol:
            pv12=pv-P*sp.linalg.spsolve(P.T*self.J*P, P.T*(self.q+self.J*pv))
            # t1=time.time()
            pv=pv12-sp.linalg.lgmres(self.J, self.q+self.J*pv12,x0=(pv-pv12)/2,tol=1e-10, maxiter=15)[0]
            # t2=time.time()
            # pv=pv12-amg.krylov.gmres(self.J, self.q+self.J*pv12,restart=5,maxiter=15)[0]
            # pv=pv12-iluJ.solve(self.q+self.J*pv12)
            # import pdb; pdb.set_trace()

            # sa=-sp.linalg.spsolve(self.J,self.q)
            # pvt=self.jacobi_iteration(sa,5)
            # import pdb; pdb.set_trace()

            # t3=time.time()

            # import pdb; pdb.set_trace()
            print(np.linalg.norm(self.q+self.J*pv)/r0,np.linalg.norm(self.q+self.J*pv12)/r0,"iterative")
            r1=np.linalg.norm(self.q+self.J*pv)
        sol=pv.copy()
        return sol

    def jacobi_iteration(self, x0, ni):
        Diag=self.J.diagonal()
        l=range(len(x0))
        D=sp.csc_matrix((Diag, (l,l)),shape=(len(x0),len(x0)))
        LU=self.J-D
        D_1=sp.csc_matrix((1/Diag, (l,l)),shape=(len(x0),len(x0)))
        x1=x0.copy()
        for i in range(ni):
            x1=D_1*((self.q+self.J*x1)+LU*x1)
            import pdb; pdb.set_trace()
        return x1

    def cpr_nu_adm(self,sol0,R,P):
        J=self.J
        # r=self.q+self.J*sol0
        nc,nf=R.shape
        nc,nf=int(nc/2),int(nf/2)
        # rp=self.r[0:nf]
        # rs=self.r[nf:]
        Jpp=J[0:nf,0:nf]
        Jps=J[0:nf,nf:]
        Jsp=J[nf:,0:nf]
        Jss=J[nf:,nf:]
        R=R[0:nc,0:nf]
        P=P[0:nf,0:nc]
        lc=range(nf)
        Dss=sp.csc_matrix((Jss.diagonal(),(lc,lc)),shape=(nf,nf))
        Dps=sp.csc_matrix((Jps.diagonal(),(lc,lc)),shape=(nf,nf))
        Jssd=Jss.diagonal()
        # Jpsd[Jpsd==0]=1e-5
        # Dpp_1=sp.csc_matrix((1/Jpp.diagonal(),(lc,lc)),shape=(nf,nf))
        Dss_1=sp.csc_matrix((1/Jssd,(lc,lc)),shape=(nf,nf))
        # Vd=1/(Jpp.diagonal()*Jss.diagonal()-Jps.diagonal()*Jsp.diagonal())
        # V_1=sp.csc_matrix((Vd,(lc,lc)),shape=(nf,nf))
        xv=sol0.copy()
        rv=self.q+self.J*sol0
        Jppx=(Jpp-Dps*Dss_1*Jsp) #matriz cpr
        # rx=rv[0:nf]-Dps*Dss_1*rv[nf:]
        # dpx=sp.linalg.spsolve(Jppx,rx) #aplicar tams
        #Modifs
        sfim=-sp.linalg.spsolve(self.J,self.q)
        sf=sfim[0:nf]


        rxvp=rv[0:nf]-Dps*Dss_1*rv[nf:]
        dxp_vs1=sol0[0:nf]-sp.linalg.lgmres(Jppx, rxvp,x0=np.zeros_like(sol0[0:nf]),tol=1e-10, maxiter=15)[0]

        sol1=dxp_vs1
        sol1=np.concatenate([sol1,sol0[nf:]])
        r1=self.q+self.J*sol1
        rp_vs1=r1[0:nf]-Dps*Dss_1*r1[nf:]
        dxp_vs1_1=dxp_vs1-P*sp.linalg.spsolve(R*Jppx*P,R*rp_vs1)
        rp_vs1_1=rp_vs1+Jppx*dxp_vs1_1

        sol2=dxp_vs1_1
        sol2=np.concatenate([sol2,sol0[nf:]])
        r2=self.q+self.J*sol2
        rp_vs1_1=r2[0:nf]-Dps*Dss_1*r2[nf:]
        dxp_vs1_1_s2=dxp_vs1_1-sp.linalg.lgmres(Jppx, rp_vs1_1,x0=np.zeros_like(dxp_vs1_1),tol=1e-10, maxiter=15)[0]
        rp_vs1_1_s2=self.q[0:nf]+Jppx*dxp_vs1_1_s2

        rv_1=rv+(self.J[:,0:nf]*dxp_vs1_1_s2)
        xv_1=xv.copy()
        xv_1[0:nf]+=dxp_vs1_1_s2
        # ilu = sp.linalg.spilu(self.J)
        sA_iLU = sp.linalg.spilu(self.J)
        M = sp.linalg.LinearOperator((nf*2,nf*2), sA_iLU.solve)
        x = sp.linalg.gmres(self.J,rv_1,M=M,maxiter=15)[0]
        xv_1=xv_1-x

        print(np.linalg.norm(sol0[0:nf]-sf)/np.linalg.norm(sf),np.linalg.norm(dxp_vs1-sf)/np.linalg.norm(sf),np.linalg.norm(dxp_vs1_1-sf)/np.linalg.norm(sf),np.linalg.norm(dxp_vs1_1_s2-sf)/np.linalg.norm(sf),np.linalg.norm(xv_1[0:nf]-sf)/np.linalg.norm(sf))

        # import pdb; pdb.set_trace()
        # xv_1=xv_1-sp.linalg.spsolve(self.J, rv_1)
        # xv_1=xv_1-sp.linalg.lgmres(self.J, rv_1,x0=np.zeros_like(xv),tol=1e-10, maxiter=30)[0]

        # print(np.linalg.norm(xv_1-sfim)/np.linalg.norm(sfim))
        # import pdb; pdb.set_trace()
        return xv_1

    def cpr_nu_adm_teste(self,sol0,R,P):
        J=self.J
        # r=self.q+self.J*sol0
        nc,nf=R.shape
        nc,nf=int(nc/2),int(nf/2)
        # rp=self.r[0:nf]
        # rs=self.r[nf:]
        Jpp=J[0:nf,0:nf]
        Jps=J[0:nf,nf:]
        Jsp=J[nf:,0:nf]
        Jss=J[nf:,nf:]
        R=R[0:nc,0:nf]
        P=P[0:nf,0:nc]
        lc=range(nf)
        Dss=sp.csc_matrix((Jss.diagonal(),(lc,lc)),shape=(nf,nf))
        Dps=sp.csc_matrix((Jps.diagonal(),(lc,lc)),shape=(nf,nf))
        Jssd=Jss.diagonal()
        # Jpsd[Jpsd==0]=1e-5
        # Dpp_1=sp.csc_matrix((1/Jpp.diagonal(),(lc,lc)),shape=(nf,nf))
        Dss_1=sp.csc_matrix((1/Jssd,(lc,lc)),shape=(nf,nf))
        # Vd=1/(Jpp.diagonal()*Jss.diagonal()-Jps.diagonal()*Jsp.diagonal())
        # V_1=sp.csc_matrix((Vd,(lc,lc)),shape=(nf,nf))
        xv=sol0.copy()
        # rv=self.q+self.J*sol0
        Jppx=(Jpp-Dps*Dss_1*Jsp) #matriz cpr
        # rx=rv[0:nf]-Dps*Dss_1*rv[nf:]
        # dpx=sp.linalg.spsolve(Jppx,rx) #aplicar tams
        #Modifs
        sfim=-sp.linalg.spsolve(self.J,self.q)
        sf=sfim[0:nf]


        rxvp=rv[0:nf]-Dps*Dss_1*rv[nf:]
        rxv=np.concatenate([rxvp,self.q[nf:]])
        # dxp_vs1=sol0[0:nf]-sp.linalg.gmres(Jppx, rxvp,x0=sol0[0:nf],tol=1e-10, maxiter=150)[0]
        #
        # sol1=dxp_vs1
        # sol1=np.concatenate([sol1,sol0[nf:]])
        # r1=self.q+self.J*sol1
        # rp_vs1=r1[0:nf]-Dps*Dss_1*r1[nf:]
        dxp_vs1_1=-P*sp.linalg.spsolve(R*Jppx*P,R*rxvp)
        rp_vs1_1=rp_vs1+Jppx*dxp_vs1_1

        sol2=dxp_vs1_1
        sol2=np.concatenate([sol2,sol0[nf:]])
        r2=self.q+self.J*sol2
        rp_vs1_1=r2[0:nf]-Dps*Dss_1*r2[nf:]
        dxp_vs1_1_s2=dxp_vs1_1-sp.linalg.lgmres(Jppx, rp_vs1_1,x0=np.zeros_like(dxp_vs1_1),tol=1e-10, maxiter=15)[0]
        rp_vs1_1_s2=self.q[0:nf]+Jppx*dxp_vs1_1_s2

        rv_1=rv+(self.J[:,0:nf]*dxp_vs1_1_s2)
        xv_1=xv.copy()
        xv_1[0:nf]+=dxp_vs1_1_s2
        # ilu = sp.linalg.spilu(self.J)
        sA_iLU = sp.linalg.spilu(self.J)
        M = sp.linalg.LinearOperator((nf*2,nf*2), sA_iLU.solve)
        x = sp.linalg.gmres(self.J,rv_1,M=M,maxiter=15)[0]
        xv_1=xv_1-x

        print(np.linalg.norm(sol0[0:nf]-sf)/np.linalg.norm(sf),np.linalg.norm(dxp_vs1-sf)/np.linalg.norm(sf),np.linalg.norm(dxp_vs1_1-sf)/np.linalg.norm(sf),np.linalg.norm(dxp_vs1_1_s2-sf)/np.linalg.norm(sf),np.linalg.norm(xv_1[0:nf]-sf)/np.linalg.norm(sf))

        import pdb; pdb.set_trace()
        # xv_1=xv_1-sp.linalg.spsolve(self.J, rv_1)
        # xv_1=xv_1-sp.linalg.lgmres(self.J, rv_1,x0=np.zeros_like(xv),tol=1e-10, maxiter=30)[0]

        # print(np.linalg.norm(xv_1-sfim)/np.linalg.norm(sfim))
        # import pdb; pdb.set_trace()
        return xv_1

    def tams_CPR(self,x0,Jx,rx, R, P, tol=1e-4):
        # pv=P*sp.linalg.spsolve(P.T*self.J*P,-P.T*self.q)
        # import pdb; pdb.set_trace()
        self.r=self.q+self.J*x0
        nc,nf=R.shape
        nc,nf=int(nc/2),int(nf/2)
        R=R[0:nc,0:nf]
        P=P[0:nf,0:nc]
        # # Modifs
        # xp=x0[0:nf]
        # dxp=-sp.linalg.lgmres(Jx, rx+Jx*pv12,x0=xp,tol=1e-10, maxiter=15)[0]
        # # dxp_vs1=xp+dxp
        # rp_vs1=(self.q+self.J*(np.concatenate([xp+dxp,self.r[nf:]]))[0:nf]
        # dxp_c=-P*sp.linalg.spsolve(P.T*Jx*P, P.T*(rp_vs1))


        #Fim modifs
        # import pdb; pdb.set_trace()
        pv=x0.copy()
        n=len(pv)
        r0=np.linalg.norm(self.q)
        r1=r0.copy()
        ct=0
        # sf=-sp.linalg.spsolve(self.J,self.q)
        while (r1/r0>tol) and (ct<30):
            pv12=pv-P*sp.linalg.spsolve(P.T*Jx*P, P.T*(rx+Jx*pv))
            pv=pv12-sp.linalg.lgmres(Jx, rx+Jx*pv12,x0=(pv-pv12)/10,tol=1e-10, maxiter=15)[0]
            pv=pv12-amg.krylov.bicgstab(Jx, rx+Jx*pv12,maxiter=15)[0]
            # print(np.linalg.norm(sf[0:nf]-pv)/np.linalg.norm(sf[0:nf]))
            print(np.linalg.norm(rx+Jx*pv)/r0,np.linalg.norm(rx+Jx*pv12)/r0,"iterative")
            ct+=1
            import pdb; pdb.set_trace()
        return pv
    # @profile
    def newton_iteration_ADM(self, p, s, time_step, cpr, rel_tol=1e-3):
        nf=len(p)
        pressure = p.copy()
        swns = s.copy()
        swn1s = s.copy()
        converged=False
        count=0
        dt=time_step
        # import pdb; pdb.set_trace()
        ep_l2=[]
        ep_li=[]
        es_l1=[]
        while not converged:
            self.proc_temp=[]
            # swns[self.Assembler.wells['ws_inj']]=1
            self.swns=swns.copy()
            self.J, self.q=self.Assembler.get_jacobian_matrix(swns, swn1s, pressure, time_step)
            n=int(len(self.q)/2)
            self.proc_temp.append(self.Assembler.time_Jpp) #prep_time1
            R, P = self.get_operators()
            t0=time.time()

            if cpr:
                sol=-P*sp.linalg.spsolve(R*self.J*P, P.T*self.q) # descomentar
                sol=self.cpr_nu_adm(sol,R,P)#comentar
            else:
                sol=-P*sp.linalg.spsolve(R*self.J*P, R*self.q) # descomentar

            # self.r=self.q+self.J*sol
            self.P_NU_ADM_orig=pressure.copy()
            sf=-sp.linalg.spsolve(self.J,self.q)
            ep_l2.append(np.linalg.norm(sf[0:nf]-sol[0:nf])/np.linalg.norm(sf[0:nf]))
            ep_li.append(abs(sf[0:nf]-sol[0:nf]).max()/abs(sf[0:nf]).max())
            es_l1.append(abs(sf[nf:]-sol[nf:]).sum()/nf)
            # import pdb; pdb.set_trace()
            # sol=-sp.linalg.spsolve(self.J,self.q)


            print(np.linalg.norm(sf[0:nf]-sol[0:nf])/np.linalg.norm(sf[0:nf]),np.linalg.norm(sf[nf:]-sol[nf:])/np.linalg.norm(sf[nf:]),"aqui")
            # import pdb; pdb.set_trace()
            # aa=self.cpr_nu_adm(sol.copy(),R,P)#descomentar

            # print(np.linalg.norm(sf-sol)/np.linalg.norm(sf),"CPR")
            dim=self.J.shape

            # fim do bacalho
            self.proc_temp.append(time.time()-t0) #time5
            self.proc_cumulative.append(self.proc_temp)
            pressure+=sol[0:n]
            # import pdb; pdb.set_trace()
            sol[n+self.wells['ws_prod']]=0
            sol[self.wells['ws_inj']]=0
            swns+=sol[n:]
            swns[self.Assembler.wells['ws_inj']]=1
            # self.cpr_nu_adm(R,P)
            converged=max(abs(sol))<rel_tol and (swns.max()<1.00000001 and swns.min()>-0.000001) #descomentar

            # if converged:
            #     sol=self.tams(sol,R,P)
            swns[swns>1]=1
            swns[swns<0]=0
            print(count,f'{ max(abs(sol)):.2e}',swns.sum()/len(swns))
            # import pdb; pdb.set_trace()
            count+=1
            # pressure[self.wells["ws_p"]]
            self.PVI=(self.swns*self.porosities).sum()/self.porosities.sum()
            if (count>7) or max(abs(sol[n:]))>3:
                print('excedded maximum number of iterations finescale')
                return False, count, pressure, swns
            swns[self.wells['ws_prod']]=swns[self.viz].sum()/len(self.viz)
            na, nf=int(self.R.shape[0]/2), int(self.R.shape[1]/2)
            if not cpr:
                OR_ADM=self.averager
                swns=OR_ADM.T*(OR_ADM.tocsr()*swns/np.array(OR_ADM.sum(axis=1)).T[0]) #ativar se for sem smoother
        # import pdb; pdb.set_trace()
        self.ep_l2=np.array(ep_l2).sum()/len(ep_l2)
        self.ep_li=np.array(ep_li).sum()/len(ep_li)
        self.es_l1=np.array(es_l1).sum()/len(es_l1)
        # import pdb; pdb.set_trace()
        return True, count, pressure, swns

    @staticmethod
    def define_NU_ADM_mesh(DUAL_1: np.ndarray, GID_0: np.ndarray, GID_1: np.ndarray, fs_vols: np.ndarray):
        levels = np.ones_like(GID_1)
        NU_ADM_ID = -levels
        levels[fs_vols] = 0
        coarse_volumes = levels == 1
        NU_ADM_ID[coarse_volumes]=GID_1[coarse_volumes]
        all_cvs=np.unique(NU_ADM_ID)
        if all_cvs.min()==-1:
            all_cvs=all_cvs[1:]
        remaining_ids=np.setdiff1d(np.unique(GID_1), all_cvs)
        nids=len(fs_vols)-len(remaining_ids)
        ids=np.concatenate([remaining_ids, GID_1.max()+np.arange(nids)+1])
        # import pdb; pdb.set_trace()
        NU_ADM_ID[fs_vols]=ids

        vertices=GID_0[DUAL_1==3]
        gid1=GID_1[vertices]
        for rgid in remaining_ids:
            if rgid in gid1:
                gid1[gid1==rgid]=NU_ADM_ID[vertices[GID_1[vertices]==rgid]]
        coarse_id_NU_ADM=gid1

        return coarse_id_NU_ADM, NU_ADM_ID
    
    @staticmethod
    def update_averager_v0(GID_0, GID_1, fs_vols, adjs):
        # self.fs_vols=fs_vols
        levels=np.ones_like(GID_1)
        NU_ADM_ID = -levels
        levels[fs_vols]=0
        gid1_adjs=GID_1[adjs]
        same_gid=gid1_adjs[:,0]==gid1_adjs[:,1]
        cc_adjs=levels[adjs].sum(axis=1)==2
        adjs2=adjs[same_gid & cc_adjs]
        fines=np.tile(fs_vols,(2,1)).T
        adjs2=np.vstack([fines,adjs2])
        adjs2=np.tile(adjs2,(2,1))
        data = np.ones(len(adjs2))
        n=len(levels)
        graph = sp.csc_matrix((data, (adjs2[:,0], adjs2[:,1])),shape=(n,n))
        n,labels=sp.csgraph.connected_components(graph)
        '''
        self.NU_ADM_ID=labels
        gid1=self.GID_1[self.GID_0[self.DUAL_1==3]]
        self.coarse_id_NU_ADM=gid1
        '''
        cols=GID_0
        lines=labels
        data=np.ones_like(cols)
        # averager=sp.csc_matrix(())
        # import pdb; pdb.set_trace()
        averager=sp.csc_matrix((data, (lines, cols)), shape=(lines.max()+1, cols.max()+1))
        return averager

    @staticmethod
    def update_NU_ADM_operators_v0(OP, levels, coarse_id_NU_ADM, GID_1, GID_0, NU_ADM_ID, fs_vols):
        l, c, d=OP
        coarse=levels[l]==1
        # import pdb; pdb.set_trace()
        # mapc = self.NU_ADM_ID[self.DUAL_1==3] #aqui trocar por linha abaixo
        mapc = coarse_id_NU_ADM
        lines = l[coarse]
        cols = mapc[c[coarse]]
        # import pdb; pdb.set_trace()
        same=GID_1[lines]==c[coarse]
        cols[same]=NU_ADM_ID[lines[same]]
        # import pdb; pdb.set_trace()
        # cols = self.NU_ADM_ID[lines]
        data = d[coarse]
        ls = fs_vols
        cs = NU_ADM_ID[fs_vols]
        ds = np.ones_like(cs)

        lines = np.concatenate([lines, ls])
        cols = np.concatenate([cols, cs])
        data = np.concatenate([data, ds])

        NU_ADM_OP = [lines, cols, data]
        # import pdb; pdb.set_trace()
        # visualize.plot_labels(self.OP[:,4].T.toarray()[0])
        # visualize.plot_labels(self.NU_ADM_OP[:,12].T.toarray()[0])
        # visualize.plot_labels(self.levels)
        # import pdb; pdb.set_trace()
        cols = GID_0
        lines = NU_ADM_ID
        data = np.ones(len(lines))
        NU_ADM_OR = [lines,cols, data]

        lp, cp, dp = NU_ADM_OP
        lr, cr, dr = NU_ADM_OR
        n_f, n_ADM=lp.max()+1, cp.max()+1

        OP_NU_ADM = sp.csc_matrix((dp, (lp, cp)), shape=(n_f, n_ADM))
        OR_NU_ADM = sp.csc_matrix((dr, (lr, cr)), shape=(n_ADM, n_f))

        return OP_NU_ADM, OR_NU_ADM, None