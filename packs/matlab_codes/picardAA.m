function [p_new,tabletol,iter,ciclos]=picardAA(M_old,RHS_old,nitpicard,tolpicard,kmap,...
    parameter,w,s,nflagface,fonte,p_old,gamma,nflagno,...
    weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface...
    ,gravresult,gravrate,gravno,gravelem,gravface)
global elem
%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);
ciclos=1;
%% inicializando dados para iteração Picard
%precondicionador ILU
[L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-9));

restarrt=7;
[p_old,fl1,rr1,it1,rv1]=gmres(M_old,RHS_old,restarrt,1e-9,1000,L,U);

[p_new,iter,res_hist,tabletol]=AndAcc(p_old,1e-8,kmap,parameter,w,s,...
    nflagface,fonte,gamma,nflagno,weightDMP,auxface,...
    wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,R0,tolpicard,...
    nitpicard,gravresult,gravrate,gravno,gravelem,gravface);
end