function [p,step]=itersecant(M_old,RHS_old,nitnewton,tolnewton,interp,kmap,...
    parameter,metodoP,auxflag,w,s,nflag,valuemin,fonte,p_old,gamma,nflagno,benchmark,M_old1,RHS_old1,p_old1)
global elem
%% inicializando dados para iteração Picard
Rkk= M_old1*p_old1-RHS_old1;
Rk= M_old*p_old-RHS_old;
erropressureiter=1;
erroresidualiter=1;
step=0;
I=eye(size(elem,1));
[J]=aproxmjacobian(Rkk,Rk,p_old1,p_old);
while tolnewton<erroresidualiter | tolnewton<erropressureiter
    step=step+1;
    %% calculo da pressão
    pr=-inv(I-J)*(p_old1-Rkk);
    p_new=p_old+pr;
    %% plotagem no visit
    postprocessor(p_new,step)
    %% Interpolação das pressões na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflag,w,s,auxflag,metodoP,parameter,kmap,interp);
    
    %% Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(valuemin,p_new,pinterp_new,gamma,nflag,nflagno...
        ,parameter,kmap,fonte,metodoP,w,benchmark);
    %%    
    Rk=Rkk;
    Rkk= M_new*p_new - RHS_new+1e-23;
    p_old=p_old1;
    p_old1=p_new;
    
    [J]=aproxmjacobian(Rkk,Rk,p_old1,p_old);
    erroresidualiter=norm(Rkk)
end
p=p_old;
end