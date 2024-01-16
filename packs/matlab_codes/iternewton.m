function [p,step]=iternewton(M_old,RHS_old,nitnewton,tolnewton,interp,kmap,...
    parameter,metodoP,auxflag,w,s,nflag,valuemin,fonte,p_old,gamma,nflagno,benchmark,M_old1,RHS_old1,p_old1) 
    global elem
    %% calculo do residuo Inicial
    R0=M_old*p_old-RHS_old;
    %%
    R01=M_old1*p_old1-RHS_old1;
    %% inicializando dados para iteração Picard
    nitnewton=1;
    erropressureiter=1;
    erroresidualiter=1;
    step=0;
    [J]=aproxmjacobian(R01,R0,p_old1,p_old,nflag,w,s,metodoP,parameter,kmap,...
    interp,nflagno,benchmark,fonte,valuemin,auxflag,gamma);
    I=eye(size(elem,1));
    while tolnewton<erroresidualiter | tolnewton<erropressureiter  
        step=step+1;
        Lk=1;
        %% calculo da pressão
        Jinv= inv(J+Lk*I);
        p_new=p_old1-Jinv*R01;
        %% precondicionador
        %[L,U]=lu(M_old);
        
        %p_new=gmres(M_old,RHS_old,size(M_old,1),1e-14,500,L,U);
        %% plotagem no visit
        postprocessor(p_new,step)
        %% Interpolação das pressões na arestas (faces)
        [pinterp_new]=pressureinterp(p_new,nflag,w,s,auxflag,metodoP,parameter,kmap,interp);
        
        %% Calculo da matriz global
        [M_new,RHS_new]=globalmatrix(valuemin,p_new,pinterp_new,gamma,nflag,nflagno...
            ,parameter,kmap,fonte,metodoP,w,benchmark);
        %
        Rnew= M_new*p_new - RHS_new;
        [J]=aproxmjacobian(Rnew,R01,p_new,p_old1,nflag,w,s,metodoP,parameter,...
            kmap,interp,nflagno,benchmark,fonte,valuemin,auxflag,gamma);
        %% erro
        if (R0er ~= 0.0)
            erroresidualiter =norm(Rnew)
        else
            erroresidualiter = 0.0; %exact
        end
        %% iteração
        nitnewton = nitnewton + 1
        R01=Rnew;
        p_old1=p_new;
    end
    p=M_old\RHS_old;
    %[L,U]=lu(M_old);
    %p1=gmres(M_old,RHS_old,1000,1e-10,5000,L,U)
end