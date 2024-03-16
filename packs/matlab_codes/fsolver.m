function [p_new,tabletol,step,ciclos]=fsolver(M_old,RHS_old,nitpicard,tolpicard,kmap,...
    parameter,w,s,nflagface,fonte,p_old,gamma,nflagno,...
    weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface)
global elem
ciclos=1;
% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);
f = @(p_old)(M_old*p_old-RHS_old);
% inicializando dados para iteração Picard
step=0;
er=1;
contador=0;
erroaux1=0;
erroaux2=1e9;
while (tolpicard<er || tolpicard==er) && (step<nitpicard)
    %% atualiza iterações
    step=step+1;
    
    %options = optimset(@fsolve,'Algorithm','trust-region','JacobPattern');
    %options = optimset('Jacobian','on');
    [p_new,] = fzero(f,p_old);
    postprocessor(p_new,step)
    %         p_max=max(p_new)
    %         p_min=min(p_new)
    % Interpolação das pressões na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflagface,nflagno,w,s,parameter,weightDMP);
    
    % Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflagface,nflagno...
        ,parameter,kmap,fonte,w,s,weightDMP,auxface,...
        wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);
    f=@(p_new)(M_new*p_new-RHS_new);
    % calculo do residuo
    R = norm(M_new*p_new - RHS_new);
    % calculo do erro
    if (R0 ~= 0.0)
        er = abs(R/R0)
    else
        er = 0.0; %exact
    end
    
    M_old=M_new;
    RHS_old=RHS_new;
    tabletol(contador+1,1:2)=[contador, er];
    % guarda o erro atual e anterior com objetivo de comparar depois 
    if step==1
    erroaux1=0;
    erroaux2=er;
    else
        erroaux1=tabletol(step-1,2);
        erroaux2=er;
    end
    % contador local para guardar erros
    contador=contador+1;
    
    
    
end
end