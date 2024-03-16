
function [p,step,ciclos,tolerancia] = JFNK1(tolnewton,kmap,...
    parameter,metodoP,w,s,nflag,fonte,gamma,...
    nflagno,benchmark,M_old1,RHS_old1,p_old1,R0,weightDMP,...
    auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface)
global elem
% inicializando
R= M_old1*p_old1-RHS_old1; %Residuo
x0=1e-1*ones(size(p_old1,1),1);
er=1;
step= 0; % iteration counter
ciclos=1;
while tolnewton<er
    step = step + 1 % update iteration counter
%     if step>3
    % O método de Newton- Krylov aqui implementado foi proposto no
    % artigo "Jacobian-free Newton–Krylov methods:a survey of approaches and applications
    % D.A. Knoll, D.E. Keyes e foi implementado no Matlab segui o site:
    % http://www.mathworks.com/matlabcentral/fileexchange/45170-jacobian-free-newton-krylov--jfnk--method
    % para sistema de equações dadas, nós adaptamos a nosso contexto.
    j_v_approx1 = @(v)JV_APPROX1(v, R,p_old1,nflag,w,s,metodoP,parameter,...
        kmap,nflagno,benchmark,fonte,gamma,weightDMP,auxface,...
        wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,M_old1,RHS_old1);
    
    % Calculo pelo método de iteração GMRES do Matlab
    [v,flag,relres,iter,resvec] = gmres(j_v_approx1, R,2,1e-5,size(elem,1)*0.25,[],[],x0); % solve for Krylov vector
    
    % Calculo pelo método de iteração BICGSTAB do Matlab
    %[v,flag,relres,iter,resvec]=bicgstab(j_v_approx1,R,1e-5,5000,[],[] ,p_old1);
    flag
    
    % Calculo da pressão na iteração k+1
    p_new = p_old1 - v;
    % garantido a positividade
    n=1;
    while min(p_new)<0
        p_new=p_old1-v*(1/2^n);
        n=2*n;
    end
%     else
%        p_new=M_old1\RHS_old1; 
%     end
    % Plotagem no visit
    %S=ones(size(p_new,1),1);
    %postprocessor(p_new,S,step)
    
    % Interpolação das pressões na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflag,nflagno,w,s,metodoP,...
        parameter,weightDMP,mobility);
    
    % Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflag,nflagno...
        ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,...
        auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
    
    % Calculo do residuo
    R= M_new*p_new - RHS_new;
    R1=norm(R);
    if (R0 ~= 0.0)
        er = abs(R1/R0)
        %er = R/norm(RHS_new)
    else
        er = 0.0; %exact
    end
    errorelativo(step)=er;
    errorelativo(step)=er;
    
    % Atualizando a pressão
    p_old1=p_new;
    M_old1=M_new;
    RHS_old1=RHS_new;
end
p=p_old1;
tolerancia=er;
end
