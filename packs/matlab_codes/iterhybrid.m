function [p,step,ciclos,tolerancia]=iterhybrid(M_old,RHS_old,tolnewton,kmap,...
    parameter,w,s,nflag,fonte,p_old,gamma,nflagno,p_old1,weightDMP,auxface)
global elem
% inicializando dados para iteração Picard
R0= M_old*p_old-RHS_old;
er=1;
step=0;
ciclos=1;
while tolnewton<er
    step=step+1;
    % calculo da pressão
    p_new=M_old\RHS_old;
    
    % plotagem no visit
    postprocessor(p_new,step)
    
    % Interpolação das pressões na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflag,w,s,parameter,weightDMP);
    
    % Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflag,nflagno...
        ,parameter,kmap,fonte,w,weightDMP,auxface);
    
    % Calculo do residuo
    R = M_new*p_new - RHS_new;
    % calculo do erro
    A=logical(norm(R0) ~= 0.0);
    B=logical(norm(R0) == 0.0);
    er=A*abs(norm(R)/norm(R0))+B*0
    errorelativo(step)=er;
    
    if er<1e-5
        
        p_old1=p_new;
        while tolnewton<er
            
            % O método de Newton- Krylov aqui implementado foi proposto no
            % artigo "Jacobian-free Newton–Krylov methods:a survey of approaches and applications
            % D.A. Knoll, D.E. Keyes e foi implementado no Matlab segui o site:
            % http://www.mathworks.com/matlabcentral/fileexchange/45170-jacobian-free-newton-krylov--jfnk--method
            % para sistema de equações dadas, nós adaptamos a nosso contexto.
            j_v_approx1 = @(v)JV_APPROX1(v, R,p_old1,nflag,w,s,parameter,...
        kmap,nflagno,fonte,auxflag,gamma,weightDMP,auxface);
            
            % Calculo pelo método de iteração GMRES do Matlab
            [v,flag,relres,iter,resvec] = gmres(j_v_approx1, R,2,1e-5,size(elem,1)*0.5); % solve for Krylov vector
            
            % Calculo pelo método de iteração BICGSTAB do Matlab
            %[v,flag]=bicgstabl(j_v_approx1,R,1e-6,200);
            flag
            
            % Calculo da pressão na iteração k+1
            p_new = p_old1 - v;
            
            % plotagem no visit
            postprocessor(p_new,step+1);
            
            % Interpolação das pressões na arestas (faces)
            [pinterp_new]=pressureinterp(p_new,nflag,w,s,parameter,weightDMP);
            
            % Calculo da matriz global
            [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflag,nflagno...
                ,parameter,kmap,fonte,w,weightDMP,auxface);
            
            % Calculo do residuo
            R= M_new*p_new - RHS_new;
            
            % Calculo do erro
            A=logical(norm(R0) ~= 0.0);
            B=logical(norm(R0) == 0.0);
            er=A*abs(norm(R)/norm(R0))+B*0
            
            errorelativo(step+1)=er;
            
            p_old1=p_new;
            step=step+1;
        end
    end
    % atualizar
    M_old=M_new;
    RHS_old=RHS_new;
end
p=p_old1;
tolerancia=er;
end