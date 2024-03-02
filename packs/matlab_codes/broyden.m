
function [x, i,ciclos,tolerancia]=broyden(M_old,RHS_old,p_old,tol,kmap,parameter,...
    w,s,nflagface,fonte,gamma,nflagno,...
    weightDMP,auxface,calnormface,wells,mobility,gravresult, gravrate,...
    gravno,gravelem,gravface,grav_elem_escalar,wg,pinterp)

global jacob
% inicializando
x = p_old;  
f = M_old*x-RHS_old;
R0 = norm(f); % residuo
if ~(size(x) == size(f))
    error('f must return a column vector of the same size as x0')
end

% montagem da matriz jacobiano na forma classica
if strcmp(jacob,'classic')
    
    J=aproxmjacobian(f,p_old,nflagface,nflagno,w,s,...
        parameter,weightDMP,kmap,fonte,0,mobility,0, 0, 0, 0, ...
        0,calnormface,wells,pinterp,gravresult,gravrate,gravno,gravelem,...
        gravface,grav_elem_escalar,wg);
else
    % montagem da matriz jacobiana usando paticoes
    % particoes
    [partic]=particoes(M_old); 
    % montagem da matriz jacobiano
    [J]=aproxmjacobianv2(f,p_old,nflagface,nflagno,w,s,parameter,weightDMP,...
        kmap,fonte,auxface,mobility,0, 0, 0, 0,0,calnormface,wells,partic,...
        M_old,pinterp,gravresult,gravrate,gravno,gravelem,gravface,...
        grav_elem_escalar,wg);
end
ciclos=1;
i=1;
er=1;
while  tol<er
    if any(isnan(J(:))) || rcond(full(J)) < 1e-15
        error('Singular jacobian at iteration %d\n',i)
    end
    
    dx=-J\f; %valor da jacobiana
    
    x  = x+dx; % perturba��o
    % calculo da interpola��o
    [pinterp_new]=pressureinterp(x,nflagface,nflagno,w,s,...
        parameter,weightDMP);
    % calculo da matriz
    [M_new,RHS_new]=globalmatrix(x,pinterp_new,gamma,nflagface,nflagno...
        ,parameter,kmap,fonte,w,s,weightDMP,auxface,...
        wells,mobility,0, 0, 0, 0, 0,calnormface,gravresult,...
        gravrate,gravno,gravelem,gravface,grav_elem_escalar,wg);
    
    
    f=M_new*x-RHS_new;
    
    J  = J + f*dx'/(dx'*dx);
    R = norm(f);
    
    if (R0 ~= 0.0)
        er = abs(R/R0)
    else
        er = 0;
    end    
    i=i+1;
end
tolerancia= er;
end
function [partic]=particoes(M_old) 

% particoes
    M_old1=full(M_old);
    % particoes
    n=size(M_old1,2); % coluna
    qq=n;
    Mnn=M_old1(:,1);
    M_old1(:,1)=0*M_old1(:,1);
    v(1,1)=1;
    nn=1;
    while qq>=n
        
        pp=1;k=1;kk=1;
        
        for j=1:n
            if max(M_old1(:,j))~=0
                % verifica produto interno
                MM=dot(Mnn,M_old1(:,j));
                % produto interno zero, entao a coluna 'j' faz parte da
                % particao
                if MM==0
                    % armazena as particoes
                    v(nn,k+1)=j;
                    % junta as colunas familia
                    Mnn(:,1)= Mnn(:,1)+ M_old1(:,j);
                    % faz zero o coluna que faz parte da familia
                    M_old1(:,j)=0*M_old1(:,j);
                    % incremeto
                    k=k+1;
                else
                    vv(1,kk)=j;
                    Mnew(:,pp)=M_old1(:,j);
                    pp=pp+1;
                    kk=kk+1;
                end
                
            end
            
        end
        
        nn=nn+1;
        % coloca o flag do primeiro coluna diferente se zero avaliado
        v(nn,1)=vv(1,1); 
        % faz zero a coluna que faz parte da particao
        M_old1(:,vv(1,1))=0*M_old1(:,vv(1,1)); 
        
        if max(max(M_old1))==0
            v(end,:)=[];
            break
        end
        Mnn=Mnew(:,1);
        clear Mnew  
    end
    partic=v;
end