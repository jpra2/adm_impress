function y = JV_APPROX1(v, R, p_old1,nflag,w,s,parameter,...
    kmap,nflagno,fonte,gamma,weightDMP,auxface,...
    wells,mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,M_new,RHS_new)


% calculate perturbation
dim = size(p_old1,1);

if norm(v, 2) > eps
    
    sum = 0;
    
    for i = 1 : dim
        
        sum = sum + sqrt(eps) * (1 + p_old1(i));
        
    end
    
    per = (1 / (dim * norm(v, 2))) * sum;
else
    
    sum = 0;
    
    for i = 1 : dim
        
        sum = sum + sqrt(eps) * (1 + p_old1(i));
        
    end
    
    per = sum / dim;
    
end
%xper1 = p_old1 - per * v; % perturbed vector
% Interpolação das pressões na arestas (faces)
%[pinterp_new]=pressureinterp(xper1,nflag,nflagno,w,s,auxflag,metodoP,parameter,weightDMP,mobility);
                                 
% Calculo da matriz global
%[M_new1,RHS_new1]=globalmatrix(xper1,pinterp_new,gamma,nflag,nflagno...
%    ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,...
%wells,mobility,Hesq, Kde, Kn, Kt, Ded,,calnormface);

% calculo do residuo
%Rper1= M_new1*xper1 - RHS_new1;

xper = p_old1 + per * v; % perturbed vector

% Interpolação das pressões na arestas (faces)
[pinterp_new]=pressureinterp(xper,nflag,nflagno,w,s,parameter,weightDMP,mobility);

% Calculo da matriz global
[M_new,RHS_new]=globalmatrix(xper,pinterp_new,gamma,nflag,nflagno...
    ,parameter,kmap,fonte,metodoP,w,s,weightDMP,auxface,wells,...
    mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);

% Calculo do residuo
Rper= M_new*xper - RHS_new;

% y = (Rper - Rper1)/(2*per);% approximation of jacobian action on krylov vector
  y = (Rper - R)/per;% 
end