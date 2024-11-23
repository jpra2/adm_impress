function [x0,t tabletol] = extrapolate(x0, k, L, method,nflagface,...
    nflagno,w,s,parameter,weightDMP,kmap,fonte,auxface,mobility,...
    Hesq, Kde, Kn, Kt,Ded,calnormface,wells,tolpicard,R0)
% https://bitbucket.org/romanz/numericodes/src/default/
% obtido de esse site
N = numel(x0);
Q = zeros(N, k+1);
switch upper(method)
    case 'MPE', method = @mpe;
    case 'RRE', method = @rre;
    case 'N/A', method = '';
    case '', warning('Extrapolate:PlainIteration', 'No extrapolation');
    otherwise, error('Extrapolate:UnknownMethod', method);
end
% Perform L cycles of extrapolation method
%% Interpolação das pressões na arestas (faces)
[pinterp_new]=pressureinterp(x0,nflagface,nflagno,w,s,...
    parameter,weightDMP);

%% Calculo da matriz global
[M_new,RHS_new]=globalmatrix(x0,pinterp_new,0,nflagface,nflagno...
    ,parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
    mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
RR = norm(M_new*x0 - RHS_new);
if (R0 ~= 0.0)
    erro = abs(RR/R0)
    
else
    erro = 0.0; %exact
end
tabletol(1,1:2)=[0, erro];
for t = 1:L %
    % Compute (k+1) vectors, in addition to x0
    
%     if rcond(full(M_new))<1e-5
         [L,U] = ilu(M_new,struct('type','ilutp','droptol',1e-6));
         %min(10^-5,10^-1*erro)
         [Q(:,1),fl1,rr1,it1,rv1]=bicgstab(M_new,RHS_new,10^-12,1000,L,U,x0);
%     else
%        [Q(:,1),fl1,rr1,it1,rv1]=bicgstab(M_new,RHS_new,min(10^-5,10^-1*erro),1000);
%    end
        
    %[L,U] = ilu(M_new,struct('type','ilutp','droptol',1e-6));
    %[Q(:,1),fl1,rr1,it1,rv1]=bicgstab(M_new,RHS_new,min(10^-5,10^-1*erro),10,L,U,x0);
    %[Q(:,1),fl1,rr1,it1,rv1]=gmres(M_new,RHS_new,10,1e-9,1000,L,U);
    %residuals(t) = norm(x0 - Q(:, 1), 2);
    %erro=residuals(t);
    % gera a sequencia de vetores para extrapolar, "Picard"
    for i = 1:k
        mm=1;
        %% Interpolação das pressões na arestas (faces)
        [pinterp_new]=pressureinterp(Q(:,i),nflagface,nflagno,w,s,...
            parameter,weightDMP);
        
        %% Calculo da matriz global
        [M_new,RHS_new]=globalmatrix(Q(:,i),pinterp_new,0,nflagface,nflagno...
            ,parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
            mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
        
        RRR = norm(M_new*Q(:,i) - RHS_new);
        
        if (R0 ~= 0.0)
            erroRRE = abs(RRR/R0);
        else
            erroRRE = 0.0; %exact
        end
        
        if erroRRE<tolpicard
            x0=Q(:,i);
            mm=0;
            erro=erroRRE
            x0=x0-min(x0,0);
            break
            
        end
        Q(:,i+1)=M_new\RHS_new;
        % Compute (k+1) vectors, in addition to x0
        %[L,U] = ilu(M_new,struct('type','ilutp','droptol',1e-6));
        %[Q(:,i+1),fl1,rr1,it1,rv1]=bicgstab(M_new,RHS_new,min(10^-5,10^-1*erro),10,L,U,Q(:,i));
        %[Q(:,i+1),fl1,rr1,it1,rv1]=gmres(M_new,RHS_new,10,1e-9,1000,L,U);
    end
    
    % Compute differences (k+1)
    if mm~=0
        for i = k:-1:1
            Q(:, i+1) = Q(:, i+1) - Q(:, i);
        end
        Q(:, 1) = Q(:, 1) - x0;
        
        % Perform QR decomposition
        [Q, R] = MGS(sparse(Q));
        % Perform extrapolation
        [gamma,] = method(R, k); % s.t. x0 = X * gamma
        xi = 1 - cumsum(gamma(1:k)); % s.t. x0' = x0 + U * xi
        eta = Q(:, 1:k)*R(1:k, 1:k) * xi(:); % since U = Q * R
        x0 = x0 +  eta; % s.t. x0' = x0 + Q * R * xi
        
        x0=x0-min(x0,0); % se algum elemento de x0 é negativo, então sera zero, creterio de LIPNIKOV.
        
        %% Interpolação das pressões na arestas (faces)
        [pinterp_new]=pressureinterp(x0,nflagface,nflagno,w,s,...
            parameter,weightDMP);
        
        %% Calculo da matriz global
        [M_new,RHS_new]=globalmatrix(x0,pinterp_new,0,nflagface,nflagno...
            ,parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
            mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
        
        
        RR = norm(M_new*x0 - RHS_new);
        
        if (R0 ~= 0.0)
            erro = abs(RR/R0)
            
        else
            erro = 0.0; %exact
        end
        if k==12
            k=k
        else
            k=k+1
        end
    end
    tabletol(t+1,1:2)=[t, erro];
    
    if erro<tolpicard
        break
    end
    
    
    %Q = zeros(N, k+1);
end
end

% Reduced Rank Extrapolation (RRE)
function [gamma] = rre(R, k)
e = ones(k+1, 1);
d = backsub(R, backsub(R', e));
lambda = 1 / sum(d);
gamma = lambda * d;
%residual = sqrt(lambda);
end

% Minimal Polynomial Extrapolation (MPE)
function [gamma] = mpe(R, k)
c = backsub(R(1:k, 1:k), -R(1:k, k+1));
c = [c; 1];
gamma = c / sum(c);
%residual = abs(gamma(end)) * R(end, end);
end




































