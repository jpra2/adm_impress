function [x,iter,res_hist,tabletol] = AndAcc(x,auxtol,kmap,...
    parameter,w,s,nflagface,fonte,gamma,nflagno,...
    weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded,...
    calnormface,R0,tolpicard,nitpicard,gravresult,gravrate,...
            gravno,gravelem,gravface)
global elem
% This performs fixed-point iteration with or without Anderson
% acceleration for a given fixed-point map g and initial
% approximate solution x.
%%
% https://users.wpi.edu/~walker/Papers/anderson_accn_algs_imps.pdf
%
% Homer Walker (walker@wpi.edu), 10/14/2011.
% Set the method parameters.
%######################################################################
if nargin < 2
    error('AndAcc requires at least two arguments.');
end
%if nargin < 3
mMax = min(4,size(x,1)); % é o "m" no artigo e
%end
%if nargin < 4
itmax = nitpicard;
%end
%if nargin < 5
%atol = 1.e-6;
atol=auxtol;
%end
%if nargin < 6
%rtol = 1.e-6;
rtol=auxtol;
%end
%if nargin < 7
droptol = 1.e10;
%end
%if nargin < 8
beta = 1;
%end
%if nargin < 9
AAstart = 6; % iterações de picard
%end
% Initialize the storage arrays.
res_hist = []; % Storage of residual history.
DG = []; % Storage of g-value differences.
% Initialize printing.
if mMax == 0
    fprintf('\n No acceleration.');
elseif mMax > 0
    fprintf('\n Anderson acceleration, mMax = %d \n',mMax);
else
    error('AndAcc.m: mMax must be non-negative.');
end
fprintf('\n iter res_norm \n');
% Initialize the number of stored residuals.
mAA = 0; % é o "m_{k}" no artigo
% Top of the iteration loop.
%% Interpolação das pressões na arestas (faces)
% Apply g and compute the current residual norm.
[pinterp_new]=pressureinterp(x,nflagface,nflagno,w,s,...
    parameter,weightDMP);
%% Calculo da matriz global
[M_new,RHS_new]=globalmatrix(x,pinterp_new,gamma,nflagface,nflagno...
    ,parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
    mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);

RR = norm(M_new*x - RHS_new);

if (R0 ~= 0.0)
    erro = abs(RR/R0)
else
    erro = 0.0; %exact
end
tabletol(1,1:2)=[0, erro];
erroaux1=0;
erroaux2=erro;
for iter = 0:itmax
    
    if erro<tolpicard
        break
    else
        
        if abs(erroaux1-erroaux2)<1e-10
            % utilizo quando a sequencia de erros se mantem quase constantes
            gval=M_new\RHS_new;
        else
            % precondicionador
            [L,U] = ilu(M_new,struct('type','ilutp','droptol',1e-6));
            restarrt=7;
            [gval,]=gmres(M_new,RHS_new,restarrt,1e-9,1000,L,U);
        end
    end
    
    fval = gval - x;
    
    if mMax == 0 || iter < AAstart,
        % Without acceleration, update x <- g(x) to obtain the next
        % approximate solution.
        x = gval;
    else
        % Apply Anderson acceleration.
        % Update the df vector and the DG array.
        if iter > AAstart,
            df = fval-f_old;
            if mAA < mMax,
                DG = [DG gval-g_old];
            else
                DG = [DG(:,2:mAA) gval-g_old];
            end
            mAA = mAA + 1;
        end
        f_old = fval;
        g_old = gval;
        if mAA == 0
            % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
            x = gval;
        else
            % If mAA > 0, solve the least-squares problem and update the
            % solution.
            if mAA == 1
                % If mAA == 1, form the initial QR decomposition.
                R(1,1) = norm(df);
                Q = R(1,1)\df;
            else
                % If mAA > 1, update the QR decomposition.
                if mAA > mMax
                    % If the column dimension of Q is mMax, delete the first column and
                    % update the decomposition.
                    [Q,R] = qrdelete(Q,R,1);
                    mAA = mAA - 1;
                    % The following treats the qrdelete quirk described below.
                    if size(R,1) ~= size(R,2),
                        Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
                    end
                    % Explanation: If Q is not square, then qrdelete(Q,R,1) reduces the
                    % column dimension of Q by 1 and the column and row
                    % dimensions of R by 1. But if Q *is* square, then the
                    % column dimension of Q is not reduced and only the column
                    % dimension of R is reduced by one. This is to allow for
                    % MATLAB’s default "thick" QR decomposition, which always
                    % produces a square Q.
                end
                % Now update the QR decomposition to incorporate the new
                % column.
                for j = 1:mAA - 1
                    R(j,mAA) = Q(:,j)'*df;
                    df = df - R(j,mAA)*Q(:,j);
                end
                R(mAA,mAA) = norm(df);
                Q = [Q,R(mAA,mAA)\df];
                
            end
            if droptol > 0
                % Drop residuals to improve conditioning if necessary.
                condDF = cond(R);
                while condDF > droptol && mAA > 1
                    fprintf(' cond(D) = %e, reducing mAA to %d \n', condDF, mAA-1);
                    [Q,R] = qrdelete(Q,R,1);
                    DG = DG(:,2:mAA);
                    mAA = mAA - 1;
                    % The following treats the qrdelete quirk described above.
                    if size(R,1) ~= size(R,2),
                        Q = Q(:,1:mAA); R = R(1:mAA,:);
                    end
                    condDF = cond(R);
                end
            end
            % Solve the least-squares problem.
            gamma_AA = R\(Q'*fval);
            % Update the approximate solution.
            x = gval - DG*gamma_AA;
            % Apply damping if beta is a function handle or if beta > 0
            % (and beta ~= 1).
            if isa(beta,'function_handle'),
                x = x - (1-beta(iter))*(fval - Q*R*gamma_AA);
            else
                if beta > 0 && beta ~= 1
                    x = x - (1-beta)*(fval - Q*R*gamma_AA);
                end
            end
        end
    end
    x=x-min(x,0);
    % Apply g and compute the current residual norm.
    [pinterp_new]=pressureinterp(x,nflagface,nflagno,w,s,parameter,weightDMP);
    %% Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(x,pinterp_new,gamma,nflagface,nflagno...
        ,parameter,kmap,fonte,w,s,weightDMP,auxface,wells,...
        mobility,Hesq, Kde, Kn, Kt, Ded,calnormface,gravresult,gravrate,...
            gravno,gravelem,gravface);
    
    RR = norm(M_new*x - RHS_new);
    
    if (R0 ~= 0.0)
        erro = abs(RR/R0)
        
    else
        erro = 0.0; %exact
    end
    
    tabletol(iter+2,1:2)=[iter+1, erro];
    % guarda o erro atual e anterior com objetivo de comparar depois
    erroaux1=tabletol(iter+1,2);
    erroaux2=erro;
    %% plotagem no visit
    % S=ones(size(x,1),1);
    % postprocessor(x,S,iter)
    
end
