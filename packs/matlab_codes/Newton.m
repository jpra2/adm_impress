function x=Newton(f,x0,tol,maxiter,nflagface,nflagno,w,s,auxflag,metodoP,...
    parameter,weightDMP,kmap,fonte,auxface,mobility,Hesq, Kde, Kn, Kt, ...
    Ded,calnormface,wells,benchmark);
% NEWTON solves a nonlinear system of equations
%x=Newton(f,x0,tol,maxiter,fp); solves the nonlinear system of
%equations f(x)=0 using Newtons methods, starting with the initial
%guess x0 up to a tolerance tol, doing at most maxit iterations. An
%analytical Jacobian can be given in the parameter fp
if nargin<4,
    maxit=100;
end;
if nargin<3,
    tol=1e6;
end;
x=x0; i=0; dx=ones(size(x0));
while norm(f)>tol & norm(dx)>tol & i<maxiter
    
        J=NumericalJacobian(f,x,nflagface,nflagno,w,s,auxflag,metodoP,...
    parameter,weightDMP,kmap,fonte,auxface,mobility,Hesq, Kde, Kn, Kt, ...
    Ded,calnormface,wells,benchmark);
    
    dx=-J\f; 
    x=x+dx;
    i=i+1;
end
if i>=maxiter
    error('Newton did not converge: maximum number of iterations exceeded');
    
end

function J=NumericalJacobian(f,x,nflagface,nflagno,w,s,auxflag,metodoP,...
    parameter,weightDMP,kmap,fonte,auxface,mobility,Hesq, Kde, Kn, Kt, ...
    Ded,calnormface,wells,benchmark);
% NUMERICALJACOBIAN computes a Jacobian numerically
%J=NumericalJacobian(f,x); computes numerically a Jacobian matrix
% for the function f at x
for i=1:length(x)
    xd=x; 
    h=sqrt(10^-8*(1+abs(xd(i)))); 
    xd(i)=xd(i)+h;
    % interpolação nos nós ou faces

    
    
end;
 %% Interpolação das pressões na arestas (faces)
    [pinterp_new]=pressureinterp(xd,nflagface,nflagno,w,s,auxflag,metodoP,...
        parameter,weightDMP);
    
    %% Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(xd,pinterp_new,0,nflagface,nflagno...
        ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,...
        mobility,Hesq, Kde, Kn, Kt, Ded,calnormface);
J=(M_new*xd-RHS_new-f)/h;



