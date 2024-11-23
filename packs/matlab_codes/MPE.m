% function [pressure]=MPE(A,b,x0)
% % MPE Minimal Polynomial Extrapolation
% % [Y,X,U,Gamma]=MPE(G,d,x0,n) 
% % performs n number of iteration steps
% % x(j+1)=G*x(j)+d starting with x0 and computes an extrapolated
% % vector y at each iteration, stored in the matrix Y. The basic
% % iteration steps and their differences are stored in X and U
% % respectively, and Gamma contains the coefficients of the
% % polynomial at each step.
% %##########################################################################
% % A=[ 0 -4 -8 -2
% % -4 -7 -7 -8
% % -9 -5 -4 -5
% % 0 -5 -9 -6 ];
% % n=max(size(A)); x=(1:n)'; b=A*x;
% % x0=x-5 ; % Starting vector
% G=eye(size(A))-A;
% 
% Gamma=zeros(size(A,1),size(A,1));
% n=length(b);
% Y=[]; U=[]; 
% %Gamma=[];
% xk=x0; 
% X=xk;
% p=0;
% for k=0:2
%     xk1=G*xk+b;
%     X=[X,xk1];
%     uk=xk1-xk;
%     U=[U,uk];
%     c=[-U(:,1:k)\uk;1];
%     y=X(:,1:k+1)*c/sum(c);
%     Gamma(:,k+1)=[c/sum(c);zeros(n-k-1,1)];
%     %Gamma(:,k+1)=c/sum(c);
%     Y=[Y,y];
%     xk=xk1;
%     p=p+1;
% end
% 
% R=b-A*Y;
% %norma=norm(R-U'*Gamma);
% u0=U(:,1);
% K=[u0 G*u0 G^2*u0 G^3*u0];
% 
% coeff=Gamma(:,p)'/Gamma(p,p);
% pressure=coeff(p:-1:1)-poly(G);
%######################################################################
