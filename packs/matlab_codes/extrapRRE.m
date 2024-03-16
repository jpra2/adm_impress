% extrap   Reduced rank exrapolation (RRE)
%
% Usage:
%
%  called from smacof_rre
%
% TOSCA = Toolbox for Surface Comparison and Analysis
% Web: http://tosca.cs.technion.ac.il
% Version: 1.0
%
% (C) Copyright Michael Bronstein, 2008
% All rights reserved.
%
% License:
%
% ANY ACADEMIC USE OF THIS CODE MUST CITE THE ABOVE REFERENCES.
% ANY COMMERCIAL USE PROHIBITED. PLEASE CONTACT THE AUTHORS FOR
% LICENSING TERMS. PROTECTED BY INTERNATIONAL INTELLECTUAL PROPERTY
% LAWS. PATENTS PENDING.

function x= extrapRRE(XX)
K = size(XX,2)-1;
dX  = XX(:,2:end) - XX(:,1:end-1);
%######################################################
%K = size(XX,2)-2;
%dX  = XX(:,3:end)-2*XX(:,2:end-1)+XX(:,1:end-2);
%#####################################################
rank(dX)
[Q,R] = mgs(dX);
y = R'\ones(K,1);
gamma = R\y;
gamma = gamma(:)/sum(gamma);

%x=sum(repmat(gamma',[size(XX,1) 1]).*XX(:,1:K),2);
%#######################################################
% alternativa RRE
 for k=1:length(gamma)
     if k==1
         ksi(1,1)=1-gamma(1,1);
     else
        ksi(k,1)=ksi(k-1,1)-gamma(k,1);
        
     end
 end
 x=XX(:,1)+Q*R*ksi; % xx é igual a x de arriba
end


% mr    Modified Gram-Schmidt. This is a more stable way to compute a
%       QR factorization
function [Q, R] = mgs(A)
% [m,n] = size(A);
% % we assume that m>= n.
% V = A;
% R = zeros(n,n);
% for i=1:n
%    R(i,i) = norm(V(:,i));
%    V(:,i) = V(:,i)/R(i,i);
%    if (i < n)
%       for j = i+1:n
%          R(i,j) = V(:,i)' * V(:,j);
%          V(:,j) = V(:,j) - R(i,j) * V(:,i);
%       end
%    end
% end
% Q = V;
%###########################################
% implemtacao feito pelo pesoal do MIT
% https://ocw.mit.edu/courses/mathematics/18-06-linear-algebra-spring-2010/related-resources/MIT18_06S10_gramschmidtmat.pdf
[m,n] = size(A);
Q=zeros(m,n); R=zeros(n,n);
for j=1:n
    v=A(:,j);
    for i=1:j-1
        R(i,j)=Q(:,i)'*A(:,j);
        v=v-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end
end


