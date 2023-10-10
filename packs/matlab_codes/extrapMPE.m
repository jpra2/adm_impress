
function x= extrapMPE(XX)

% método MPE
%########################################################
%https://en.wikipedia.org/wiki/Minimum_polynomial_extrapolation
U=XX(:,2:end-1)-XX(:,1:end-2);

c=-pinv(U)*(XX(:,end)-XX(:,end-1));

c(end+1,1)=1;
x=(XX(:,2:end)*c)/sum(c);
%######################################################

end

% mr    Modified Gram-Schmidt. This is a more stable way to compute a
%       QR factorization
function [Q, R] = mgs(A);
[m,n] = size(A);
% we assume that m>= n.
V = A;
R = zeros(n,n);
for i=1:n,
   R(i,i) = norm(V(:,i));
   V(:,i) = V(:,i)/R(i,i);
   if (i < n)
      for j = i+1:n
         R(i,j) = V(:,i)' * V(:,j);
         V(:,j) = V(:,j) - R(i,j) * V(:,i);
      end
   end
end
Q = V;
end

