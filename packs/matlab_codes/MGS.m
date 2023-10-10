% Performs in-place QR using Modified Gram-Schmidt algorithm
%   [Q, R] = MGS(A), such that Q*R = A, Q'*Q = I
%   and R is upper-diagonal.
%   Note that we replace A by Q, thus saving memory.
function [Q, R] = MGS(Q)
%% Formato I-variante1.
[N, M] = size(Q);
R = zeros(M);
for k = 1:M
    for j = 1:k-1
        rjk = Q(:, j)' * Q(:, k);
        Q(:, k) = Q(:, k) - rjk * Q(:, j);
        R(j, k) = rjk;
    end
    R(k, k) = norm(Q(:, k));
    Q(:, k) = Q(:, k) / R(k, k);
end
%% FormatoII-variante1: Este algoritmo produzeos mesmos resultados que do Formato I
% [N, M] = size(A);
% Q=zeros(N,M);
% R=zeros(M,M);
% for k = 1:M
%     v=A(:,k);
%     for j = 1:k-1
%         R(j,k)=Q(:, j)' * A(:, k);
%         
%        v = v - R(j,k) * Q(:, j);
%       
%     end
%     R(k, k) = norm(v);
%     Q(:, k) = v / R(k, k);
% end


%% Formato III-variante2
% [N, M] = size(A);
% R=zeros(M,M);
% for k = 1:M
%     R(k, k) = norm(A(:,k));
%     Q(:, k) = A(:,k) / R(k, k);
%     
%     if k<M
%         R(k,k+1:M)=Q(:, k)' * A(:, k+1:M);
%         
%        A(:, k+1:M) = A(:, k+1:M) - Q(:,k) * R(k, k+1:M);
%       
%     end
% end

%% Formator IV-variante3: com reortogonalização
% [N, M] = size(Q);
% R=zeros(M);
% for k=1:M
%     tt = 0;
%     t = norm(Q(:,k));
%     reorth = 1;
%     while reorth,
%         % orth:
%         for i=1:k-1
%             s = Q(:,i)'*Q(:,k);
%             if tt==0, 
%                 R(i,k) = s; 
%             end;
%             Q(:,k) = Q(:,k)-s*Q(:,i);
%         end
%         tt = norm(Q(:,k));
%         reorth=0;
%         if tt<t/10, 
%             t=tt; 
%             reorth=1;
%         end
%     end
%     R(k,k) = tt;
%     Q(:,k) = Q(:,k)/R(k,k);
% end































