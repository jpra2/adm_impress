load transition_matrix.dat
 % Convert text to a sparse matrix
 A = spconvert(transition_matrix);
 % find the matrix dimension
 n = max(size(A));
 % split into 3 vectors indicating all non-zero parts
 [i,j,s] = find(A);
 % rebuild the matrix with correct dimensions
 A = sparse(i,j,s,n,n);
 % create source of rank
 e = ones([n,1]) * 1.0/n;

 % --- Repeated MPE iteration --- %
 m = 30;
 p = 0.1;
 delta = Inf;
 epsilon = 0.00000001;
 i = 0;
 r = e;
 while delta > epsilon
 [v] = MPE(A, r, e, p, m);
 r_new = v;
 % r_new may be negative
 r_new = sign(sum(r_new)) * r_new;
 r_new = r_new *( 1 / norm(r_new, 1));
 
 % transform eigenvector back to default basis
 r = r_new;
 % one iteration of the power method
 y = (1-p)*(A*r);
 % norm(r,1) should be 1
 % Note: w > 0, because p > 0 (should trump any numerical errors)
 w = 1 - norm(y,1);
 % w > 0, so r_new > 0
 r_new = y + w*e;
 delta = norm(r - r_new, 1);
 i = i + 1;
 r = r_new;
 end
 csvwrite('MPEIteration_PageRank.dat', r);
 
 function [ v ] = MPE(A, r, e, p, maxit)
 n = max(size(e));
 R = zeros(n,maxit);
 % perform power iteration maxit times
 for i = 1:maxit
 y = (1-p)*(A*r);
 w = 1 - norm(y,1);
 r = y + w*e;
 R(:,i) = r;
 end
 U=R(:,2:end-1)-R(:,1:end-2);
 c= U \ -(R(:,end) - R(:,end-1));
 c(end+1,1)=1;
 v=(R(:,2:end)*c);
 v = bsxfun(@max, v, 0);
 v = v ./(norm(v, 1));
 end