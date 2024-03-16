% BACKSUB Performs backsubstitution to solve Ax = b.
%   x = BACKSUB(A, b), where A must be 
%   lower/upper-diagonal matrix.
function x = backsub(A, b)
% 
% N = numel(b);
% x = zeros(N, 1);
% 
% [J, I] = meshgrid(1:N, 1:N);
% K = J > I;
% is_lower = ~any(any(A(K)));
% is_upper = ~any(any(A(K')));
% if ~is_lower && ~is_upper
%     error('Invalid matrix!');
% end
%     
% I = 1:N; % is_lower
% if is_upper
%     I = fliplr(I);
% end
x=A\b;

% for i = I
%     x(i) = b(i) / A(i, i);
% 
%     b = b - x(i) * A(:, i);
% end
