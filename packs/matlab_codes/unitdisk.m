function [c,ceq] = unitdisk(t,Raux,media)

c = (t(1)-media(1,1))^2 + (t(2)-media(1,2))^2 - Raux^2; % Compute nonlinear inequalities at x.
ceq = [];% Compute nonlinear equalities at x.
end