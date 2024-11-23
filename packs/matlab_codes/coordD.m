function [xD,yD] = coordD(xJ,xI,coordv1,coordv2)

% {a1,a2} - point v1 vector
% {b1,b2} - (v2 - v1) vector
% {c1,c2} - collocation xI vector
% {d1,d2} - (xJ - xI) vector

a1 = coordv1(1);
a2 = coordv1(2);

b1 = coordv2(1) - a1;
b2 = coordv2(2) - a2;

c1 = xI(1);
c2 = xI(2);

d1 = xJ(1) - c1;
d2 = xJ(2) - c2;

m = -((a2*d1 - c2*d1 - a1*d2 + c1*d2)/(b2*d1 - b1*d2));

xD = a1 + m*b1;
yD = a2 + m*b2;

end