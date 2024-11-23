function [xM,xD] = auxiliar_points(inedge,xt,coord)

xM = zeros(size(inedge,1),3);
xD = zeros(size(inedge,1),3);

for i = 1:size(inedge,1)
    
    v1 = inedge(i,1);
    v2 = inedge(i,2);
    
    Iele = inedge(i,3);
    Jele = inedge(i,4);
    
    %midle-point
    xM(i,:) = (coord(v1,:)+coord(v2,:))/2;
    
    % {a1,a2} - point v1 vector
    % {b1,b2} - (v2 - v1) vector
    % {c1,c2} - collocation xI vector
    % {d1,d2} - (xJ - xI) vector
    
    a1 = coord(v1,1);
    a2 = coord(v1,2);
    
    b1 = coord(v2,1) - a1;
    b2 = coord(v2,2) - a2;
    
    c1 = xt(Iele,1);
    c2 = xt(Iele,2);
    
    d1 = xt(Jele,1) - c1;
    d2 = xt(Jele,2) - c2;
    
    m = -((a2*d1 - c2*d1 - a1*d2 + c1*d2)/(b2*d1 - b1*d2));
    
    xD(i,:) = [a1 a2 0] + m*[b1 b2 0];
    
end
end