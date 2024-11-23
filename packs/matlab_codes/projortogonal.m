function [Dorthog]=projortogonal(ielem,no1,no2)
global coord centelem

centerx=centelem(ielem,1);
centery=centelem(ielem,2);

point11=coord(no1,1);
point12=coord(no1,2);

point21=coord(no2,1);
point22=coord(no2,2);

pendent=(point22-point12)/(point21-point11+1e-16);
A=-pendent;
B=1;
C=-(point12-pendent*point11);

Dorthog=abs(A*centerx+B*centery+C)/sqrt(A^2+B^2);

end