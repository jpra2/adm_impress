% This routine generates an nm by nm mesh, commonly known as a
% Shestakov mesh, with 0 < r < rm and 0 < z < zm. The randomness
% is controlled by parameter "a", where 0 <= a <= 0.5. A value of 
% a = 0.5 gives a rectangular mesh. Boomerang zones are allowed 
% but bowties are not allowed.
%
% Matlab version by Jean C. Ragusa, Texas A&M University, April 1, 2013
%
% based on the original Fortran program by Alek Shestakov; that program 
% is found in: Nuclear Science & Engineering, Vol 105, pp.88-104 (1990),
% "Test Problems in Radiative Transfer Calculations", by A. Shestakov, 
% D. Kershaw, and G. Zimmerman
clear all; close all; clc;
% number of subdivisions of the original rectangle
nc =4;
% random parameter
a  = 0.25;
% rectangle dimensions
L = 1;
rm = L;
zm = L;
% allocate memory
nm = 2^nc + 1;
ind = nm-1;
r=zeros(nm,nm);
z=zeros(nm,nm);
% initialize 4 corners
for i = 0:1
    k = 1 + i*ind;
    for j = 0:1
        l = 1 + j*ind;
        r(k,l) = L*(1-i);
        z(k,l) = L*j;
    end
end
% fill in the rest of the points
for nl = 0:nc-1
    nn = 2^nl;
    inc = ind/nn;
    inch = inc/2;
    for k = 1:nn
        k1 = 1+(k-1)*inc;
        k2 = k1+inc;
        k3 = k1+inch;
        for l = 1:nn
            l1 = 1+(l-1)*inc;
            l2 = l1+inc;
            l3 = l1+inch;
            if (l==1)
                ar = a+rand(1,1)*(1-2*a);
                r(k3,l1) = ar*r(k1,l1)+(1-ar)*r(k2,l1);
                z(k3,l1) = ar*z(k1,l1)+(1-ar)*z(k2,l1);
            end
            ar = a+rand(1,1)*(1-2*a);
            r(k2,l3) = ar*r(k2,l1)+(1-ar)*r(k2,l2);
            z(k2,l3) = ar*z(k2,l1)+(1-ar)*z(k2,l2);
            ar = a+rand(1,1)*(1-2*a);
            r(k3,l2) = ar*r(k1,l2)+(1-ar)*r(k2,l2);
            z(k3,l2) = ar*z(k1,l2)+(1-ar)*z(k2,l2);
            if (k==1)
                ar = a+rand(1,1)*(1-2*a);
                r(k1,l3) = ar*r(k1,l1)+(1-ar)*r(k1,l2);
                z(k1,l3) = ar*z(k1,l1)+(1-ar)*z(k1,l2);
            end
            ar = a+rand(1,1)*(1-2*a);
            br = a+rand(1,1)*(1-2*a);
            r1 = r(k1,l1);
            r2 = r(k2,l1);
            r3 = r(k2,l2);
            r4 = r(k1,l2);
            z1 = z(k1,l1);
            z2 = z(k2,l1);
            z3 = z(k2,l2);
            z4 = z(k1,l2);
            % check for boomerang zones
            det2 = (r2-r1)*(z3-z2)-(r3-r2)*(z2-z1);
            det3 = (r3-r2)*(z4-z3)-(r4-r3)*(z3-z2);
            det4 = (r4-r3)*(z1-z4)-(r1-r4)*(z4-z3);
            det1 = (r1-r4)*(z2-z1)-(r2-r1)*(z1-z4);
            if (det2>0)
                d = (r4-r3)*(z2-z1)-(r2-r1)*(z4-z3);
                r3p = ((r2-r1)*(r4*z3-r3*z4)-(r4-r3)*(r2*z1-r1*z2))/d;
                z3p = ((z2-z1)*(r4*z3-r3*z4)-(z4-z3)*(r2*z1-r1*z2))/d;
                d = (r4-r1)*(z2-z3)-(r2-r3)*(z4-z1);
                r1p = ((r2-r3)*(r4*z1-r1*z4)-(r4-r1)*(r2*z3-r3*z2))/d;
                z1p = ((z2-z3)*(r4*z1-r1*z4)-(z4-z1)*(r2*z3-r3*z2))/d;
                r3 = r3p;
                z3 = z3p;
                r1 = r1p;
                z1 = z1p;
            elseif (det3>0)
                d = (r1-r4)*(z3-z2)-(r3-r2)*(z1-z4);
                r4p = ((r3-r2)*(r1*z4-r4*z1)-(r1-r4)*(r3*z2-r2*z3))/d;
                z4p = ((z3-z2)*(r1*z4-r4*z1)-(z1-z4)*(r3*z2-r2*z3))/d;
                d = (r1-r2)*(z3-z4)-(r3-r4)*(z1-z2);
                r2p = ((r3-r4)*(r1*z2-r2*z1)-(r1-r2)*(r3*z4-r4*z3))/d;
                z2p = ((z3-z4)*(r1*z2-r2*z1)-(z1-z2)*(r3*z4-r4*z3))/d;
                r4 = r4p;
                z4 = z4p;
                r2 = r2p;
                z2 = z2p;
            elseif (det4>0)
                d = (r2-r1)*(z4-z3)-(r4-r3)*(z2-z1);
                r1p = ((r4-r3)*(r2*z1-r1*z2)-(r2-r1)*(r4*z3-r3*z4))/d;
                z1p = ((z4-z3)*(r2*z1-r1*z2)-(z2-z1)*(r4*z3-r3*z4))/d;
                d = (r2-r3)*(z4-z1)-(r4-r1)*(z2-z3);
                r3p = ((r4-r1)*(r2*z3-r3*z2)-(r2-r3)*(r4*z1-r1*z4))/d;
                z3p = ((z4-z1)*(r2*z3-r3*z2)-(z2-z3)*(r4*z1-r1*z4))/d;
                r1 = r1p;
                z1 = z1p;
                r3 = r3p;
                z3 = z3p;
            elseif (det1>0)
                d = (r3-r2)*(z1-z4)-(r1-r4)*(z3-z2);
                r2p = ((r1-r4)*(r3*z2-r2*z3)-(r3-r2)*(r1*z4-r4*z1))/d;
                z2p = ((z1-z4)*(r3*z2-r2*z3)-(z3-z2)*(r1*z4-r4*z1))/d;
                d = (r3-r4)*(z1-z2)-(r1-r2)*(z3-z4);
                r4p = ((r1-r2)*(r3*z4-r4*z3)-(r3-r4)*(r1*z2-r2*z1))/d;
                z4p = ((z1-z2)*(r3*z4-r4*z3)-(z3-z4)*(r1*z2-r2*z1))/d;
                r2 = r2p;
                z2 = z2p;
                r4 = r4p;
                z4 = z4p;
            end
            r(k3,l3) = ar*br*r1 + ar*(1-br)*r4 + (1-ar)*br*r2 + (1-ar)*(1-br)*r3;
            z(k3,l3) = ar*br*z1 + ar*(1-br)*z4 + (1-ar)*br*z2 + (1-ar)*(1-br)*z3;
        end
    end
end
% plotting in matlab
figure(1)
surf(r,z,ones(nm,nm))
view(0,-90)
