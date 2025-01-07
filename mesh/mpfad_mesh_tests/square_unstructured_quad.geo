Lx = 100;
Ly = 100;
cf = 1;
k1 = 3;
h1 = cf;
hx = Lx/k1;
hy = Ly/k1;
h = h1;
Point(1) = {0, 0, 0, h};
Point(2) = {Lx, 0, 0, h};
Point(3) = {Lx, Ly, 0, h};
Point(4) = {0, Ly, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Transfinite Curve {1, 3} = 3 + 1 Using Progression 1;
Transfinite Curve {2, 4} = 3 + 1 Using Progression 1;
Transfinite Surface {6};
Recombine Surface {6};

