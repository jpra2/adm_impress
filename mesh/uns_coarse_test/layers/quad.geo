Lx = 1.0;
Ly = 1.0;
h = 0.1;
hx = 0.1;
hy = 0.1;
Point(1) = {0, 0, 0};
Point(2) = {Lx, 0, 0};
Point(3) = {Lx, Ly, 0};
Point(4) = {0, Ly, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Transfinite Curve {1} = 3 + 1 Using Progression 1;
Transfinite Curve {2} = 3 + 1 Using Progression 1;
Transfinite Curve {3} = 3 + 1 Using Progression 1;
Transfinite Curve {4} = 3 + 1 Using Progression 1;
Transfinite Surface {6};
Recombine Surface {6};

