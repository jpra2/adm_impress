Lx = 160;
Ly = Lx/2;
gridsize = 0.7;
gridsizee = gridsize/3;

Point(1) = {0, 0, 0, gridsize};
Point(2) = {Lx, 0, 0, gridsize};
Point(3) = {Lx, Ly, 0, gridsize};
Point(4) = {0, Ly, 0, gridsize};

Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line Loop(9) = {5, 6, 7, 8};

Plane Surface(10) = {9};

Physical Line("Inflow") = {8};
Physical Line("Outflow") = {6};
Physical Line("Walls") = {5, 7};
Physical Line("Boundary") = {5, 6, 7, 8};

Physical Surface("Volume") = {10};