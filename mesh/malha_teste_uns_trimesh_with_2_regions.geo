Lx = 160;
Ly = Lx/2;

Lxc = Lx/2;
Lyc = Ly/2;
Lxx = Lx/6;
Lyy = Ly/6;

// mesh4 2 regions
// fine mesh 1 gridsize = 0.7
// coarse mesh 1 gridsize = 10

gridsize = 12;
gridsizee = gridsize/3;

Point(1) = {0, 0, 0, gridsize};
Point(2) = {Lx, 0, 0, gridsize};
Point(3) = {Lx, Ly, 0, gridsize};
Point(4) = {0, Ly, 0, gridsize};
Point(5) = {Lxc - Lxx, Lyc-Lyy, 0, gridsize};
Point(6) = {Lxc + Lxx, Lyc-Lyy, 0, gridsize};
Point(7) = {Lxc+Lxx, Lyc+Lyy, 0, gridsize};
Point(8) = {Lxc-Lxx, Lyc+Lyy, 0, gridsize};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};


Line Loop(1) = {1, 2, 3, 4, -5, -6, -7, -8};
Line Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line("Inflow") = {4};
Physical Line("Outflow") = {2};
Physical Line("Walls") = {1, 3};
Physical Line("Boundary") = {1, 2, 3, 4};

Physical Surface("Region1") = {1};
Physical Surface("Region2") = {2};
