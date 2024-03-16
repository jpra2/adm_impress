boxdim = 1;
n = 16;
gridsize = boxdim/n;

Point(1) = {0, 0, 0, gridsize};
Point(2) = {1, 0, 0, gridsize};
Point(3) = {1, 1, 0, gridsize};
Point(4) = {0, 1, 0, gridsize};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Physical Line(101) = {1, 2, 3, 4};
Physical Surface(1) = {1};