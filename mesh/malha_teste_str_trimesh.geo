boxdim = 2;
n = 8;
gridsize = boxdim/n;

Point(1) = {-1, -1, 0, gridsize};
Point(2) = {1, -1, 0, gridsize};
Point(3) = {1, 1, 0, gridsize};
Point(4) = {-1, 1, 0, gridsize};

Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line Loop(9) = {8, 5, 6, 7};

Plane Surface(10) = 9;

Transfinite Line{5, 6, 7, 8} = n + 1;
Transfinite Surface{10};
