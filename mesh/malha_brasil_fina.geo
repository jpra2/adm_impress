Lx = 1.5;
Ly = 1;
R = 0.2;

d1 = 1.2;
d2 = 0.75;

x1 = (Lx - d1)/2;
x2 = x1 + d1/2;
x3 = x1 + d1;
x4 = x2;
x5 = Lx/2;
x6 = x5+R;
x7 =x5-R;

y1 = Ly/2;
y2 = (Ly-d2)/2 + d2;
y3 = y1;
y4 = y2 - d2;
y5 = Ly/2;
y6 = y5;
y7 = y5;

tr1 = 2;
tr2 = 2;
tr3 = 4;


gridsize = 0.5;
gridsize2 = 0.36;
gridsize3 = 0.7;

Point(1) = {0, 0, 0, gridsize};
Point(2) = {Lx, 0, 0, gridsize};
Point(3) = {Lx, Ly, 0, gridsize};
Point(4) = {0, Ly, 0, gridsize};

Point(5) = {x1, y1, 0, gridsize2};
Point(6) = {x2, y2, 0, gridsize2};
Point(7) = {x3, y3, 0, gridsize2};
Point(8) = {x4, y4, 0, gridsize2};

Point(9) = {x5, y5, 0, gridsize2};
Point(10) = {x6, y6, 0, gridsize2};
Point(11) = {x7, y7, 0, gridsize2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Circle(9) = {10, 9, 11};
Circle(10) = {11, 9, 10};


Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};
Curve Loop(3) = {9, 10};

Plane Surface(1) = {1, 2};
Plane Surface(2) = {2, 3};

//Transfinite Curve {1, 2, 3, 4} = tr1 + 1 Using Progression 1;
//Transfinite Curve {5, 6, 7, 8} = tr2 + 1 Using Progression 1;
//Transfinite Curve {9, 10} = tr3 + 1 Using Progression 1;

//Transfinite Surface {1};
//Transfinite Surface {2};

Recombine Surface {1, 2};

Physical Line("Inflow") = {4};
Physical Line("Outflow") = {2};
Physical Line("Walls") = {1, 3};
Physical Line("Boundary") = {1, 2, 3, 4, 9, 10};
Physical Line("External_Boundary") = {1, 2, 3, 4};
Physical Line("Internal_Boundary") = {9, 10};

Physical Surface("Region1") = {1};
Physical Surface("Region2") = {2};
