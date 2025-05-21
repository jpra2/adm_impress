Mesh.Algorithm = 2;

Lx = 1.0;
Ly = 1.0;
h = 0.1;
hx = 0.1;
hy = 0.1;
cl1 = 0.015;

Point(1) = {0, 0, 0, cl1};
Point(2) = {0.5, 0, 0, cl1};
Point(3) = {Lx, 0, 0, cl1};
Point(4) = {Lx, 0.5, 0, cl1};
Point(5) = {Lx, Ly, 0, cl1};
Point(6) = {0.5, Ly, 0, cl1};
Point(7) = {0, Ly, 0, cl1};
Point(8) = {0, 0.5, 0, cl1};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {2, 8};
Line(10) = {3, 7};
Line(11) = {4, 6};

Line Loop(1) = {1, 9, 8};
Line Loop(2) = {2, 10, 7, -9};
Line Loop(3) = {3, 11, 6, -10};
Line Loop(4) = {4, 5, -11};


Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

//Transfinite Curve {1, 8, 4, 5} = 32 + 1 Using Progression 1;
//Transfinite Curve {9, 11} = 45 + 1 Using Progression 1;
//Transfinite Curve {10} = 85 + 1 Using Progression 1;
//Transfinite Curve {2, 3, 6, 7} = 32 + 1 Using Progression 1;
//Transfinite Curve {4} = 3 + 1 Using Progression 1;
//Transfinite Surface {6};
//Recombine Surface {1};
//Recombine Surface {2};
//Recombine Surface {3};
//Recombine Surface {4};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Point(201) = {1};
Physical Point(202) = {5};

