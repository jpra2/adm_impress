Mesh.Algorithm = 2;

Lx = 1.0;
Ly = 1.0;
h = 0.1;
hx = 0.1;
hy = 0.1;
kx=0.15;
cl1 = kx;
cl2 = kx;
cl3 = kx;

Point(1) = {0, 0, 0, cl2};
Point(2) = {0.5, 0, 0, cl1};
Point(3) = {Lx, 0, 0, cl1};
Point(4) = {Lx, 0.5, 0, cl1};
Point(5) = {Lx, Ly, 0, cl2};
Point(6) = {0.5, Ly, 0, cl1};
Point(7) = {0, Ly, 0, cl3};
Point(8) = {0, 0.5, 0, cl1};
Point(9) = {0.25, 0.25, 0, cl2};
Point(10) = {0.75, 0.25, 0, cl2};
Point(11) = {0.75, 0.75, 0, cl2};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {1, 9};
Line(10) = {9, 8};
Line(11) = {2, 9};
Line(12) = {9, 10};
Line(13) = {10, 7};
Line(14) = {3, 10};
Line(15) = {11, 10};
Line(16) = {11, 6};
Line(17) = {4, 11};
Line(18) = {5,11};


Line Loop(1) = {8, 9, 10};
Line Loop(2) = {1, 11, -9};
Line Loop(3) = {12, 13, 7, -10};
Line Loop(4) = {2, 14, -12, -11};
Line Loop(5) = {3, 17, 15, -14};
Line Loop(6) = {16, 6, -13, -15};
Line Loop(7) = {4, 18, -17};
Line Loop(8) = {-18, 5, -16};


Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};

//Transfinite Curve {1, 8, 4, 5} = 4 + 1 Using Progression -1;
//Transfinite Curve {9, 11} = 8 + 1 Using Progression 1;
Transfinite Curve {13} = 8 + 1 Using Progression 1;
//Transfinite Curve {2, 3, 6, 7} = 7 + 1 Using Progression 1;
//Transfinite Curve {4} = 3 + 1 Using Progression 1;
//Transfinite Surface {6};
Recombine Surface {1};
Recombine Surface {2};
Recombine Surface {3};
Recombine Surface {4};
Recombine Surface {5};
Recombine Surface {6};
Recombine Surface {7};
Recombine Surface {8};
Physical Surface(1) = {1, 2};
Physical Surface(2) = {3, 4};
Physical Surface(3) = {5, 6};
Physical Surface(4) = {7, 8};
Physical Point(201) = {1};
Physical Point(202) = {5};

