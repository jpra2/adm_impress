Mesh.Algorithm = 2;

cl__1 = 0.3;
cl__2 = 0.6;
cl__3 = 0.01;
cl__4 = 0.08;

Lx = 1;
Ly = 1;

Point(1) = {0, 0, 0, cl__1};
Point(2) = {Lx, 0, 0, cl__2};
Point(3) = {Lx, Ly, 0, cl__1};
Point(4) = {0, Ly, 0, cl__2};

//Point(5) = {0.25, 0.25, 0, cl__3};
//Point(6) = {0.75, 0.25, 0, cl__2};
//Point(7) = {0.75, 0.75, 0, cl__3};
//Point(8) = {0.25, 0.75, 0, cl__2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//Line(5) = {5, 6};
//Line(6) = {6, 7};
//Line(7) = {7, 8};
//Line(8) = {8, 5};

Curve Loop(1) = {1, 2, 3, 4};
//Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1};
//Plane Surface(2) = {2};

//Transfinite Curve {1, 2} = 30 + 1 Using Progression 1;
//Transfinite Curve {3, 4} = 30 + 1 Using Progression 1;
//Transfinite Curve {5, 7} = 35 + 1 Using Progression 1;
//Transfinite Curve {6, 8} = 35 + 1 Using Progression 1;
//Transfinite Surface {1};
//Transfinite Surface {1};
//Transfinite Surface {3};
Recombine Surface{1};
//Recombine Surface{2};

//Line {5} In Surface {1};
//Line {6} In Surface {1};
//Line {7} In Surface {1};
//Line {8} In Surface {1};
//Line {9} In Surface {1};
//Line {10} In Surface {1};
//Line {11} In Surface {1};
//Line {12} In Surface {1};
//Physical Point(201) = {1, 2, 3, 4};
//Physical Curve(201) = {8};
//Physical Curve(202) = {4};
//Physical Curve(203) = {1, 2, 3, 5, 6, 7};
//Physical Surface(1) = {1};
//Physical Surface(2) = {2};
//Physical Surface(3) = {3};

//Mesh 2;
