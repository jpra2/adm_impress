Mesh.Algorithm = 2;

//cl__1 = 3.5;
//cl__2 = 3.5;

cl__1 = 1;
cl__2 = 1;

Lx = 100;
Ly = 100;

Point(1) = {0, 0, 0, cl__1};
Point(2) = {Lx, 0, 0, cl__1};
Point(3) = {Lx, Ly, 0, cl__1};
Point(4) = {0, Ly, 0, cl__1};

Point(5) = {40, 80, 0, cl__2};
Point(6) = {44, 80, 0, cl__2};
Point(7) = {44, 0, 0, cl__2};
Point(8) = {40, 0, 0, cl__2};

Point(9) = {56, Ly, 0, cl__2};
Point(10) = {60, Ly, 0, cl__2};
Point(11) = {60, 20, 0, cl__2};
Point(12) = {56, 20, 0, cl__2};

Line(1) = {1, 8};
Line(2) = {8, 7};
Line(3) = {7, 2};
Line(4) = {2, 3};
Line(5) = {3, 10};
Line(6) = {10, 9};
Line(7) = {9, 4};
Line(8) = {4, 1};

Line(9) = {7, 6};
Line(10) = {6, 5};
Line(11) = {5, 8};

Line(12) = {10, 11};
Line(13) = {11, 12};
Line(14) = {12, 9};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
//Curve Loop(2) = {11, 10, 9, 2};
Curve Loop(2) = {2, 11, 10, 9};
Curve Loop(3) = {-6, 12, 13, 14};

Plane Surface(1) = {1, 3, 2};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
//Transfinite Curve {4, 8} = 40 + 1 Using Progression 1;
//Transfinite Curve {11, 9, 14, 12} = 10 + 1 Using Progression 1;
//Transfinite Surface {1};
//Transfinite Surface {2};
//Transfinite Surface {3};

//Line {5} In Surface {1};
//Line {6} In Surface {1};
//Line {7} In Surface {1};
//Line {8} In Surface {1};
//Line {9} In Surface {1};
//Line {10} In Surface {1};
//Line {11} In Surface {1};
//Line {12} In Surface {1};
//Physical Point(201) = {1, 2, 3, 4};
Physical Curve(201) = {8};
Physical Curve(202) = {4};
Physical Curve(203) = {1, 2, 3, 5, 6, 7};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
