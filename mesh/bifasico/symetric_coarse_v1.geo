Mesh.Algorithm = 6;

cl__1 = 0.01;
cl__2 = 0.02;
cl__3 = 0.5;
cl__4 = 0.5;

Lx = 1;
Ly = 1;

Point(1) = {0, 0, 0, cl__3};
Point(2) = {0.25, 0, 0, cl__3};
Point(3) = {0.5, 0, 0, cl__3};
Point(4) = {0.75, 0, 0, cl__4};
Point(5) = {1, 0, 0, cl__4};
Point(6) = {1, 0.25, 0, cl__4};
Point(7) = {1, 0.5, 0, cl__3};
Point(8) = {1, 0.75, 0, cl__3};
Point(9) = {1, 1, 0, cl__3};
Point(10) = {0.75, 1, 0, cl__3};
Point(11) = {0.5, 1, 0, cl__3};
Point(12) = {0.25, 1, 0, cl__3};
Point(13) = {0, 1, 0, cl__3};
Point(14) = {0, 0.75, 0, cl__3};
Point(15) = {0, 0.5, 0, cl__3};
Point(16) = {0, 0.25, 0, cl__3};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 1};
Line(18) = {3, 15};
Line(19) = {4, 14};
Line(20) = {6, 12};
Line(21) = {7, 11};

//Line(5) = {5, 6};
//Line(6) = {6, 7};
//Line(7) = {7, 8};
//Line(8) = {8, 5};

Curve Loop(1) = {1, 2, 18, 15, 16};
Curve Loop(2) = {3, 19, 14, -18};
Curve Loop(3) = {4, 5, 20, 12, 13, -19};
Curve Loop(4) = {6, 21, 11, -20};
Curve Loop(5) = {7, 8, 9, 10, -21};
//Curve Loop(6) = {-21, 7, 22, 10};
//Curve Loop(7) = {-22, 8, 9};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
//Plane Surface(6) = {6};
//Plane Surface(7) = {7};

Transfinite Curve {18} = 3 + 1 Using Progression 1;
Transfinite Curve {19} = 3 + 1 Using Progression 1;
Transfinite Curve {20} = 3 + 1 Using Progression 1;
Transfinite Curve {21} = 3 + 1 Using Progression 1;
//Transfinite Surface {1};
//Transfinite Surface {1};
//Transfinite Surface {3};
Recombine Surface{1};
Recombine Surface{2};
Recombine Surface{3};
Recombine Surface{4};
Recombine Surface{5};
//Recombine Surface{6};
//Recombine Surface{7};

//Line {5} In Surface {1};
//Line {6} In Surface {1};
//Line {7} In Surface {1};
//Line {8} In Surface {1};
//Line {9} In Surface {1};
//Line {10} In Surface {1};
//Line {11} In Surface {1};
//Line {12} In Surface {1};
//Physical Point(201) = {1};
//Physical Point(202) = {9};
//Physical Curve(201) = {8};
//Physical Curve(202) = {4};
//Physical Curve(203) = {1, 2, 3, 5, 6, 7};
//Physical Surface(1) = {1};
//Physical Surface(2) = {2};
//Physical Surface(3) = {3};
//Physical Surface(4) = {4};
//Physical Surface(5) = {5};
//Physical Surface(6) = {6};
//Physical Surface(7) = {7};

//Mesh 2;
