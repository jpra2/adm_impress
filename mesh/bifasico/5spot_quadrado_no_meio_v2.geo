//SetFactory("OpenCASCADE");
Mesh.Algorithm = 6;
//Mesh.MeshSizeExtendFromBoundary = 0;
//Mesh.MeshSizeFromPoints = 0;
//Mesh.MeshSizeFromCurvature = 0;
//Mesh.Smoothing = 100;
//Mesh.SubdivisionAlgorithm = 5;
//Mesh.MeshSizeMin = 9.0;
Mesh.MeshSizeMax = 9.375;
Mesh.Color.Points = {255, 0, 0};

//cl__1 = 18.75;
//cl__2 = 18.75;

cl__1 = 9.375;
cl__2 = 9.375;

Lx = 300;
Ly = 300;

Point(1) = {0, 0, 0, cl__1};
Point(2) = {Lx, 0, 0, cl__1};
Point(3) = {Lx, Ly, 0, cl__1};
Point(4) = {0, Ly, 0, cl__1};
Point(10) = {4, 0, 0, cl__1};
Point(11) = {0, 4, 0, cl__1};
Point(12) = {Lx-4, Ly, 0, cl__1};
Point(13) = {Lx, Ly-4, 0, cl__1};


Point(5) = {93.75, 93.75, 0, cl__2};
Point(6) = {206.5, 93.75, 0, cl__2};
Point(7) = {206.5, 206.5, 0, cl__2};
Point(8) = {93.75, 206.5, 0, cl__2};

Line(1) = {1, 10};
Line(2) = {10, 2};
Line(3) = {2, 13};
Line(4) = {13, 3};
Line(5) = {3, 12};
Line(6) = {12, 4};
Line(7) = {4, 11};
Line(8) = {11, 1};

Line(9) = {5, 6};
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 5};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Curve Loop(2) = {9, 10, 11, 12};

Plane Surface(1) = {1,-2};
Plane Surface(2) = {2};
//Plane Surface(3) = {3};
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
//Physical Curve(201) = {8};
//Physical Curve(202) = {4};
//Physical Curve(203) = {1, 2, 3, 5, 6, 7};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
//Physical Surface(3) = {3};
