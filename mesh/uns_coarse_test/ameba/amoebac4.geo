// Gmsh project created on Mon Jun 22 20:44:21 2020
cl1 = 0.3;
cl2 = 0.009;
cl3 = 0.009;
cl4 = 0.009;
cl2 = 1.0;
cl3 = 1.0;
cl4 = 1.0;
// 0.015
r1 = 0.015;
r2 = 0.015;
r3 = 0.015;
x1 = 0.3;
y1 = 0.3;
x2 = 1.7;
y2 = 0.5;
x3 = 0.9;
y3 = 0.45;

// 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
Mesh.Algorithm =8;
//Mesh.CharacteristicLengthFactor = 5.23;
Point(1) = {0.18, 0.13, 0, cl1};
Point(2) = {0.04, 0.3, 0, cl1};
Point(3) = {0.1, 0.49, 0, cl1};
Point(4) = {0.34, 0.64, 0, cl1};
Point(5) = {0.53, 0.6, 0, cl1};
Point(6) = {0.8, 0.66, 0, cl1};
//Point(7) = {0.4, 0.6, 0, cl1};
Point(7) = {0.96, 0.77, 0, cl1};
Point(8) = {1.14, 0.94, 0, cl1};
Point(9) = {1.34, 0.93, 0, cl1};
Point(10) = {1.48, 0.83, 0, cl1};
Point(11) = {1.62, 0.69, 0, cl1};
Point(12) = {1.81, 0.64, 0, cl1};
Point(13) = {1.91, 0.48, 0, cl1};
Point(14) = {1.71, 0.3, 0, cl1};

Point(15) = {1.45, 0.27, 0, cl1};
Point(16) = {1.11, 0.24, 0, cl1};
Point(17) = {0.9, 0.08, 0, cl1};
Point(18) = {0.71, 0.02, 0, cl1};
Point(19) = {0.51, 0.03, 0, cl1};
Point(20) = {0.3, 0.07, 0, cl1};

Point(21) = {x1, y1, 0, cl2};
Point(22) = {x2, y2, 0, cl2};
// Point(23) = {x3, y3, 0, cl2};



Point(24) = {x1 + r1, y1, 0, cl2};
Point(25) = {x1, y1 + r1, 0, cl2};
Point(26) = {x1 -r1, y1, 0, cl2};
Point(27) = {x1, y1 -r1, 0, cl2};
Circle(2) = {24, 21,27};
Circle(3) = {27, 21,26};
Circle(4) = {26, 21,25};
Circle(5) = {25, 21,24};

Point(28) = {x2 + r2, y2, 0, cl3};
Point(29) = {x2, y2 + r2, 0, cl3};
Point(30) = {x2 -r2, y2, 0, cl3};
Point(31) = {x2, y2 -r2, 0, cl3};
Circle(6) = {28, 22,29};
Circle(7) = {29, 22,30};
Circle(8) = {30, 22,31};
Circle(9) = {31, 22,28};

// Point(32) = {x3 + r3, y3, 0, cl4};
// Point(33) = {x3, y3 + r3, 0, cl4};
// Point(34) = {x3 -r3, y3, 0, cl4};
// Point(35) = {x3, y3 -r3, 0, cl4};
// Circle(10) = {32, 23,33};
// Circle(11) = {33, 23,34};
// Circle(12) = {34, 23,35};
// Circle(13) = {35, 23,32};

//+
BSpline(1) = {1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20,1};
//BSpline(1) = {1,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1}; 



Line Loop(1) = {1};



// Line Loop(2) = {2,3,4,5};
// Line Loop(2) = {2,3};
// Line Loop(3) = {6,7,8,9};
// Line Loop(4) = {10,11,12,13};

// Plane Surface(1) = {1,2,3,4};
Plane Surface(1) = {1};
Transfinite Line {1} = 35 Using Progression 1;
// Transfinite Line {2,3,4,5,6,7,8,9,10,11,12,13} = 2 Using Progression 1;
// Transfinite Line {2,3,6,7,8,9,10,11,12,13} = 2 Using Progression 1;

//Plane Surface(2) = {2};
//Plane Surface(3) = {3};
//Plane Surface(4) = {4};
//Plane Surface(1) = {1,2,3,4};
//Plane Surface(2) = {2};
Plane Surface(1) = {1};
Recombine Surface{1};

Point(50) = {x1+0.02, y1+0.02, 0, cl1};
Point{50} In Surface {1};

// Point(51) = {x3+0.05, y3, 0, cl1};
// Point{51} In Surface {1};

Point(52) = {x2-0.05, y2, 0, cl1};
Point{52} In Surface {1};

// Point(53) = {x2-0.05, y2+0.02, 0, cl1};
// Point{53} In Surface {1};


// Point{21} In Surface {1};
// Point{22} In Surface {1};
// Point{23} In Surface {1};


//Circle{2} In Surface {1};
//Circle{3} In Surface {1};
//Circle{4} In Surface {1};
//Circle{5} In Surface {1};

//Recombine Surface{1};
//Recombine Surface{2};
//Recombine Surface{3};
//Recombine Surface{4};


// Line{2,3,4,5,6,7,8,9,10,11,12,13} In Surface {1};
//Physical Point(201) = {21,22,23};

Physical Line(201) = {1};
Physical Surface(1) = {1};
//Physical Surface(2) = {1};
